"""Conformational Benchmark and Scientific Validation Suite.

Provides standard test molecular systems with theoretical / experimental literature references:
  - n-Butano: Anti vs Gauche torsional equilibrium (steric repulsion & 1,4 interactions)
  - 1,2-Etandiolo: Gauche vs Anti equilibrium (intramolecular hydrogen bonding)
  - Cicloesano: Chair vs Twist-Boat conformational equilibrium (ring strain)
"""

import time
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
from rdkit import Chem

from energyinator.config import CalculationConfig
from energyinator.pipeline import run_pipeline
from energyinator.models import CalculationResult
from energyinator.geometry import parse_xyz, calc_dihedral

logger = logging.getLogger("energyinator")

# Literature-curated benchmark reference data with full physical and bibliographic annotations
BENCHMARK_SYSTEMS = {
    "butane": {
        "id": "butane",
        "name": "n-Butano",
        "smiles": "CCCC",
        "process_description": "Equilibrio torsionale attorno al legame C2-C3: Anti (esteso) -> Gauche (ripiegato)",
        "conformation_a": "Anti (diedro C-C-C-C ~ 180°, minimo assoluto a minore energia)",
        "conformation_b": "Gauche (diedro C-C-C-C ~ ±60°, conformero secondario)",
        "reference_delta_e_kcalmol": 0.86,  # Gauche is ~0.86 kcal/mol higher than Anti in gas phase
        "reference_level": "CCSD(T)/CBS // Exp. Spettroscopico",
        "citation": "Dixon et al., J. Phys. Chem. A 2000, 104, 3438–3443; Allinger et al., JACS 1990, 112, 8293",
        "scientific_significance": "Valuta la capacità del metodo di riprodurre l'ingombro sterico e le interazioni repulsive 1,4 tra metili.",
        "tolerance_kcalmol": 0.40,
        "smarts_dihedral": "CCCC",
        "expected_order": "a_lowest"  # Conformation A is lowest, Delta E = E(B) - E(A)
    },
    "ethylene_glycol": {
        "id": "ethylene_glycol",
        "name": "1,2-Etandiolo (Glicole)",
        "smiles": "OCCO",
        "process_description": "Equilibrio torsionale O-C-C-O: Gauche (ponte H chelato) -> Anti (aperto non chelato)",
        "conformation_a": "Gauche (diedro O-C-C-O ~ 60°, stabilizzato da legame a idrogeno intramolecolare OH···O)",
        "conformation_b": "Anti (diedro O-C-C-O ~ 180°, catena aperta priva di ponte a idrogeno)",
        "reference_delta_e_kcalmol": 1.95,  # Anti is ~1.95 kcal/mol higher than Gauche in gas phase
        "reference_level": "CCSD(T)/aug-cc-pVTZ // Spettroscopia MW",
        "citation": "Howard & Jardon, J. Phys. Chem. A 1998, 102, 6084–6089; Freeman et al., J. Mol. Struct. 2001",
        "scientific_significance": "Verifica se il metodo quantistico / campo di forze descrive accuratamente i ponti a idrogeno intramolecolari deboli.",
        "tolerance_kcalmol": 0.60,
        "smarts_dihedral": "OCCO",
        "expected_order": "a_lowest"  # Conformation A (Gauche) is lowest, Delta E = E(Anti) - E(Gauche)
    },
    "cyclohexane": {
        "id": "cyclohexane",
        "name": "Cicloesano",
        "smiles": "C1CCCCC1",
        "process_description": "Inversione e flessione d'anello: Sedia (Chair) -> Barca-Ritorta (Twist-Boat)",
        "conformation_a": "Sedia (Chair, geometria rigida priva di eclissamenti torsionali, minimo globale)",
        "conformation_b": "Barca-Ritorta / Twist-Boat (minimo secondario flessibile ad alta energia)",
        "reference_delta_e_kcalmol": 5.50,  # Twist-boat is ~5.50 kcal/mol higher than Chair
        "reference_level": "Sperimentale Gas-Phase / CBS-Q",
        "citation": "Squillacote et al., J. Am. Chem. Soc. 1975, 97, 3299–3300; NIST Chemistry WebBook",
        "scientific_significance": "Quantifica l'accuratezza nella stima della tensione d'anello (ring strain) e delle interazioni 1,3-diassiali.",
        "tolerance_kcalmol": 1.20,
        "smarts_dihedral": None,
        "expected_order": "a_lowest"
    }
}


def _analyze_conformations_by_geometry(
    system_key: str,
    res: CalculationResult,
    sys_info: Dict[str, Any]
) -> Tuple[float, str, str]:
    """Identifies conformational states via geometric torsional analysis.
    
    Returns:
        (computed_delta_e, status_note, geometric_details)
    """
    conformers = res.conformers
    if not conformers:
        return 0.0, "Nessun conformero generato", "Dati non disponibili"

    if len(conformers) == 1:
        return 0.0, "Solo 1 conformero isolato", "Un solo conformero rilevato entro la finestra energetica"

    smarts = sys_info.get("smarts_dihedral")
    if not smarts:
        # Generic approach: difference between lowest energy and first higher conformer
        computed_delta = conformers[1].delta_dft_kcalmol or 0.0
        return (
            computed_delta,
            "Confronto tra Minimo Globale (Conf 0) e primo conformero superiore (Conf 1)",
            f"Conf 0 vs Conf 1: Delta E = {computed_delta:.2f} kcal/mol"
        )

    # Substructure search for dihedral definition
    mol = Chem.MolFromSmiles(sys_info["smiles"])
    if not mol:
        computed_delta = conformers[1].delta_dft_kcalmol or 0.0
        return computed_delta, "Analisi generica", f"Delta E = {computed_delta:.2f} kcal/mol"

    query = Chem.MolFromSmarts(smarts)
    match = mol.GetSubstructMatch(query)
    if len(match) < 4:
        computed_delta = conformers[1].delta_dft_kcalmol or 0.0
        return computed_delta, "Analisi generica", f"Delta E = {computed_delta:.2f} kcal/mol"

    a0, a1, a2, a3 = match[:4]

    # Measure dihedral for each conformer
    dihedrals = []
    for c in conformers:
        if not c.xyz:
            dihedrals.append(None)
            continue
        try:
            _, coords = parse_xyz(c.xyz)
            d_val = calc_dihedral(coords[a0], coords[a1], coords[a2], coords[a3])
            dihedrals.append(d_val)
        except Exception:
            dihedrals.append(None)

    if system_key == "butane":
        # Anti: |dihedral| >= 120°; Gauche: 30° <= |dihedral| < 120°
        anti_confs = [(c, d) for c, d in zip(conformers, dihedrals) if d is not None and abs(d) >= 120.0]
        gauche_confs = [(c, d) for c, d in zip(conformers, dihedrals) if d is not None and 30.0 <= abs(d) < 120.0]

        if anti_confs and gauche_confs:
            best_anti, d_anti = anti_confs[0]
            best_gauche, d_gauche = gauche_confs[0]
            e_anti = best_anti.dft_energy_kcalmol or (best_anti.dft_energy_hartree * 627.509)
            e_gauche = best_gauche.dft_energy_kcalmol or (best_gauche.dft_energy_hartree * 627.509)
            computed_delta = round(e_gauche - e_anti, 2)
            geo_desc = f"Anti (diedro {d_anti:.1f}°) vs Gauche (diedro {d_gauche:.1f}°)"
            return computed_delta, "Equilibrio Anti -> Gauche identificato geometricamente", geo_desc

    elif system_key == "ethylene_glycol":
        # Gauche (H-bonded): |dihedral| < 120°; Anti (open): |dihedral| >= 120°
        gauche_confs = [(c, d) for c, d in zip(conformers, dihedrals) if d is not None and abs(d) < 120.0]
        anti_confs = [(c, d) for c, d in zip(conformers, dihedrals) if d is not None and abs(d) >= 120.0]

        if gauche_confs and anti_confs:
            best_gauche, d_gauche = gauche_confs[0]
            best_anti, d_anti = anti_confs[0]
            e_gauche = best_gauche.dft_energy_kcalmol or (best_gauche.dft_energy_hartree * 627.509)
            e_anti = best_anti.dft_energy_kcalmol or (best_anti.dft_energy_hartree * 627.509)
            computed_delta = round(e_anti - e_gauche, 2)
            geo_desc = f"Gauche chelato (diedro {d_gauche:.1f}°) vs Anti aperto (diedro {d_anti:.1f}°)"
            return computed_delta, "Rottura ponte H (Gauche -> Anti) identificata", geo_desc

    # Fallback to delta between conformer 1 and conformer 0
    computed_delta = conformers[1].delta_dft_kcalmol or 0.0
    d0 = dihedrals[0] if dihedrals and dihedrals[0] is not None else 0.0
    d1 = dihedrals[1] if len(dihedrals) > 1 and dihedrals[1] is not None else 0.0
    geo_desc = f"Minimo 0 ({d0:.1f}°) vs Conformero 1 ({d1:.1f}°)"
    return computed_delta, "Confronto minimi consecutivi", geo_desc


def run_benchmark(
    systems: Optional[List[str]] = None,
    method: str = "b3lyp",
    basis: str = "6-31G(d)",
    force_field: str = "mmff94",
    solvent: Optional[str] = None,
    sp_only: bool = True
) -> Dict[str, Any]:
    """Executes the standard conformational benchmark suite and computes statistical accuracy metrics."""
    target_keys = systems or list(BENCHMARK_SYSTEMS.keys())
    results_rows = []
    display_rows = []
    start_total = time.time()
    n_passed = 0

    solv_label = solvent if solvent else "Fase Gassosa (Vuoto)"

    for key in target_keys:
        if key not in BENCHMARK_SYSTEMS:
            continue
        item = BENCHMARK_SYSTEMS[key]
        cfg = CalculationConfig(
            smiles=item["smiles"],
            n_conformers=15,
            rmsd_threshold=0.20,
            force_field=force_field,
            psi4_method=method,
            psi4_basis=basis,
            solvent=solvent,
            sp_only=sp_only,
            psi4_threads=2,
            psi4_memory="2 GB"
        )

        t0 = time.time()
        status = "ACCURATO (PASS)"
        computed_delta = 0.0
        geo_details = ""
        notes = "Calcolo completato con successo"
        res: Optional[CalculationResult] = None

        try:
            res = run_pipeline(cfg)
            elapsed = round(time.time() - t0, 1)

            if len(res.conformers) >= 2:
                computed_delta, notes, geo_details = _analyze_conformations_by_geometry(key, res, item)
                err = round(abs(computed_delta - item["reference_delta_e_kcalmol"]), 2)
                if err <= item["tolerance_kcalmol"]:
                    status = "ACCURATO (PASS)"
                    n_passed += 1
                elif err <= item["tolerance_kcalmol"] * 1.5:
                    status = "ACCETTABILE (LIEVE DEVIAZIONE)"
                    notes = f"Scostamento moderato ({err:.2f} kcal/mol vs tolleranza {item['tolerance_kcalmol']:.2f})"
                else:
                    status = "DEVIAZIONE ELEVATA (WARN)"
                    notes = f"Deviazione significativa ({err:.2f} kcal/mol vs tolleranza {item['tolerance_kcalmol']:.2f})"
            else:
                computed_delta = 0.0
                err = item["reference_delta_e_kcalmol"]
                status = "NON ISOLATO (WARN)"
                notes = "Solo 1 conformero isolato entro la finestra energetica"
                geo_details = "Nessun conformero secondario rilevato"

            results_rows.append({
                "Sistema": item["name"],
                "SMILES": item["smiles"],
                "Processo_Conformazionale": item["process_description"],
                "Conformero_Rif_A": item["conformation_a"],
                "Conformero_Rif_B": item["conformation_b"],
                "Metodo_Calcolo": f"{method.upper()}/{basis} ({solv_label})",
                "Campo_Forze": force_field.upper(),
                "Delta_E_Calc_kcal_mol": computed_delta,
                "Delta_E_Ref_kcal_mol": item["reference_delta_e_kcalmol"],
                "Livello_Riferimento": item["reference_level"],
                "Errore_Assoluto_kcal_mol": err,
                "Tolleranza_kcal_mol": item["tolerance_kcalmol"],
                "Esito_Validazione": status,
                "Tempo_s": elapsed,
                "Note_Geometriche": geo_details,
                "Citazione_Letteratura": item["citation"],
                "Significato_Scientifico": item["scientific_significance"]
            })

            display_rows.append({
                "Sistema": item["name"],
                "Processo": item["process_description"].split(":")[0],
                "ΔE Calc (kcal/mol)": f"{computed_delta:.2f}",
                "ΔE Rif (kcal/mol)": f"{item['reference_delta_e_kcalmol']:.2f}",
                "Errore (kcal/mol)": f"{err:.2f}",
                "Tolleranza": f"±{item['tolerance_kcalmol']:.2f}",
                "Esito": "✅ PASS" if "PASS" in status else ("⚠️ DEVIAZIONE" if "DEVIAZIONE" in status else "⚠️ WARN")
            })

        except Exception as e:
            elapsed = round(time.time() - t0, 1)
            results_rows.append({
                "Sistema": item["name"],
                "SMILES": item["smiles"],
                "Processo_Conformazionale": item["process_description"],
                "Conformero_Rif_A": item["conformation_a"],
                "Conformero_Rif_B": item["conformation_b"],
                "Metodo_Calcolo": f"{method.upper()}/{basis} ({solv_label})",
                "Campo_Forze": force_field.upper(),
                "Delta_E_Calc_kcal_mol": 0.0,
                "Delta_E_Ref_kcal_mol": item["reference_delta_e_kcalmol"],
                "Livello_Riferimento": item["reference_level"],
                "Errore_Assoluto_kcal_mol": None,
                "Tolleranza_kcal_mol": item["tolerance_kcalmol"],
                "Esito_Validazione": "FALLITO (ERRORE)",
                "Tempo_s": elapsed,
                "Note_Geometriche": f"Eccezione: {e}",
                "Citazione_Letteratura": item["citation"],
                "Significato_Scientifico": item["scientific_significance"]
            })
            display_rows.append({
                "Sistema": item["name"],
                "Processo": item["process_description"].split(":")[0],
                "ΔE Calc (kcal/mol)": "-",
                "ΔE Rif (kcal/mol)": f"{item['reference_delta_e_kcalmol']:.2f}",
                "Errore (kcal/mol)": "-",
                "Tolleranza": f"±{item['tolerance_kcalmol']:.2f}",
                "Esito": "❌ ERROR"
            })

    total_time = round(time.time() - start_total, 1)
    df_results = pd.DataFrame(results_rows)
    df_display = pd.DataFrame(display_rows)
    mae = df_results["Errore_Assoluto_kcal_mol"].dropna().mean() if not df_results["Errore_Assoluto_kcal_mol"].dropna().empty else 0.0

    bench_dict = {
        "df_results": df_results,
        "df_display": df_display,
        "mean_absolute_error_kcalmol": round(float(mae), 3),
        "total_elapsed_s": total_time,
        "method": method.upper(),
        "basis": basis,
        "force_field": force_field.upper(),
        "solvent": solv_label,
        "n_systems": len(target_keys),
        "n_passed": n_passed,
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }

    # Generate standardized CSV and Markdown exports
    bench_dict["csv_content"] = generate_benchmark_csv(bench_dict)
    bench_dict["markdown_report"] = generate_benchmark_markdown(bench_dict)

    return bench_dict


def generate_benchmark_csv(bench_dict: Dict[str, Any]) -> str:
    """Produces an informative, standardized CSV report for the benchmark."""
    df: pd.DataFrame = bench_dict["df_results"]
    mae = bench_dict["mean_absolute_error_kcalmol"]
    n_sys = bench_dict["n_systems"]
    n_pass = bench_dict["n_passed"]
    verdict = "VALIDATO (ELEVATA ACCURATEZZA)" if mae <= 0.50 else "DEVIAZIONE MODERATA"

    header_lines = [
        "# ==============================================================================",
        "# ENERGYINATOR - REPORT DI VALIDAZIONE METODOLOGICA E BENCHMARK CONFORMAZIONALE",
        f"# Data Esecuzione: {bench_dict.get('timestamp', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))}",
        f"# Metodo Quantistico: {bench_dict.get('method', 'B3LYP')}/{bench_dict.get('basis', '6-31G(d)')}",
        f"# Campo di Forze Preliminare: {bench_dict.get('force_field', 'MMFF94')}",
        f"# Fase / Solvente: {bench_dict.get('solvent', 'Fase Gassosa (Vuoto)')}",
        f"# Errore Assoluto Medio (MAE): {mae:.2f} kcal/mol",
        f"# Sistemi Testati: {n_sys} | Superati (Entro Tolleranza): {n_pass}",
        f"# Verdetto Complessivo: {verdict}",
        "# =============================================================================="
    ]
    csv_body = df.to_csv(index=False)
    return "\n".join(header_lines) + "\n" + csv_body


def generate_benchmark_markdown(bench_dict: Dict[str, Any]) -> str:
    """Generates an academic, publication-ready Markdown report for the benchmark."""
    df: pd.DataFrame = bench_dict["df_results"]
    mae = bench_dict["mean_absolute_error_kcalmol"]
    timestamp = bench_dict.get("timestamp", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    method = bench_dict.get("method", "B3LYP")
    basis = bench_dict.get("basis", "6-31G(d)")
    ff = bench_dict.get("force_field", "MMFF94")
    solv = bench_dict.get("solvent", "Fase Gassosa (Vuoto)")

    md = []
    md.append(f"# Report di Validazione Metodologica & Benchmark Conformazionale")
    md.append(f"**Software:** Energyinator Molecular Suite &bull; **Data:** {timestamp}  ")
    md.append(f"**Protocollo Computazionale:** `{method}/{basis}` ({solv}) &bull; **Campo di Forze:** `{ff}`  ")
    md.append(f"**Errore Medio Assoluto (MAE):** `{mae:.2f} kcal/mol`  ")
    md.append("")
    md.append("## 1. Sommario Tabellare dei Risultati")
    md.append("")
    md.append("| Sistema | Processo Valutato | ΔE Calc (kcal/mol) | ΔE Rif (kcal/mol) | Livello Riferimento | Errore (kcal/mol) | Esito |")
    md.append("| :--- | :--- | :---: | :---: | :--- | :---: | :---: |")

    for _, r in df.iterrows():
        err_str = f"{r['Errore_Assoluto_kcal_mol']:.2f}" if r['Errore_Assoluto_kcal_mol'] is not None else "-"
        md.append(
            f"| **{r['Sistema']}** | {r['Processo_Conformazionale'].split(':')[0]} | "
            f"{r['Delta_E_Calc_kcal_mol']:.2f} | {r['Delta_E_Ref_kcal_mol']:.2f} | "
            f"{r['Livello_Riferimento']} | {err_str} | {r['Esito_Validazione']} |"
        )

    md.append("")
    md.append("## 2. Dettaglio dei Sistemi Molecolari & Rilevanza Chimico-Fisica")
    md.append("")

    for _, r in df.iterrows():
        md.append(f"### {r['Sistema']}")
        md.append(f"- **Descrizione Processo:** {r['Processo_Conformazionale']}")
        md.append(f"- **Stato di Minimo (A):** {r['Conformero_Rif_A']}")
        md.append(f"- **Stato Superiore (B):** {r['Conformero_Rif_B']}")
        md.append(f"- **Dati Geometrici Rilevati:** {r['Note_Geometriche']}")
        md.append(f"- **Significato Chimico:** {r['Significato_Scientifico']}")
        md.append(f"- **Fonte Bibliografica:** *{r['Citazione_Letteratura']}*")
        err_str = f"{r['Errore_Assoluto_kcal_mol']:.2f}" if r['Errore_Assoluto_kcal_mol'] is not None else "N/A"
        md.append(f"- **Accuratezza Metodo:** Errore = `{err_str} kcal/mol` (Tolleranza ammissibile: `±{r['Tolleranza_kcal_mol']:.2f} kcal/mol`)")
        md.append("")

    md.append("## 3. Considerazioni Metodologiche per la Tesi / Pubblicazione")
    if mae <= 0.50:
        md.append(
            f"Il protocollo `{method}/{basis}` dimostra un'elevata accuratezza conformazionale "
            f"(MAE = {mae:.2f} kcal/mol < 0.50 kcal/mol), confermando l'affidabilità metodologica "
            f"per l'esplorazione conformazionale e la stima delle popolazioni di Boltzmann."
        )
    else:
        md.append(
            f"Il protocollo `{method}/{basis}` mostra una deviazione media di {mae:.2f} kcal/mol. "
            f"Per una maggiore accuratezza su sistemi flessibili, si consiglia l'impiego di una base "
            f"con funzioni di polarizzazione e diffusione (es. 6-311+G(d,p) o def2-TZVP) o l'inclusione "
            f"di correzioni di dispersione empirica (es. D3BJ)."
        )

    return "\n".join(md)


def main():
    print("=" * 70)
    print(" ENERGYINATOR - Conformational Benchmark & Validation Suite")
    print("=" * 70)
    bench = run_benchmark(method="b3lyp", basis="3-21G", sp_only=True)
    print("\nRISULTATI DEL BENCHMARK:")
    print(bench["df_display"].to_string(index=False))
    print("-" * 70)
    print(f"Errore Medio Assoluto (MAE): {bench['mean_absolute_error_kcalmol']} kcal/mol")
    print(f"Tempo totale esecuzione:    {bench['total_elapsed_s']} s")
    print("=" * 70)


if __name__ == "__main__":
    main()
