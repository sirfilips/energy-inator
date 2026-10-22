import os
import sys
import argparse
import logging
from energyinator.config import CalculationConfig
from energyinator.pipeline import run_pipeline

# Initialize console logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="energyinator",
        description="Energyinator: Framework ibrido per campionamento conformazionale e chimica quantistica (RDKit + Psi4)."
    )
    parser.add_argument(
        "smiles",
        nargs="?",
        default=None,
        help="Stringa SMILES della molecola bersaglio (se omessa, verrà richiesta interattivamente)."
    )
    parser.add_argument(
        "-f", "--input-file",
        type=str,
        default=None,
        help="File contenente lista di stringhe SMILES per elaborazione batch (un SMILES per riga)."
    )
    parser.add_argument(
        "-n", "--n-conformers",
        type=int,
        default=50,
        help="Numero di conformeri iniziali da generare con ETKDG (default: 50)."
    )
    parser.add_argument(
        "-r", "--rmsd-threshold",
        type=float,
        default=0.2,
        help="Soglia RMSD (in Å) per il filtraggio classico dei duplicati (default: 0.2)."
    )
    parser.add_argument(
        "--post-dft-rmsd",
        type=float,
        default=0.15,
        help="Soglia RMSD (in Å) post-DFT per deduplicare conformeri collassati (default: 0.15)."
    )
    parser.add_argument(
        "--ff",
        type=str,
        choices=["mmff94", "mmff94s", "uff"],
        default="mmff94",
        help="Campo di forze classico per il pre-rilassamento (default: mmff94)."
    )
    parser.add_argument(
        "-w", "--energy-window",
        type=float,
        default=10.0,
        help="Finestra energetica FF in kcal/mol: scarta conformeri ad alta tensione prima di Psi4 (default: 10.0, 0 per disattivare)."
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Numero massimo di conformeri unici da calcolare con Psi4 (default: calcola tutti quelli entro la finestra)."
    )
    parser.add_argument(
        "-m", "--method",
        type=str,
        default="b3lyp",
        help="Metodo quantomeccanico DFT / ab-initio per Psi4 (default: b3lyp; opzioni: wb97x-d, pbe, etc.)."
    )
    parser.add_argument(
        "-b", "--basis",
        type=str,
        default="6-31G(d)",
        help="Basis set per i calcoli quantistici Psi4 (default: 6-31G(d))."
    )
    parser.add_argument(
        "--solvent",
        type=str,
        default=None,
        help="Solvente implicito (modello PCM, es. water, chloroform, dmso, methanol, toluene)."
    )
    parser.add_argument(
        "--frequencies",
        action="store_true",
        help="Calcola frequenze vibrazionali e correzioni termiche (Gibbs, ZPVE, Entropia) sul ground state."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=min(8, os.cpu_count() or 4),
        help=f"Numero di thread CPU per calcoli paralleli (default: {min(8, os.cpu_count() or 4)})."
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="8 GB",
        help="Memoria RAM massima allocata a Psi4 (default: '8 GB')."
    )
    parser.add_argument(
        "--sp-only",
        action="store_true",
        help="Esegue solo calcolo Single Point invece dell'ottimizzazione completa di geometria."
    )
    parser.add_argument(
        "--charge",
        type=int,
        default=None,
        help="Forza una specifica carica formale (se non specificata, calcolata automaticamente da RDKit)."
    )
    parser.add_argument(
        "--multiplicity",
        type=int,
        default=None,
        help="Forza una specifica molteplicità di spin (se non specificata, calcolata automaticamente)."
    )
    parser.add_argument(
        "--temp",
        type=float,
        default=298.15,
        help="Temperatura in Kelvin per la distribuzione di Boltzmann (default: 298.15 K)."
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default=None,
        help="Directory personalizzata dove salvare i file generati."
    )
    return parser.parse_args()


def print_summary(res):
    print("\n" + "=" * 65)
    print(f"RISULTATI OTTIMIZZAZIONE QUANTISTICA PER: {res.smiles}")
    print("=" * 65)
    print(f"Stato: Carica = {res.formal_charge}, Molteplicità = {res.spin_multiplicity}")
    print(f"Livello teorico: {res.method.upper()} / {res.basis}")
    if res.solvent:
        print(f"Solvatazione (PCM): {res.solvent.capitalize()}")
    print(f"Campo di forza: {res.force_field.upper()}")
    print(f"Conformeri generati: {res.total_generated} | Filtrati FF: {res.unique_ff} | Calcolati DFT: {res.calculated_dft} | Unici DFT: {res.unique_dft}")
    print("-" * 65)
    print(f"Energia MINIMA (Minimo globale):  {res.min_energy_hartree:.6f} Hartree")
    print(f"                             --> {res.min_energy_kcalmol:.2f} kcal/mol")
    print(f"                             --> {res.min_energy_ev:.2f} eV")
    print(f"Geometria minima salvata in:      {res.best_output_path}")
    print("-" * 65)
    print(f"Energia MASSIMA:                 {res.max_energy_hartree:.6f} Hartree")
    print(f"                             --> {res.max_energy_kcalmol:.2f} kcal/mol")
    print(f"                             --> {res.max_energy_ev:.2f} eV")
    print(f"Delta E max-min:                 {res.max_energy_kcalmol - res.min_energy_kcalmol:.2f} kcal/mol")
    print("-" * 65)
    if res.excel_path:
        print(f"Tabella completa Excel:          {res.excel_path}")
    print(f"Tempo totale impiegato:          {res.elapsed_time:.2f} secondi")
    print("=" * 65 + "\n")


def main():
    args = parse_args()
    
    # Process batch file if specified
    if args.input_file:
        if not os.path.isfile(args.input_file):
            print(f"Errore: File '{args.input_file}' non trovato.")
            sys.exit(1)
        with open(args.input_file, "r", encoding="utf-8") as f:
            smiles_list = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        
        print(f"Avvio elaborazione batch per {len(smiles_list)} molecole...")
        for idx, smi in enumerate(smiles_list, 1):
            print(f"\n[{idx}/{len(smiles_list)}] Elaborazione SMILES: {smi}")
            config = CalculationConfig(
                smiles=smi,
                n_conformers=args.n_conformers,
                rmsd_threshold=args.rmsd_threshold,
                post_dft_rmsd=args.post_dft_rmsd,
                force_field=args.ff,
                ff_energy_window=args.energy_window if args.energy_window > 0 else None,
                top_n=args.top_n,
                psi4_memory=args.memory,
                psi4_threads=args.threads,
                psi4_method=args.method,
                psi4_basis=args.basis,
                sp_only=args.sp_only,
                solvent=args.solvent,
                calc_frequencies=args.frequencies,
                charge=args.charge,
                multiplicity=args.multiplicity,
                temperature_k=args.temp,
                output_dir=args.output_dir
            )
            try:
                res = run_pipeline(config)
                print_summary(res)
            except Exception as e:
                print(f"Errore durante l'elaborazione di {smi}: {e}")
        return

    # Single molecule run
    smiles = args.smiles
    if not smiles:
        smiles = input("Inserisci il codice SMILES della molecola: ").strip()
        if not smiles:
            print("Errore: Nessun codice SMILES fornito.")
            sys.exit(1)

    config = CalculationConfig(
        smiles=smiles,
        n_conformers=args.n_conformers,
        rmsd_threshold=args.rmsd_threshold,
        post_dft_rmsd=args.post_dft_rmsd,
        force_field=args.ff,
        ff_energy_window=args.energy_window if args.energy_window > 0 else None,
        top_n=args.top_n,
        psi4_memory=args.memory,
        psi4_threads=args.threads,
        psi4_method=args.method,
        psi4_basis=args.basis,
        sp_only=args.sp_only,
        solvent=args.solvent,
        calc_frequencies=args.frequencies,
        charge=args.charge,
        multiplicity=args.multiplicity,
        temperature_k=args.temp,
        output_dir=args.output_dir
    )

    try:
        res = run_pipeline(config)
        print_summary(res)
    except Exception as e:
        print(f"\nErrore durante il calcolo: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
