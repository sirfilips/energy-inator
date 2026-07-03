import os
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem
import psi4
import pandas as pd
import time
import logging
from tqdm import tqdm
import sys

# =========================
# Configurazione globale Psi4
# =========================
PSI4_MEMORY = '8 GB'
PSI4_THREADS = 8
PSI4_METHOD = 'b3lyp'
PSI4_BASIS = '6-31G(d)'

# Usa script_dir già presente nello script
script_dir = os.path.dirname(os.path.abspath(__file__))
log_path = os.path.join(script_dir, "energyLOG.log")

logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(log_path, encoding="utf-8"),
        logging.StreamHandler(sys.stdout)
    ]
)

def genera_conformeri(smiles, n_conformers=50, rmsd_threshold=0.2):
    molecola = Chem.MolFromSmiles(smiles)
    if molecola is None:
        raise ValueError("SMILES non valido o impossibile da convertire in molecola.")
    molecola = Chem.AddHs(molecola)
    logging.info(f"Generazione di {n_conformers} conformeri per SMILES: {smiles}")
    ids = AllChem.EmbedMultipleConfs(
        molecola,
        numConfs=int(n_conformers),
        pruneRmsThresh=rmsd_threshold,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True
    )
    # EmbedMultipleConfs può restituire una tupla o lista di id; assicurati di usare la sua lunghezza reale
    if not ids:
        logging.error("Impossibile generare conformeri 3D della molecola.")
        raise ValueError("Impossibile generare conformeri 3D della molecola.")
    logging.info(f"Generati {len(ids)} conformeri (richiesti: {n_conformers})")
    return molecola, ids

def psi4_optimize_worker(args):
    conf_xyz, conf_id, output_dir, main_timestamp, smiles = args
    # Usa il modulo psi4 già importato a livello globale; non sovrascrivere memoria/threads qui
    # psi4.set_memory e psi4.set_num_threads sono definiti globalmente in calcola_energia_psi4
    try:
        mol = psi4.geometry(conf_xyz)
        energia = psi4.optimize(PSI4_METHOD, molecule=mol)
        mol = psi4.core.get_active_molecule()
        if mol is None:
            logging.debug(f"Psi4 non ha restituito una molecola per conf {conf_id}")
            psi4.core.clean()
            return None
        xyz_filename = f"geom_ottimizzata_{main_timestamp}_{conf_id}.xyz"
        output_path = os.path.join(output_dir, xyz_filename)
        natoms = mol.natom()
        bohr_to_angstrom = 0.52917721092
        coords = []
        for i in range(natoms):
            symbol = mol.symbol(i)
            x = mol.x(i) * bohr_to_angstrom
            y = mol.y(i) * bohr_to_angstrom
            z = mol.z(i) * bohr_to_angstrom
            coords.append(f"{symbol} {x:.8f} {y:.8f} {z:.8f}")
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(f"{natoms}\n")
            f.write(f"{smiles}\n")
            for line in coords:
                f.write(line + "\n")
        hartree_to_kcalmol = 627.509
        hartree_to_ev = 27.2114
        energia_kcalmol = energia * hartree_to_kcalmol
        energia_ev = energia * hartree_to_ev
        xyz_ottimizzato = f"{natoms}\n{smiles}\n" + "\n".join(coords)
        psi4.core.clean()
        return {
            "Conformer_ID": conf_id,
            "Energy_Hartree": energia,
            "Energy_kcalmol": energia_kcalmol,
            "Energy_eV": energia_ev,
            "xyz": conf_xyz,
            "xyz_ottimizzato": xyz_ottimizzato,
            "output_path": output_path
        }
    except Exception as e:
        logging.debug(f"Errore Psi4 conf {conf_id}: {repr(e)}")
        logging.debug(f"Input XYZ:\n{conf_xyz}")
        psi4.core.clean()
        return None

def ottimizza_conformeri(molecola, ids, output_dir, main_timestamp, top_n=10, rmsd_threshold=0.2):
    min_energia = None
    max_energia = None
    best_xyz = None
    worst_xyz = None
    best_output_path = None
    worst_output_path = None
    best_energia_kcalmol = None
    best_energia_ev = None
    best_xyz_ottimizzato = None
    worst_energia_kcalmol = None
    worst_energia_ev = None
    worst_xyz_ottimizzato = None
    conformers_data = []

    # 1. Ottimizzazione UFF/MMFF e calcolo energia (ESEGUIUTO SEQUENZIALMENTE per evitare race condition)
    def uff_optimize(conf_id):
        try:
            # UFFOptimizeMolecule modifica la molecola in-place per il confId passato
            AllChem.UFFOptimizeMolecule(molecola, confId=conf_id)
            props = AllChem.UFFGetMoleculeForceField(molecola, confId=conf_id)
            energia = props.CalcEnergy()
            return (conf_id, energia)
        except Exception as e:
            logging.warning(f"Errore UFFOptimizeMolecule per conf {conf_id}: {e}")
            return None

    results = []
    for conf_id in ids:
        r = uff_optimize(conf_id)
        if r is not None:
            results.append(r)

    # 2. Filtra duplicati RMSD
    unique_ids = []
    ref_confs = []
    for conf_id, uff_energy in sorted(results, key=lambda x: x[1]):
        is_unique = True
        for ref_id in ref_confs:
            rms = AllChem.GetConformerRMS(molecola, conf_id, ref_id)
            if rms < rmsd_threshold:
                is_unique = False
                break
        if is_unique:
            unique_ids.append((conf_id, uff_energy))
            ref_confs.append(conf_id)
        if len(unique_ids) >= top_n:
            break

    # 3. Prepara input per Psi4 (sequenziale)
    psi4_inputs = []
    # Usa lo SMILES originario senza H espliciti per la label del file
    smiles_label = Chem.MolToSmiles(Chem.RemoveHs(molecola))
    for conf_id, uff_energy in unique_ids:
        conf = molecola.GetConformer(conf_id)
        num_atomi = molecola.GetNumAtoms()
        xyz = f"{num_atomi}\nGenerated by Energyinator\n"
        for i in range(num_atomi):
            posizione = conf.GetAtomPosition(i)
            simbolo = molecola.GetAtomWithIdx(i).GetSymbol()
            xyz += f"{simbolo} {posizione.x:.6f} {posizione.y:.6f} {posizione.z:.6f}\n"
        psi4_inputs.append((xyz, conf_id, output_dir, main_timestamp, smiles_label))

    # 4. Ottimizzazione Psi4 sequenziale
    psi4_results = []
    for args in tqdm(psi4_inputs, desc="Psi4 (sequenziale)", unit="conf"):
        psi4_results.append(psi4_optimize_worker(args))

    # 5. Raccolta risultati
    for idx, result in enumerate(psi4_results):
        if result is None:
            continue
        conf_id = result["Conformer_ID"]
        uff_energy = unique_ids[idx][1]
        energia = result["Energy_Hartree"]
        energia_kcalmol = result["Energy_kcalmol"]
        energia_ev = result["Energy_eV"]
        xyz = result["xyz"]
        xyz_ottimizzato = result["xyz_ottimizzato"]
        output_path = result["output_path"]

        conformers_data.append({
            "Conformer_ID": conf_id,
            "UFF_Energy": uff_energy,
            "Energy_Hartree": energia,
            "Energy_kcalmol": energia_kcalmol,
            "Energy_eV": energia_ev
        })

        if (min_energia is None) or (energia < min_energia):
            min_energia = energia
            best_xyz = xyz
            best_output_path = output_path
            best_energia_kcalmol = energia_kcalmol
            best_energia_ev = energia_ev
            best_xyz_ottimizzato = xyz_ottimizzato

        if (max_energia is None) or (energia > max_energia):
            max_energia = energia
            worst_xyz = xyz
            worst_output_path = output_path
            worst_energia_kcalmol = energia_kcalmol
            worst_energia_ev = energia_ev
            worst_xyz_ottimizzato = xyz_ottimizzato

    return (
        min_energia, best_energia_kcalmol, best_energia_ev, best_xyz_ottimizzato, best_output_path,
        max_energia, worst_energia_kcalmol, worst_energia_ev, worst_xyz_ottimizzato, worst_output_path,
        conformers_data
    )

def calcola_energia_psi4(smiles, n_conformers=50, rmsd_threshold=0.2):
    start_time = time.time()
    # Impostazioni globali per psi4: queste saranno usate da tutti i worker
    psi4.set_memory(PSI4_MEMORY)
    psi4.set_num_threads(PSI4_THREADS)
    psi4.set_options({'basis': PSI4_BASIS})

    script_dir = os.path.dirname(os.path.abspath(__file__))
    main_timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(script_dir, "output", f"result_{main_timestamp}")
    os.makedirs(output_dir, exist_ok=True)

    molecola, ids = genera_conformeri(smiles, n_conformers, rmsd_threshold)
    (
        min_energia, best_energia_kcalmol, best_energia_ev, best_xyz_ottimizzato, best_output_path,
        max_energia, worst_energia_kcalmol, worst_energia_ev, worst_xyz_ottimizzato, worst_output_path,
        conformers_data
    ) = ottimizza_conformeri(molecola, ids, output_dir, main_timestamp)

    if min_energia is None:
        raise RuntimeError("Nessun conformero è stato ottimizzato con successo.")

    if conformers_data:
        df = pd.DataFrame(conformers_data)
        # Ordina per Conformer_ID crescente
        df = df.sort_values(by="Conformer_ID")
        # Porta la colonna UFF_Energy in fondo
        cols = [col for col in df.columns if col != "UFF_Energy"] + ["UFF_Energy"]
        df = df[cols]
        excel_path = os.path.join(output_dir, f"conformers_energies_{main_timestamp}.xlsx")
        df.to_excel(excel_path, index=False)

    elapsed_time = time.time() - start_time

    return (
        min_energia, best_energia_kcalmol, best_energia_ev, best_xyz_ottimizzato, best_output_path,
        max_energia, worst_energia_kcalmol, worst_energia_ev, worst_xyz_ottimizzato, worst_output_path,
        elapsed_time
    )



if __name__ == "__main__":
    codice_smiles = input("Inserisci il codice SMILES della molecola: ")
    try:
        (
            energia_min, energia_min_kcalmol, energia_min_ev, xyz_min_ottimizzata, output_path_min,
            energia_max, energia_max_kcalmol, energia_max_ev, xyz_max_ottimizzata, output_path_max,
            tempo_totale
        ) = calcola_energia_psi4(codice_smiles)
        print(f"L'energia MINIMA della molecola è: {energia_min:.6f} Hartree")
        print(f"In kcal/mol: {energia_min_kcalmol:.2f} kcal/mol")
        print(f"In eV: {energia_min_ev:.2f} eV")
        print("\nGeometria ottimizzata (min, XYZ):\n")
        print(xyz_min_ottimizzata)
        print(f"\nGeometria ottimizzata (min) salvata in: {output_path_min}")

        print("\n" + "="*60 + "\n")

        print(f"L'energia MASSIMA della molecola è: {energia_max:.6f} Hartree")
        print(f"In kcal/mol: {energia_max_kcalmol:.2f} kcal/mol")
        print(f"In eV: {energia_max_ev:.2f} eV")
        print("\nGeometria ottimizzata (max, XYZ):\n")
        print(xyz_max_ottimizzata)
        print(f"\nGeometria ottimizzata (max) salvata in: {output_path_max}")

        print(f"\nTempo totale impiegato per i calcoli: {tempo_totale:.2f} secondi")

    except Exception as e:
        logging.error("Errore durante il calcolo o la scrittura della geometria:")
        logging.error(e)
