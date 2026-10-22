# Energyinator

**Framework Ibrido ad Alta Efficienza per il Campionamento Conformazionale e l'Ottimizzazione Quantomeccanica**

Energyinator è un framework Python progettato per generare conformeri molecolari 3D partendo da codici SMILES e calcolarne le energie tramite metodi di chimica quantistica. Il protocollo combina la velocità e l'efficacia di **RDKit** (algoritmo ETKDGv3 e campi di forza MMFF94/UFF) con l'accuratezza ab-initio/DFT di **Psi4**.

Include un'interfaccia grafica moderna in **Streamlit** che consente di visualizzare e sovrapporre conformeri in 3D, misurare parametri geometrici (senza dover ricorrere a software esterni come Avogadro), calcolare la curva continua delle barriere di potenziale torsionali (1D PES Scan) e analizzare le popolazioni termodinamiche di Boltzmann e le proprietà elettroniche molecolari.

---

## Caratteristiche Principali

- **Generazione Stocastica Multithread**: Campionamento conformazionale accelerato con algoritmo ETKDGv3 (`rdkit`).
- **Supporto Campi di Forza Avanzati**: Pre-rilassamento geometrico con **MMFF94** (o MMFF94s), con fallback automatico su **UFF**.
- **Screening Energetico a Imbuto (Energy Window)**: Scarto preventivo delle conformazioni ad alta tensione classica per abbattere i tempi computazionali quantistici.
- **Deduplicazione RMSD con Simmetria**: Riconoscimento rigoroso dell'automorfismo chimico (`GetBestRMS`) per eliminare strutture spazialmente ridondanti.
- **Riconoscimento Automatico di Carica e Molteplicità**: Rilevamento dello stato di carica formale e del numero di elettroni spaiati direttamente dalla struttura molecolare.
- **Calcolo Quantomeccanico con Psi4**: Supporto per DFT ibrida (`B3LYP`, funzionali con dispersione empirica come `wB97X-D`, `PBE`, ecc.) e ab-initio (`HF`, `MP2`).
- **Proprietà Elettroniche e Polarità**: Estrazione automatica della funzione d'onda SCF del **momento di dipolo elettrico ($\mu$ in Debye)**, delle energie orbitali **HOMO e LUMO** e del relativo **gap energetico** ($\Delta E_{gap}$).
- **Scansione Torsionale 1D (PES Scan)**: Calcolo continuo della superficie di energia potenziale da 0° a 360° attorno a qualsiasi angolo diedro selezionato per quantificare l'altezza della barriera rotazionale ($\Delta E^\ddagger$) sia classica (`MMFF94`/`UFF`) che quantistica (`Psi4`).
- **Sovrapposizione 3D (Conformer Overlay)**: Allineamento geometrico tramite algoritmo di rotazione di Kabsch (SVD) per sovrapporre qualsiasi conformero eccitato al ground state e visualizzare a colpo d'occhio i legami ruotati.
- **Analisi Termodinamica e Popolazioni di Boltzmann**: Calcolo automatico dei $\Delta E$ (in kcal/mol e kJ/mol), identificazione degli stati degeneri ($g_i$) e calcolo della distribuzione di Boltzmann a 298.15 K ($P_i$).
- **Garanzia di Comparabilità (Profilo Legacy)**: Profili di calcolo pre-configurati (tra cui il **Profilo Legacy v1.0**) per garantire che ogni nuova elaborazione sia rigorosamente identica e confrontabile con le serie storiche pubblicate.
- **Architettura Stabile a Processo Isolato**: Esecuzione del motore quantistico in un processo worker dedicato per prevenire conflitti OpenMP e access violation su Windows.
- **Esportazione Chiara e Dettagliata**: Geometrie 3D ottimizzate in formato `.xyz`, foglio `.xlsx` per l'analisi tabellare e file strutturato `.json` per il ripristino istantaneo delle sessioni.

---

## Requisiti di Sistema

Psi4 richiede l'installazione tramite l'ecosistema Conda.

| Requisito | Minimo | Consigliato |
| :--- | :--- | :--- |
| **Python** | 3.10 | 3.10 – 3.14 |
| **Gestore Pacchetti** | Miniconda o Anaconda | Miniconda |
| **RAM** | 8 GB | 16 GB |
| **CPU** | Quad-Core | 8+ Core (per multithreading Psi4 e RDKit) |
| **Spazio su Disco** | 5 GB | 10 GB |

---

## Installazione Rapida

### 1. Clona il repository
```bash
git clone https://github.com/sirfilips/energy-inator.git
cd energy-inator
```

### 2. Crea l'ambiente Conda dedicato
Puoi creare direttamente l'ambiente locale utilizzando il file `environment.yml` fornito nel repository:

```bash
conda env create --prefix .\.venv -f environment.yml
conda activate .\.venv
```

Oppure in un ambiente denominato globale:
```bash
conda env create -f environment.yml
conda activate energyinator
```

---

## Utilizzo

### Interfaccia Grafica Web (Streamlit Studio)
Puoi avviare l'applicazione grafica completa semplicemente facendo doppio clic su:
`avvia_gui.bat`

Oppure da terminale:
```bash
streamlit run app_gui.py
# oppure con ambiente locale
.\.venv\python.exe -m streamlit run app_gui.py
```

L'interfaccia è strutturata in 5 ambienti specializzati:
1. **Struttura 3D e Geometria**:
   - Visualizzatore molecolare WebGL interattivo (rotazione, zoom e traslazione).
   - Modalità **Sovrapponi al Minimo Globale (Overlay)** con allineamento Kabsch (ground state in ciano trasparente vs conformero attivo in CPK).
   - **Ispettore Geometrico (Anti-Avogadro)**: misurazione istantanea di distanze di legame (Å), angoli di legame (°) e angoli diedri/torsionali (°).
   - Lettura del momento di dipolo ($\mu$), HOMO, LUMO e gap energetico del conformero attivo.
2. **Scansione Torsionale (PES)**:
   - Ambiente a schermo intero per calcolare la barriera di rotazione continua attorno a un angolo diedro a 4 atomi.
   - Grafico continuo della curva di potenziale $\Delta E(\theta)$ con indicazione del minimo e del picco della barriera ($\Delta E^\ddagger$).
   - Slider interattivo per scorrere lungo la curva e vedere la molecola 3D ruotare in diretta a ogni angolo.
3. **Analisi Energetica e Popolazioni**:
   - Metriche di sintesi ($E_{min}$, $E_{max}$, $\Delta E_{max-min}$, conformeri unici DFT).
   - Grafico dello spettro delle energie relative $\Delta E$ in kcal/mol.
   - Diagramma a ciambella delle popolazioni di Boltzmann all'equilibrio termodinamico a 298.15 K.
   - Grafico comparativo del momento di dipolo elettrico nei conformeri.
4. **Tabella Dati ed Esportazione**:
   - Tabella riassuntiva di tutti i conformeri con dati quantistici, geometrici ed elettronici.
   - Download con 1 clic del report Excel (`.xlsx`) e del file coordinate `.xyz` del minimo globale.
   - Tasto per aprire direttamente la cartella dei risultati in Windows Explorer.
5. **Log di Calcolo**:
   - Visualizzatore del registro di calcolo (`run.log`) e dell'output scientifico completo di Psi4 (`psi4_output.dat`).

---

### Modalità Interattiva da Terminale
Se lanciato senza argomenti, il programma chiederà la stringa SMILES:
```bash
python Energyinator.py
```

### Modalità Riga di Comando (CLI)
Puoi specificare direttamente la molecola e personalizzare tutti i parametri di calcolo:

```bash
# Esempio su butano con B3LYP e base 6-31G(d)
python Energyinator.py "CCCC"

# Esempio con campo di forza MMFF94, finestra energetica di 8 kcal/mol e 8 thread
python Energyinator.py "CCCCCC" --ff mmff94 -w 8.0 -t 8

# Screening rapido Single-Point (senza ottimizzazione geometrica)
python Energyinator.py "c1ccccc1" --sp-only

# Utilizzo di funzionale con correzione di dispersione nativa
python Energyinator.py "c1ccccc1-c2ccccc2" -m wb97x-d -b def2-SVP
```

### Elaborazione Batch (Multi-molecola)
Puoi processare una serie di molecole partendo da un file `.smi` o `.txt` contenente uno SMILES per riga:
```bash
python Energyinator.py -f molecole.smi --n-conformers 30 --ff mmff94
```

### Utilizzo Programmatico in Python

```python
from energyinator import CalculationConfig, run_pipeline, run_dihedral_scan

# 1. Pipeline di conformer search e chimica quantistica
config = CalculationConfig(
    smiles="CCCC",              # Butano
    n_conformers=30,
    rmsd_threshold=0.2,
    force_field="mmff94",
    ff_energy_window=10.0,      # Filtra conformeri FF con Delta E > 10 kcal/mol
    psi4_method="b3lyp",
    psi4_basis="6-31G(d)",
    psi4_threads=8,
    temperature_k=298.15        # Calcolo popolazioni di Boltzmann
)

result = run_pipeline(config)

print(f"Minimo globale: {result.min_energy_hartree:.6f} Eh ({result.min_energy_kcalmol:.2f} kcal/mol)")
for conf in result.conformers:
    print(f"Conf {conf.conf_id}: Delta E = {conf.delta_dft_kcalmol:.2f} kcal/mol | Pop. = {conf.boltzmann_pop:.1f}% | Dipolo = {conf.dipole_debye:.2f} D")

# 2. Scansione della barriera torsionale 1D (PES Scan)
pes = run_dihedral_scan(
    smiles="CCCC",
    xyz_str=result.best_xyz_optimized,
    atom_indices=(0, 1, 2, 3),  # Carboni C0-C1-C2-C3
    step_deg=15,
    engine="mmff94"
)
print(f"Barriera rotazionale massima: {pes['max_barrier_kcalmol']:.2f} kcal/mol (Minimo a {pes['min_angle']}°)")
```

---

## Opzioni della Riga di Comando

| Parametro | Breve | Default | Descrizione |
| :--- | :--- | :--- | :--- |
| `smiles` | Posizionale | `None` | Stringa SMILES della molecola bersaglio |
| `--input-file` | `-f` | `None` | File con elenco di stringhe SMILES (batch) |
| `--n-conformers` | `-n` | `50` | Numero di conformeri stocastici iniziali |
| `--rmsd-threshold`| `-r` | `0.2` | Soglia RMSD (Å) per filtraggio di simmetria |
| `--ff` | | `mmff94` | Campo di forze: `mmff94`, `mmff94s`, `uff` |
| `--energy-window` | `-w` | `10.0` | Finestra energetica FF in kcal/mol (0 per disattivare) |
| `--top-n` | | `None` | Limite massimo di conformeri unici per Psi4 |
| `--method` | `-m` | `b3lyp` | Metodo quantistico (`b3lyp`, `wb97x-d`, `pbe`, ecc.) |
| `--basis` | `-b` | `6-31G(d)` | Basis set quantistico per Psi4 |
| `--threads` | `-t` | `8` | Thread CPU per calcoli paralleli |
| `--memory` | | `'8 GB'` | Memoria RAM allocata a Psi4 |
| `--sp-only` | | `False` | Esegue solo Single-Point senza ottimizzazione geometrica |
| `--charge` | | `Auto` | Carica formale forzata (di default calcolata da RDKit) |
| `--multiplicity` | | `Auto` | Molteplicità di spin forzata (di default calcolata da RDKit) |
| `--temp` | | `298.15` | Temperatura in K per la statistica di Boltzmann |
| `--output-dir` | `-o` | `output/` | Directory personalizzata di salvataggio |

---

## Struttura dei Risultati

I risultati vengono organizzati automaticamente in una directory con marcatura oraria `output/result_YYYYMMDD_HHMMSS/`:
```text
output/result_20261021_174202/
├── geom_ottimizzata_20261021_174202_0.xyz   # Coordinate del conformero 0
├── geom_ottimizzata_20261021_174202_1.xyz   # Coordinate del minimo globale (anti)
├── geom_ottimizzata_20261021_174202_2.xyz   # Coordinate del conformero 2
├── conformers_energies_20261021_174202.xlsx # Tabella riassuntiva Excel completa
├── result.json                              # Dati serializzati per la GUI o per analisi esterne
├── run.log                                  # Registro di calcolo e avanzamento
└── psi4_output.dat                          # Output scientifico grezzo di Psi4
```

### Contenuto del File Excel
Il report generato contiene le seguenti colonne per ciascun conformero, ordinate dal minimo globale agli stati eccitati:
- `Conformer_ID`: Identificativo numerico del conformero.
- `Degeneracy`: Grado di degenerazione / molteplicità statistica dello stato ($g_i$).
- `Energy_Hartree`: Energia elettronica totale in Hartree ($E_h$).
- `Energy_kcalmol`: Energia in kcal/mol.
- `Energy_eV`: Energia in elettronvolt (eV).
- `Delta_E_kcalmol`: Energia relativa $\Delta E$ rispetto al minimo globale ($0.00$ per il ground state).
- `Delta_E_kJmol`: $\Delta E$ in kJ/mol.
- `Boltzmann_Pop_%`: Popolazione conformazionale all'equilibrio termodinamico a 298.15 K ($P_i$).
- `Dipole_Debye`: Momento di dipolo elettrico $\mu$ in Debye.
- `HOMO_eV`: Energia dell'orbitale molecolare occupato più alto (HOMO) in eV.
- `LUMO_eV`: Energia dell'orbitale molecolare virtuale più basso (LUMO) in eV.
- `Gap_eV`: Gap energetico HOMO-LUMO ($\Delta E_{gap}$) in eV.
- `FF_Energy_kcalmol`: Energia classica del campo di forze pre-ottimizzato.
- `Delta_FF_kcalmol`: Energia classica relativa nel campo di forze.
- `XYZ_File`: Nome del file `.xyz` corrispondente.

---

## Esecuzione dei Test

Il progetto include una suite di test unitari automatizzati (`pytest`) per verificare la validità dei modelli, del campionamento conformazionale, dell'analisi di Boltzmann, della scansione PES e della pipeline completa:

```bash
pytest tests/
# oppure con l'ambiente conda locale:
.\.venv\python.exe -m pytest tests/
```

---

## Licenza

Questo progetto è rilasciato sotto licenza **MIT**. Consulta il file `LICENSE` per maggiori dettagli.
