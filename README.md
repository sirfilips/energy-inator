# Energyinator

**Framework Ibrido ad Alta Efficienza per il Campionamento Conformazionale e l'Ottimizzazione Quantomeccanica**

Energyinator è un framework Python progettato per generare conformeri molecolari 3D partendo da codici SMILES e calcolarne le energie tramite metodi di chimica quantistica. Il protocollo combina la velocità e l'efficacia di **RDKit** (algoritmo ETKDGv3 e campi di forza MMFF94/UFF) con l'accuratezza ab-initio/DFT di **Psi4**.

Include un'interfaccia grafica moderna in **Streamlit** che consente di visualizzare e sovrapporre conformeri in 3D, misurare parametri geometrici (senza dover ricorrere a software esterni come Avogadro), calcolare la curva continua delle barriere di potenziale torsionali (1D PES Scan rilassata e mappe 2D) e analizzare le popolazioni termodinamiche di Boltzmann e le proprietà elettroniche molecolari.

---

## Caratteristiche Principali

- **Generazione Stocastica Multithread**: Campionamento conformazionale accelerato con algoritmo ETKDGv3 (`rdkit`).
- **Supporto Campi di Forza Avanzati**: Pre-rilassamento geometrico con **MMFF94** (o MMFF94s), con fallback automatico su **UFF**.
- **Screening Energetico a Imbuto (Energy Window)**: Scarto preventivo delle conformazioni ad alta tensione classica per abbattere i tempi computazionali quantistici.
- **Deduplicazione RMSD Doppia (Pre & Post-DFT)**: Eliminazione rigorosa delle ridondanze geometriche classiche e post-ottimizzazione quantistica dei conformeri collassati nello stesso minimo locale.
- **Riconoscimento Automatico di Carica e Molteplicità**: Rilevamento dello stato di carica formale e del numero di elettroni spaiati direttamente dalla struttura molecolare.
- **Calcolo Quantomeccanico con Psi4**: Supporto per DFT ibrida (`B3LYP`, funzionali con dispersione empirica come `wB97X-D`, `PBE`, ecc.) e ab-initio (`HF`, `MP2`).
- **Modello di Solvatazione Implicita (IEFPCM)**: Simulazione in fase liquida con dielettrico continuo per acqua e solventi organici (cloroformio, DMSO, metanolo, toluene, acetonitrile, DCM).
- **Termochimica Completa & Popolazioni di Gibbs**: Calcolo opzionale delle frequenze armoniche per ottenere ZPVE, Entropia ($S$) e l'energia libera di Gibbs reale ($\Delta G^\circ_{298}$), utilizzata per la corretta statistica di Boltzmann.
- **Proprietà Elettroniche e Polarità**: Estrazione automatica della funzione d'onda SCF del **momento di dipolo elettrico ($\mu$ in Debye)**, delle energie orbitali **HOMO e LUMO** e del relativo **gap energetico** ($\Delta E_{gap}$).
- **Scansione Torsionale 1D Rilassata (Relaxed PES Scan)**: Calcolo continuo della superficie di energia potenziale da 0° a 360° attorno a qualsiasi angolo diedro con rilassamento geometrico vincolato per eliminare picchi sterici artificiali.
- **Mappe di Potenziale Torsionale 2D (Ramachandran-like)**: Mappe di contorno 2D e superfici 3D interattive per esplorare coppie di rotori accoppiati.
- **Ricerca Nomi PubChem Integrata**: Risoluzione istantanea da nome comune (es. *aspirin*, *paracetamol*) in SMILES canonico tramite le API REST di PubChem.
- **Sovrapposizione 3D (Conformer Overlay)**: Allineamento geometrico tramite algoritmo di rotazione di Kabsch (SVD) per sovrapporre qualsiasi conformero eccitato al ground state e visualizzare a colpo d'occhio i legami ruotati.
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

---

## Utilizzo

### Interfaccia Grafica Web (Streamlit Studio)
Puoi avviare l'applicazione grafica completa semplicemente facendo doppio clic su:
`avvia_gui.bat`

Oppure da terminale:
```bash
streamlit run app_gui.py
```

L'interfaccia è strutturata in 5 ambienti specializzati:
1. **Struttura 3D e Geometria**:
   - Visualizzatore molecolare WebGL interattivo (rotazione, zoom e traslazione).
   - Modalità **Sovrapponi al Minimo Globale (Overlay)** con allineamento Kabsch.
   - **Ispettore Geometrico (Anti-Avogadro)**: misurazione istantanea di distanze di legame (Å), angoli di legame (°) e angoli diedri/torsionali (°).
   - Lettura del momento di dipolo ($\mu$), HOMO, LUMO, gap energetico e proprietà di Gibbs del conformero attivo.
2. **Scansione Torsionale (PES 1D e Mappa 2D)**:
   - Scansione 1D con opzione rilassata (vincolata) o rigida con slider continuo e molecola 3D animata.
   - Mappa conformazionale 2D (contour plot e superficie 3D) per coppie di diedri accoppiati.
3. **Analisi Energetica e Popolazioni**:
   - Spettro $\Delta E$ e $\Delta G^\circ$, diagramma a ciambella delle popolazioni di Boltzmann a 298.15 K e polarità dipolare.
4. **Tabella Dati ed Esportazione**:
   - Tabella riassuntiva di tutti i conformeri con dati quantistici, geometrici ed elettronici.
   - Download del report Excel (`.xlsx`), file coordinate `.xyz` e apertura rapida della cartella.
5. **Log di Calcolo**:
   - Visualizzatore del registro di calcolo (`run.log`) e dell'output scientifico di Psi4 (`psi4_output.dat`).

---

## Opzioni della Riga di Comando

| Parametro | Breve | Default | Descrizione |
| :--- | :--- | :--- | :--- |
| `smiles` | Posizionale | `None` | Stringa SMILES della molecola bersaglio |
| `--input-file` | `-f` | `None` | File con elenco di stringhe SMILES (batch) |
| `--n-conformers` | `-n` | `50` | Numero di conformeri stocastici iniziali |
| `--rmsd-threshold`| `-r` | `0.2` | Soglia RMSD classica (Å) per deduplicazione pre-DFT |
| `--post-dft-rmsd` | | `0.15` | Soglia RMSD (Å) post-DFT per deduplicare conformeri collassati |
| `--ff` | | `mmff94` | Campo di forze: `mmff94`, `mmff94s`, `uff` |
| `--energy-window` | `-w` | `10.0` | Finestra energetica FF in kcal/mol (0 per disattivare) |
| `--top-n` | | `None` | Limite massimo di conformeri unici per Psi4 |
| `--method` | `-m` | `b3lyp` | Metodo quantistico (`b3lyp`, `wb97x-d`, `pbe`, ecc.) |
| `--basis` | `-b` | `6-31G(d)` | Basis set quantistico per Psi4 |
| `--solvent` | | `None` | Solvente implicito PCM (`water`, `chloroform`, `dmso`, ecc.) |
| `--frequencies` | | `False` | Calcola frequenze vibrazionali e $\Delta G^\circ$ di Gibbs |
| `--threads` | `-t` | `8` | Thread CPU per calcoli paralleli |
| `--memory` | | `'8 GB'` | Memoria RAM allocata a Psi4 |
| `--sp-only` | | `False` | Esegue solo Single-Point senza ottimizzazione geometrica |
| `--charge` | | `Auto` | Carica formale forzata |
| `--multiplicity` | | `Auto` | Molteplicità di spin forzata |
| `--temp` | | `298.15` | Temperatura in K per la statistica di Boltzmann |
| `--output-dir` | `-o` | `output/` | Directory personalizzata di salvataggio |
