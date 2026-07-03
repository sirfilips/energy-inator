# Energyinator

**Conformer Generation and Quantum Chemistry Energy Calculator**

Energyinator is a Python tool for generating 3D molecular conformers from SMILES strings and calculating their energies using quantum chemistry methods. It combines RDKit for conformer generation and Psi4 for high-accuracy energy calculations.

---

## ✨ Features

- **Conformer Generation**: Creates multiple 3D conformers from SMILES using RDKit's ETKDG method
- **Force Field Optimization**: Optimizes conformers using UFF/MMFF force fields
- **Quantum Chemistry**: Calculates single-point energies and optimizes geometries using Psi4
- **RMSD Filtering**: Removes duplicate conformers based on RMSD threshold
- **Energy Analysis**: Provides energy values in Hartree, kcal/mol, and eV
- **Result Export**: Saves optimized geometries as XYZ files and energy data as Excel spreadsheets
- **Comprehensive Logging**: Detailed logging for debugging and progress tracking

---

## 📋 Requirements

### Python Dependencies

Create a virtual environment and install the required packages:

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install required packages
pip install rdkit>=2023.03.0 psi4>=1.8.0 pandas>=1.5.0 tqdm>=4.64.0
```

### System Requirements


| Requirement    | Minimum   | Recommended                    |
| -------------- | --------- | ------------------------------ |
| **Python**     | 3.8       | 3.10+                          |
| **Memory**     | 8 GB RAM  | 16 GB RAM                      |
| **CPU**        | Dual-core | Multi-core (for parallel Psi4) |
| **Disk Space** | 2 GB      | 5 GB                           |


### Psi4 System Dependencies

Psi4 requires additional system libraries. Installation methods vary by platform:

#### **Ubuntu/Debian**

```bash
sudo apt-get update
sudo apt-get install -y python3-dev cmake gfortran libopenblas-dev liblapack-dev
```

#### **CentOS/RHEL**

```bash
sudo yum install -y python3-devel cmake gcc-gfortran openblas-devel lapack-devel
```

#### **macOS (with Homebrew)**

```bash
brew install cmake gcc openblas lapack
```

#### **Windows**

Use **Windows Subsystem for Linux (WSL)** or install [MSYS2](https://www.msys2.org/):

```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gfortran mingw-w64-x86_64-cmake
```

---

## 🚀 Installation

### Option 1: Clone from GitHub

```bash
git clone https://github.com/your-username/Energinator.git
cd Energyinator
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### Option 2: Manual Setup

1. Download `Energyinator.py`
2. Create `requirements.txt` with:
  ```
   rdkit>=2023.03.0
   psi4>=1.8.0
   pandas>=1.5.0
   tqdm>=4.64.0
  ```
3. Install dependencies: `pip install -r requirements.txt`

---

## 📖 Usage

### Quick Start

```bash
# Activate virtual environment
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Run the script
python Energyinator.py
```

When prompted, enter a SMILES string. Examples:

- Water: `O`
- Methane: `C`
- Ethanol: `CCO`
- Benzene: `c1ccccc1`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`

### Programmatic Usage

```python
from Energyinator import calcola_energia_psi4

smiles = "CCO"  # Ethanol
result = calcola_energia_psi4(smiles, n_conformers=30, rmsd_threshold=0.2)

min_energy_hartree, min_energy_kcalmol, min_energy_ev, best_xyz, best_path, \
max_energy_hartree, max_energy_kcalmol, max_energy_ev, worst_xyz, worst_path, \
total_time = result

print(f"Minimum energy: {min_energy_kcalmol:.2f} kcal/mol")
print(f"Calculation time: {total_time:.2f} seconds")
```

---

## ⚙️ Configuration

Modify the global configuration at the top of `Energyinator.py`:

```python
# Memory allocation for Psi4
PSI4_MEMORY = '8 GB'

# Number of threads for parallel computation
PSI4_THREADS = 8

# Quantum chemistry method and basis set
PSI4_METHOD = 'b3lyp'
PSI4_BASIS = '6-31G(d)'
```

### Configuration Options


| Parameter        | Default      | Description                                                          |
| ---------------- | ------------ | -------------------------------------------------------------------- |
| `PSI4_MEMORY`    | `'8 GB'`     | Memory allocation for Psi4 calculations                              |
| `PSI4_THREADS`   | `8`          | Number of CPU threads for parallel computation                       |
| `PSI4_METHOD`    | `'b3lyp'`    | Quantum chemistry method (e.g., `'hf'`, `'b3lyp'`, `'mp2'`, `'pbe'`) |
| `PSI4_BASIS`     | `'6-31G(d)'` | Basis set for quantum calculations                                   |
| `n_conformers`   | `50`         | Number of conformers to generate                                     |
| `rmsd_threshold` | `0.2`        | RMSD threshold (Å) for duplicate conformer filtering                 |
| `top_n`          | `10`         | Number of unique conformers to optimize with Psi4                    |


### Recommended Settings


| Molecule Size        | Conformers | RMSD Threshold | Method | Basis Set | Memory  |
| -------------------- | ---------- | -------------- | ------ | --------- | ------- |
| Small (< 20 atoms)   | 50-100     | 0.1-0.2        | b3lyp  | 6-31G(d)  | 8 GB    |
| Medium (20-50 atoms) | 30-50      | 0.2-0.5        | pbe    | 6-31G(d)  | 8-16 GB |
| Large (> 50 atoms)   | 20-30      | 0.5-1.0        | b3lyp  | 3-21G     | 16+ GB  |


---

## 📁 Output Files

The script creates the following structure:

```
output/
└── result_YYYYMMDD_HHMMSS/
    ├── geom_ottimizzata_YYYYMMDD_HHMMSS_0.xyz
    ├── geom_ottimizzata_YYYYMMDD_HHMMSS_1.xyz
    ├── ...
    └── conformers_energies_YYYYMMDD_HHMMSS.xlsx
```

### XYZ File Format

```
6                    # Number of atoms
CCO                  # SMILES string
C 0.00000000 0.00000000 0.00000000
C 1.52000000 0.00000000 0.00000000
O 2.28000000 1.40000000 0.00000000
H 0.50000000 0.96000000 0.00000000
H 0.50000000 -0.96000000 0.00000000
H 2.28000000 -0.40000000 0.00000000
```

### Excel File Contents

- `Conformer_ID`: Unique identifier
- `Energy_Hartree`: Energy in Hartree
- `Energy_kcalmol`: Energy in kcal/mol
- `Energy_eV`: Energy in electron volts
- `UFF_Energy`: Force field energy

---

## 🔬 Supported Molecule Types

### ✅ Supported

- Organic molecules (hydrocarbons, alcohols, amines, etc.)
- Heterocyclic compounds
- Small to medium-sized molecules (up to ~50 heavy atoms)
- Neutral molecules

### ⚠️ Limitations

- **Large molecules**: May require significant computational resources
- **Transition metals**: Limited support in standard basis sets
- **Charged species**: May require different methods/basis sets
- **Radicals**: May need specialized methods (e.g., UHF, UKS)

---

## 💡 Best Practices

### Performance Tips

- **Reduce conformers** for large molecules
- **Increase RMSD threshold** to filter more duplicates
- **Use faster methods** like `pbe` instead of `b3lyp`
- **Use smaller basis sets** like `3-21G` for quick tests

### Accuracy Tips

- **Use larger basis sets** like `6-311+G(d,p)` for publication-quality results
- **Increase conformer count** for flexible molecules
- **Decrease RMSD threshold** for more diverse conformers
- **Use correlated methods** like `mp2` for higher accuracy

---

## 🐛 Troubleshooting

### 1. **Invalid SMILES**

```
ValueError: SMILES non valido o impossibile da convertire in molecola.
```

**Solution**: Verify your SMILES using [PubChem](https://pubchem.ncbi.nlm.nih.gov/).

---

### 2. **Cannot Generate 3D Conformers**

```
ValueError: Impossibile generare conformeri 3D della molecola.
```

**Solution**:

- Check molecule size/complexity
- Reduce `n_conformers`
- Ensure RDKit is properly installed

---

### 3. **Psi4 Installation Issues**

```
ModuleNotFoundError: No module named 'psi4'
```

**Solution**:

- Install system dependencies first
- Try: `conda install -c conda-forge psi4`
- Check [Psi4 installation guide](https://psicode.org/psi4manual/stable/installation.html)

---

### 4. **Memory Errors**

```
MemoryError: Not enough memory
```

**Solution**:

- Reduce `PSI4_MEMORY`
- Use smaller basis sets
- Process fewer conformers

---

### 5. **Slow Performance**

**Solution**:

- Reduce `n_conformers`
- Increase `rmsd_threshold`
- Use faster method (e.g., `pbe`)

---

## 📚 References

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [Psi4 Documentation](https://psicode.org/psi4manual/stable/)
- [SMILES Notation](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)

---

## 🤝 Contributing

Contributions are welcome!

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## 📜 License

**MIT License**

Copyright (c) 2026 Filippo Zanardi

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

## 📧 Contact

**Filippo Zanardi**

Project Link: [https://github.com/your-username/Energinator](https://github.com/your-username/Energinator)
