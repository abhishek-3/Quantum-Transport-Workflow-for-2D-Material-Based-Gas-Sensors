# Transport Calculations — TranSIESTA/TBtrans 4.1-b4

A complete, reusable workflow for computing quantum transport properties
of 2D systems using SIESTA, TranSIESTA, and TBtrans, with
post-processing via sisl.

## Repository Structure

```
Transport/
├── electrode/
│   └── electrode.fdf          # SIESTA electrode SCF input
├── scattering/
│   ├── scattering.fdf         # TranSIESTA scattering region input
│   └── tbtrans.fdf            # TBtrans post-processing input
├── analysis/
│   └── extract_transport.py   # sisl post-processing and plotting
└── README.md
```

---

## Requirements

### Software
| Package | Version | Purpose |
|---|---|---|
| SIESTA | 4.1-b4 | Electrode SCF calculation |
| TranSIESTA | 4.1-b4 | NEGF transport calculation |
| TBtrans | 4.1-b4 | Transmission post-processing |
| Python | ≥ 3.8 | Analysis scripts |
| sisl | ≥ 0.12 | TBT.nc file parsing |
| numpy | any | Numerical operations |
| matplotlib | any | Plotting |
| scipy | any | Integration utilities |

### Python Installation
```bash
pip install sisl numpy matplotlib scipy
```

---

## Workflow Overview

```
Step 1: Electrode SCF
    electrode.fdf  →  siesta  →  electrode.TSHS

Step 2: Scattering Region (per bias point)
    scattering.fdf  →  transiesta  →  scattering.TSHS
                                    →  scattering.TSDE
                                    →  scattering.DM

Step 3: Transmission Post-processing (per bias point)
    tbtrans.fdf  →  tbtrans  →  scattering.TBT.nc
                              →  scattering.AVTRANS_Left_Right

Step 4: Analysis and Plotting
    extract_transport.py  →  figures/
```

---

## Step-by-Step Instructions

### Step 1 — Electrode Calculation

The electrode is the bulk-periodic unit cell of your lead material.
It must have **no vacuum** in the transport direction, and atoms must
interact only with their nearest periodic image (one principal layer).

**Edit `electrode/electrode.fdf`:**

1. Set `NumberOfSpecies` and `NumberOfAtoms` for your system.
2. Update `ChemicalSpeciesLabel` with your species.
3. Replace `LatticeVectors` with your electrode unit cell.
   - Transport direction: **no vacuum**
   - Cell length in transport direction ≈ one principal layer thickness
4. Replace `AtomicCoordinatesAndAtomicSpecies` with your coordinates.
5. Set the k-grid:
   - Transport direction: dense (9 or more)
   - Lateral direction: moderate (9)
   - Vacuum direction: always 1

**Run:**
```bash
cd electrode/
mpirun -np <N> siesta electrode.fdf > electrode.out
```

**Check output:**
```bash
grep "SCF Converged" electrode.out        # confirms SCF converged
ls -lh electrode.TSHS                    # must exist before next step
```

**Copy to scattering directory:**
```bash
cp electrode.TSHS ../scattering/
```

---

### Step 2 — Scattering Region Calculation

The scattering cell contains:
```
[ left electrode buffer | device region | right electrode buffer ]
```

The left and right electrode buffer atoms must be **geometrically
identical** to the electrode unit cell (same fractional coordinates,
same atom ordering). The device region sits between them and may
contain any adsorbed molecules or defects.

**Verify electrode-scattering consistency** before running:
- Left electrode atoms (1 to `used-atoms`) must match `electrode.fdf`
  coordinates exactly (within 0.001 Å).
- Right electrode atoms (last `used-atoms` atoms) must match the
  electrode after subtracting the X offset.
- Lattice Y and Z vectors must be identical in both files.

**Edit `scattering/scattering.fdf`:**

1. Set `NumberOfAtoms` (left buffer + device + right buffer).
2. Update `ChemicalSpeciesLabel` — must match `electrode.fdf`.
3. Set `used-atoms` in both electrode blocks to the electrode atom count.
4. Set `semi-inf-dir` to match your transport direction (`-A`/`+A` etc.).
5. Replace `LatticeVectors` with your scattering cell vectors.
   - Transport direction = N × electrode transport length
6. Replace `AtomicCoordinatesAndAtomicSpecies` with your coordinates.
7. Set `TS.Voltage = 0.0 eV` for the first run.

**Run at zero bias (first run):**
```bash
cd scattering/
mpirun -np <N> transiesta scattering.fdf > scattering_V0.0.out
```

**Check output:**
```bash
grep "SCF Converged\|transiesta: SCF" scattering_V0.0.out
ls -lh scattering.TSHS scattering.TSDE scattering.DM
```

---

### Step 2b — Bias Chaining (finite bias)

For each subsequent bias point, create a new directory, copy the
previous `.DM` and `.TSDE` files, and update `TS.Voltage`.

**Recommended directory structure:**
```
scattering/
├── bias_0.0eV/
│   ├── scattering.fdf    (TS.Voltage 0.0 eV, TS.DMUpdate.Init true)
│   ├── electrode.TSHS
│   └── scattering.TSHS / .TSDE / .DM   ← written by TranSIESTA
├── bias_0.1eV/
│   ├── scattering.fdf    (TS.Voltage 0.1 eV, TS.DE.Read true)
│   ├── electrode.TSHS
│   ├── scattering.DM     ← copied from bias_0.0eV/
│   ├── scattering.TSDE   ← copied from bias_0.0eV/
│   └── scattering.TSHS / .TSDE / .DM   ← written by TranSIESTA
...
```

**Flags to change between bias points:**

| Flag | V = 0.0 eV | V > 0.0 eV |
|---|---|---|
| `TS.Voltage` | `0.0 eV` | `0.1 eV`, `0.2 eV`, ... |
| `TS.DMUpdate.Init` | `true` | `false` |
| `DM.UseSaveDM` | `false` | `true` |
| `TS.DE.Read` | `false` | `true` |
| `TS.DE.Save` | `true` | `true` |
| `chem-pot` in electrode blocks | absent | add `chem-pot Left/Right` |
| `%block TS.ChemPots` | commented out | uncomment |

**Example bash loop for bias scan:**
```bash
for V in 0.0 0.1 0.2 0.3 0.4 0.5; do
    mkdir -p bias_${V}eV
    cp electrode.TSHS bias_${V}eV/
    # copy .DM and .TSDE from previous step if V > 0
    cd bias_${V}eV
    sed "s/TS.Voltage.*$/TS.Voltage ${V} eV/" ../scattering.fdf > scattering.fdf
    mpirun -np <N> transiesta scattering.fdf > scattering_V${V}.out
    cd ..
done
```

---

### Step 3 — TBtrans Post-processing

Run TBtrans in the **same directory** as the corresponding TranSIESTA
output (electrode.TSHS and scattering.TSHS must both be present).

**Edit `scattering/tbtrans.fdf`:**

1. Set `TBT.Voltage` to match the TranSIESTA bias point.
2. At zero bias: ChemPots blocks are not needed (commented out).
3. At finite bias: uncomment ChemPots blocks and add `chem-pot`
   lines to electrode blocks.
4. Set `TBT.Atoms.Device` to the device atom range (1-indexed).
5. Set `TBT.T.Eig` to the number of eigenchannels to compute.

**Run:**
```bash
cd bias_0.0eV/
cp ../tbtrans.fdf .
mpirun -np <N> tbtrans tbtrans.fdf > tbtrans.out
```

**Optional — run BTD analyze step first (recommended):**
```bash
tbtrans -fdf TBT.Analyze tbtrans.fdf > analyze.out
# Read analyze.out and pick the best pivot method
# Then set TBT.BTD.Pivot.Device in tbtrans.fdf
```

**Check output:**
```bash
ls -lh scattering.TBT.nc                  # main output
ls scattering.AVTRANS_Left_Right           # plain-text transmission
```

---

### Step 4 — Analysis and Plotting

**Edit `analysis/extract_transport.py` USER CONFIGURATION section:**

```python
SYSTEM_LABEL       = "scattering"      # match SystemLabel in fdf files
ELEC_LEFT          = "Left"            # match TBT.Elecs names
ELEC_RIGHT         = "Right"
DEVICE_ATOM_START  = 120               # 0-indexed device start atom
DEVICE_ATOM_END    = 240               # 0-indexed device end atom
VOLTAGES           = [0.0, 0.1, ...]   # all computed bias points
BIAS_DIR_TEMPLATE  = "bias_{:.1f}eV"  # directory naming convention
ZERO_BIAS_DIR      = "bias_0.0eV"     # zero-bias directory
E_MIN, E_MAX       = -3.0, 3.0        # energy range for plots
N_EIGCHANNELS      = 4                 # number of eigenchannels
```

**Run:**
```bash
cd analysis/
python extract_transport.py
```

**Output figures in `analysis/figures/`:**

| File | Description |
|---|---|
| `transmission.png` | Zero-bias T(E) spectrum |
| `eigenchannels.png` | T(E) decomposed into eigenchannels |
| `spectral_dos.png` | Spectral DOS on device atoms |
| `transmission_contour.png` | 2D T(E, V) colormap |
| `iv_rv_curves.png` | I-V and R-V curves |
| `transmission_comparison.png` | Pristine vs adsorbed T(E) |

---

## Common Errors and Fixes

| Error | Cause | Fix |
|---|---|---|
| `Electrode connectivity is not perfect` | Basis orbital range exceeds electrode cell | Double electrode cell in transport direction |
| `Semi-infinite direction not understood` | Wrong `semi-inf-dir` format | Use `-A`, `+A`, `-B`, `+B`, `-C`, `+C` |
| `Electrode coordinates do not overlap` | Atom ordering mismatch | Rebuild scattering POSCAR by tiling electrode POSCAR |
| `Principal cell extending out @R=-2` | Electrode too thin for basis cutoff | Use 2×1×1 or 3×1×1 electrode |
| `no unit specified for SCF.H.Tolerance` | Missing unit | Change to `SCF.H.Tolerance 1.0e-4 eV` or remove line |
| `Minimum split_norm parameter` warning | `PAO.SplitNorm` too small | Increase to 0.25 or 0.30 |
| `TBT.nc not found` in Python script | TBtrans did not complete | Check `tbtrans.out` for errors |

---

## Transport Direction Reference

| Transport along | `semi-inf-dir Left` | `semi-inf-dir Right` | Dense k in |
|---|---|---|---|
| X (A vector) | `-A` | `+A` | A |
| Y (B vector) | `-B` | `+B` | B |
| Z (C vector) | `-C` | `+C` | C |

---

## Citation

If you use this workflow, please cite:

- SIESTA: J. M. Soler et al., *J. Phys.: Condens. Matter* 14, 2745 (2002)
- M. Brandbyge, J.-L. Mozos, P. Ordejón, J. Taylor, and K. Stokbro, Phys. Rev. B 65, 165401 (2002).
- TranSIESTA: N. Papior et al., *Comput. Phys. Commun.* 212, 8 (2017)
- TBtrans: N. Papior et al., *Comput. Phys. Commun.* 212, 8 (2017)
- sisl: N. R. Papior, *sisl* (2018), https://doi.org/10.5281/zenodo.597181

---

## Author

Abhishek
abhishek.3@iitj.ac.in
