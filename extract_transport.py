"""
extract_transport.py
====================
Extracts and plots transport properties from TBtrans NetCDF (.TBT.nc)
output files using sisl.

Capabilities
------------
1. Transmission spectrum T(E) at a single bias point
2. Transmission eigenchannels T_n(E)
3. Spectral DOS (ADOS) on device atoms
4. I-V curve and R-V curve from multiple bias points
5. Transmission contour map T(E, V)

Requirements
------------
    pip install sisl numpy matplotlib scipy

Usage
-----
    python extract_transport.py

Edit the USER CONFIGURATION section below to match your system.

Author  : [Your Name]
System  : Mo2CO2 MXene — TranSIESTA/TBtrans 4.1-b4
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import warnings
warnings.filterwarnings('ignore')

try:
    import sisl
except ImportError:
    raise ImportError(
        "sisl is required. Install with: pip install sisl"
    )

# ============================================================
# USER CONFIGURATION — edit this section to match your system
# ============================================================

# SystemLabel used in TBtrans (must match tbtrans.fdf SystemLabel)
SYSTEM_LABEL = "scattering"

# Electrode names (must match TBT.Elecs block in tbtrans.fdf)
ELEC_LEFT  = "Left"
ELEC_RIGHT = "Right"

# Device atom range (0-indexed for sisl)
# If your device atoms are 121-240 (1-indexed), set:
DEVICE_ATOM_START = 120     # 0-indexed = atom 121
DEVICE_ATOM_END   = 240     # 0-indexed = atom 240 (exclusive in Python slicing)

# Bias voltages calculated (in Volts) — list all bias points
# Must have corresponding .TBT.nc files in their directories
VOLTAGES = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1.0]

# Directory structure for bias calculations
# Each bias point should have its TBtrans output in a subdirectory.
# Format: bias_0.0eV/, bias_0.1eV/, etc.
# Set to "." if all files are in the current directory (single bias).
BIAS_DIR_TEMPLATE = "bias_{:.1f}eV"

# Zero-bias directory (for single-point transmission plot)
ZERO_BIAS_DIR = "bias_0.0eV"

# Energy range for plots (eV) — should match TBT.Contour range
E_MIN = -3.0
E_MAX =  3.0

# Number of eigenchannels to plot
N_EIGCHANNELS = 4

# Output directory for figures
OUTPUT_DIR = "figures"

# Figure format
FIG_FORMAT = "png"   # "png", "pdf", or "svg"
FIG_DPI    = 300

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_nc_path(bias_dir):
    """Return the path to the TBT.nc file for a given bias directory."""
    return os.path.join(bias_dir, f"{SYSTEM_LABEL}.TBT.nc")


def load_tbt(bias_dir):
    """Load a TBTrans NetCDF file and return the sisl sile object."""
    path = get_nc_path(bias_dir)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"TBT.nc file not found: {path}\n"
            f"Check that TBtrans completed successfully in {bias_dir}/"
        )
    return sisl.get_sile(path)


def make_output_dir():
    """Create the output directory for figures if it does not exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Figures will be saved to: {OUTPUT_DIR}/")


# ============================================================
# PLOT 1 — Transmission spectrum T(E) at zero bias
# ============================================================

def plot_transmission(bias_dir=ZERO_BIAS_DIR, label=""):
    """
    Plot the k-averaged transmission spectrum T(E).

    Parameters
    ----------
    bias_dir : str
        Directory containing the TBT.nc file.
    label : str
        Optional label for the plot title.
    """
    print(f"\n[1] Plotting transmission spectrum from {bias_dir}/")
    tbt = load_tbt(bias_dir)

    E = tbt.E                                         # energy grid (eV)
    T = tbt.transmission(ELEC_LEFT, ELEC_RIGHT)       # k-averaged T(E)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(E, T, color='steelblue', linewidth=1.8,
            label=f'T(E){" — " + label if label else ""}')
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8,
               label=r'$E_F$')
    ax.set_xlabel(r'$E - E_F$ (eV)', fontsize=13)
    ax.set_ylabel(r'Transmission $T(E)$', fontsize=13)
    ax.set_xlim(E_MIN, E_MAX)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11)
    ax.set_title('Zero-Bias Transmission Spectrum', fontsize=13)
    plt.tight_layout()

    fname = os.path.join(OUTPUT_DIR, f"transmission.{FIG_FORMAT}")
    plt.savefig(fname, dpi=FIG_DPI)
    print(f"   Saved: {fname}")
    plt.close()


# ============================================================
# PLOT 2 — Transmission eigenchannels T_n(E)
# ============================================================

def plot_eigenchannels(bias_dir=ZERO_BIAS_DIR, label=""):
    """
    Plot total transmission and individual eigenchannel contributions.

    Parameters
    ----------
    bias_dir : str
        Directory containing the TBT.nc file.
    label : str
        Optional label for the plot title.
    """
    print(f"\n[2] Plotting eigenchannels from {bias_dir}/")
    tbt = load_tbt(bias_dir)

    E    = tbt.E
    T    = tbt.transmission(ELEC_LEFT, ELEC_RIGHT)
    Teig = tbt.transmission_eig(ELEC_LEFT, ELEC_RIGHT)  # shape: (n_E, n_channels)

    n_ch = min(N_EIGCHANNELS, Teig.shape[1])
    colors = plt.cm.tab10(np.linspace(0, 0.5, n_ch))

    fig, ax = plt.subplots(figsize=(7, 5))

    # Total transmission
    ax.semilogy(E, T, color='black', linewidth=2.0,
                label='Total T(E)', zorder=5)

    # Individual eigenchannels
    for i in range(n_ch):
        ax.semilogy(E, Teig[:, i], linewidth=1.2,
                    color=colors[i], linestyle='--',
                    label=f'Channel {i+1}', alpha=0.85)

    ax.axvline(0, color='gray', linestyle=':', linewidth=0.8)
    ax.set_xlabel(r'$E - E_F$ (eV)', fontsize=13)
    ax.set_ylabel(r'Transmission $T_n(E)$', fontsize=13)
    ax.set_xlim(E_MIN, E_MAX)
    ax.set_ylim(1e-6, None)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_title(
        f'Transmission Eigenchannels{" — " + label if label else ""}',
        fontsize=13
    )
    plt.tight_layout()

    fname = os.path.join(OUTPUT_DIR, f"eigenchannels.{FIG_FORMAT}")
    plt.savefig(fname, dpi=FIG_DPI)
    print(f"   Saved: {fname}")
    plt.close()


# ============================================================
# PLOT 3 — Spectral DOS on device atoms
# ============================================================

def plot_spectral_dos(bias_dir=ZERO_BIAS_DIR, label=""):
    """
    Plot the spectral (electrode-projected) DOS on the device region.

    Parameters
    ----------
    bias_dir : str
        Directory containing the TBT.nc file.
    label : str
        Optional label for the plot title.
    """
    print(f"\n[3] Plotting spectral DOS from {bias_dir}/")
    tbt = load_tbt(bias_dir)

    E = tbt.E
    device_atoms = range(DEVICE_ATOM_START, DEVICE_ATOM_END)

    ADOS_L = tbt.ADOS(ELEC_LEFT,  atoms=device_atoms)
    ADOS_R = tbt.ADOS(ELEC_RIGHT, atoms=device_atoms)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(E, ADOS_L, color='steelblue', linewidth=1.5,
            label=f'ADOS from {ELEC_LEFT}')
    ax.plot(E, ADOS_R, color='tomato',    linewidth=1.5,
            label=f'ADOS from {ELEC_RIGHT}', linestyle='--')
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8)
    ax.set_xlabel(r'$E - E_F$ (eV)', fontsize=13)
    ax.set_ylabel(r'Spectral DOS (eV$^{-1}$)', fontsize=13)
    ax.set_xlim(E_MIN, E_MAX)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11)
    ax.set_title(
        f'Spectral DOS — Device Region{" — " + label if label else ""}',
        fontsize=13
    )
    plt.tight_layout()

    fname = os.path.join(OUTPUT_DIR, f"spectral_dos.{FIG_FORMAT}")
    plt.savefig(fname, dpi=FIG_DPI)
    print(f"   Saved: {fname}")
    plt.close()


# ============================================================
# PLOT 4 — Transmission contour map T(E, V)
# ============================================================

def plot_transmission_contour():
    """
    Plot a 2D colormap of T(E) vs bias voltage V.
    Requires TBT.nc files for all voltages in VOLTAGES list.
    """
    print(f"\n[4] Plotting transmission contour map T(E, V)")
    T_all = []
    E     = None
    valid_voltages = []

    for V in VOLTAGES:
        bias_dir = BIAS_DIR_TEMPLATE.format(V)
        try:
            tbt = load_tbt(bias_dir)
            if E is None:
                E = tbt.E
            T_all.append(tbt.transmission(ELEC_LEFT, ELEC_RIGHT))
            valid_voltages.append(V)
        except FileNotFoundError as e:
            print(f"   WARNING: skipping V={V:.1f} eV — {e}")

    if len(valid_voltages) < 2:
        print("   ERROR: need at least 2 bias points for contour map. Skipping.")
        return

    T_all = np.array(T_all)                  # shape: (n_V, n_E)
    V_arr = np.array(valid_voltages)

    fig, ax = plt.subplots(figsize=(7, 7))

    im = ax.pcolormesh(
        E, V_arr, T_all,
        cmap='viridis',
        shading='auto',
        vmin=0,
        vmax=np.percentile(T_all, 98)        # clip colorbar at 98th percentile
    )

    # Bias window triangle
    ax.plot( V_arr / 2, V_arr, 'k-', linewidth=1.2, label=r'$\pm V/2$ window')
    ax.plot(-V_arr / 2, V_arr, 'k-', linewidth=1.2)
    ax.axvline(0, color='white', linestyle='--', linewidth=0.8)

    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(r'$T(E)$', fontsize=12)

    ax.set_xlabel(r'$E - E_F$ (eV)', fontsize=13)
    ax.set_ylabel('Bias Voltage (V)', fontsize=13)
    ax.set_xlim(E_MIN, E_MAX)
    ax.set_ylim(V_arr.min(), V_arr.max())
    ax.legend(fontsize=10, loc='upper left')
    ax.set_title('Transmission Contour Map', fontsize=13)
    plt.tight_layout()

    fname = os.path.join(OUTPUT_DIR, f"transmission_contour.{FIG_FORMAT}")
    plt.savefig(fname, dpi=FIG_DPI)
    print(f"   Saved: {fname}")
    plt.close()


# ============================================================
# PLOT 5 — I-V curve and R-V curve
# ============================================================

def plot_iv_curve():
    """
    Extract current at each bias point using sisl and plot I-V and R-V curves.
    Current is computed by integrating T(E) over the bias window using the
    Landauer-Büttiker formula (sisl handles this internally).
    """
    print(f"\n[5] Plotting I-V and R-V curves")

    currents       = []
    valid_voltages = []

    for V in VOLTAGES:
        bias_dir = BIAS_DIR_TEMPLATE.format(V)
        try:
            tbt = load_tbt(bias_dir)
            # sisl.current integrates T(E) * [f_L(E) - f_R(E)] dE
            I = tbt.current(ELEC_LEFT, ELEC_RIGHT)
            currents.append(I)
            valid_voltages.append(V)
            print(f"   V = {V:.2f} V  →  I = {I*1e6:.4f} μA")
        except FileNotFoundError as e:
            print(f"   WARNING: skipping V={V:.1f} eV — {e}")
        except Exception as e:
            print(f"   WARNING: current extraction failed at V={V:.1f} eV — {e}")

    if len(valid_voltages) < 2:
        print("   ERROR: need at least 2 bias points for I-V curve. Skipping.")
        return

    V_arr = np.array(valid_voltages)
    I_arr = np.array(currents) * 1e6         # convert A → μA

    # Resistance R = V / I  (skip V=0 to avoid division by zero)
    mask  = V_arr > 0
    V_res = V_arr[mask]
    I_res = I_arr[mask]
    R_arr = V_res / (I_res * 1e-6) * 1e-3   # convert to kΩ

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- I-V panel ---
    ax = axes[0]
    ax.plot(V_arr, I_arr, 'o-', color='steelblue',
            linewidth=2.0, markersize=6, label='I-V')
    ax.set_xlabel('Bias Voltage (V)', fontsize=13)
    ax.set_ylabel(r'Current ($\mu$A)', fontsize=13)
    ax.set_title('I-V Curve', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.5)

    # --- R-V panel ---
    ax = axes[1]
    ax.plot(V_res, R_arr, 's-', color='tomato',
            linewidth=2.0, markersize=6, label='R-V')
    ax.set_xlabel('Bias Voltage (V)', fontsize=13)
    ax.set_ylabel(r'Resistance (k$\Omega$)', fontsize=13)
    ax.set_title('R-V Curve', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.5)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"iv_rv_curves.{FIG_FORMAT}")
    plt.savefig(fname, dpi=FIG_DPI)
    print(f"   Saved: {fname}")
    plt.close()


# ============================================================
# PLOT 6 — Transmission comparison: pristine vs doped/adsorbed
# ============================================================

def plot_transmission_comparison(
    dir_pristine,
    dir_doped,
    label_pristine="Pristine",
    label_doped="With adsorbate"
):
    """
    Compare transmission spectra between two systems at the same bias.
    Useful for showing how an adsorbed molecule suppresses T(E).

    Parameters
    ----------
    dir_pristine : str
        Directory containing TBT.nc for the pristine system.
    dir_doped : str
        Directory containing TBT.nc for the adsorbed/doped system.
    label_pristine : str
        Legend label for the pristine curve.
    label_doped : str
        Legend label for the doped/adsorbed curve.
    """
    print(f"\n[6] Plotting transmission comparison: {label_pristine} vs {label_doped}")

    tbt_p = load_tbt(dir_pristine)
    tbt_d = load_tbt(dir_doped)

    E   = tbt_p.E
    T_p = tbt_p.transmission(ELEC_LEFT, ELEC_RIGHT)
    T_d = tbt_d.transmission(ELEC_LEFT, ELEC_RIGHT)

    Teig_p = tbt_p.transmission_eig(ELEC_LEFT, ELEC_RIGHT)
    Teig_d = tbt_d.transmission_eig(ELEC_LEFT, ELEC_RIGHT)
    n_ch   = min(N_EIGCHANNELS, Teig_p.shape[1], Teig_d.shape[1])

    fig = plt.figure(figsize=(14, 5))
    gs  = gridspec.GridSpec(1, 2, figure=fig)

    # --- Left panel: total transmission ---
    ax1 = fig.add_subplot(gs[0])
    ax1.semilogy(E, T_p, color='steelblue', linewidth=2.0,
                 label=label_pristine)
    ax1.semilogy(E, T_d, color='tomato',    linewidth=2.0,
                 label=label_doped, linestyle='--')
    ax1.axvline(0, color='black', linestyle=':', linewidth=0.8)
    ax1.set_xlabel(r'$E - E_F$ (eV)', fontsize=13)
    ax1.set_ylabel(r'$T(E)$', fontsize=13)
    ax1.set_xlim(E_MIN, E_MAX)
    ax1.set_ylim(1e-4, None)
    ax1.legend(fontsize=11)
    ax1.set_title('Total Transmission', fontsize=13)

    # --- Right panel: eigenchannel decomposition ---
    ax2 = fig.add_subplot(gs[1])
    blue_shades = plt.cm.Blues(np.linspace(0.5, 0.9, n_ch))
    red_shades  = plt.cm.Reds( np.linspace(0.5, 0.9, n_ch))

    for i in range(n_ch):
        lbl_p = f'{label_pristine} Ch.{i+1}' if i == 0 else ''
        lbl_d = f'{label_doped} Ch.{i+1}'    if i == 0 else ''
        ax2.semilogy(E, Teig_p[:, i], color=blue_shades[i],
                     linewidth=1.2, label=lbl_p)
        ax2.semilogy(E, Teig_d[:, i], color=red_shades[i],
                     linewidth=1.2, linestyle='--', label=lbl_d)

    ax2.axvline(0, color='black', linestyle=':', linewidth=0.8)
    ax2.set_xlabel(r'$E - E_F$ (eV)', fontsize=13)
    ax2.set_ylabel(r'Eigenchannel $T_n(E)$', fontsize=13)
    ax2.set_xlim(E_MIN, E_MAX)
    ax2.set_ylim(1e-4, None)
    ax2.legend(fontsize=9)
    ax2.set_title('Eigenchannel Decomposition', fontsize=13)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"transmission_comparison.{FIG_FORMAT}")
    plt.savefig(fname, dpi=FIG_DPI)
    print(f"   Saved: {fname}")
    plt.close()


# ============================================================
# MAIN — runs all plots in sequence
# ============================================================

if __name__ == "__main__":

    make_output_dir()

    print("=" * 60)
    print("  TBtrans Post-Processing with sisl")
    print("  System Label :", SYSTEM_LABEL)
    print("  Electrodes   :", ELEC_LEFT, "→", ELEC_RIGHT)
    print("=" * 60)

    # --- Single-bias plots (zero bias by default) ---
    if os.path.exists(ZERO_BIAS_DIR):
        plot_transmission(bias_dir=ZERO_BIAS_DIR,   label="V = 0.0 eV")
        plot_eigenchannels(bias_dir=ZERO_BIAS_DIR,  label="V = 0.0 eV")
        plot_spectral_dos(bias_dir=ZERO_BIAS_DIR,   label="V = 0.0 eV")
    else:
        print(f"\nWARNING: Zero-bias directory '{ZERO_BIAS_DIR}' not found.")
        print("         Skipping single-bias plots.")

    # --- Multi-bias plots ---
    plot_transmission_contour()
    plot_iv_curve()

    # --- Comparison plot (pristine vs adsorbed) ---
    # Uncomment and set the correct directories to use this plot.
    # plot_transmission_comparison(
    #     dir_pristine = "pristine/bias_0.0eV",
    #     dir_doped    = "doped/bias_0.0eV",
    #     label_pristine = "Pristine Mo₂CO₂",
    #     label_doped    = "Mo₂CO₂ + molecule"
    # )

    print("\n" + "=" * 60)
    print("  All plots complete. Check the figures/ directory.")
    print("=" * 60)
