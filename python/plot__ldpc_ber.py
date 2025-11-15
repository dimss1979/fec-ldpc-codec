"""
plot_ldpc_ber.py
===========================================================
BER Curves for LDPC (SPA Decoder) over AWGN (BPSK)
-----------------------------------------------------------
Automatically loads:
    results/ldpc_ber_N{N}_wc{wc}_wr{wr}_iter{iter}_data.csv

Output:
    images/ldpc_ber_graph.png
    images/ldpc_ber_graph.svg

Notes:
    - Academic-quality plotting (Times New Roman + stix)
    - Auto-extracts (N, wc, wr, iter) from filename
    - Auto-displays N, wc, wr, R, iter inside the graph
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import re


# =====================================================================
#  Find latest CSV file matching pattern:
#      ldpc_ber_N{N}_wc{wc}_wr{wr}_iter{iter}_data.csv
# =====================================================================
def find_latest_ldpc_csv():
    results_dir = "results"
    if not os.path.exists(results_dir):
        raise FileNotFoundError("results/ directory not found.")

    pattern = re.compile(r"ldpc_ber_N(\d+)_wc(\d+)_wr(\d+)_iter(\d+)_data\.csv")

    candidates = []
    for f in os.listdir(results_dir):
        m = pattern.match(f)
        if m:
            candidates.append((f, m))

    if not candidates:
        raise FileNotFoundError(
            "No file matching ldpc_ber_N*_wc*_wr*_iter*_data.csv found."
        )

    # 最新更新ファイルを選択
    candidates.sort(key=lambda x: os.path.getmtime(os.path.join(results_dir, x[0])))
    filename, match = candidates[-1]
    return os.path.join(results_dir, filename), match


# =====================================================================
#  Load CSV (auto-detected)
# =====================================================================
csv_path, meta = find_latest_ldpc_csv()
N = int(meta.group(1))
wc = int(meta.group(2))
wr = int(meta.group(3))
iter_spa = int(meta.group(4))

R = (N - (N * wc) / wr) / N  # systematic LDPC: M = N*wc/wr, K=N-M

print(f"Loaded file: {csv_path}")
print(f"N={N}, wc={wc}, wr={wr}, iter={iter_spa}, R={R:.4f}")

df = pd.read_csv(csv_path)

EbN0 = df["EbN0_dB"]
ber_info = df["BER_info"]
ber_bpsk = df["BER_bpsk"]

# =====================================================================
#  Output directory
# =====================================================================
os.makedirs("images", exist_ok=True)

# =====================================================================
#  Matplotlib settings
# =====================================================================
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 14

# =====================================================================
#  Create BER figure
# =====================================================================
plt.figure(figsize=(7.5, 6))

# --- LDPC (SPA) BER curve ---
plt.semilogy(
    EbN0,
    ber_info,
    marker="o",
    markersize=8,
    markerfacecolor="none",
    markeredgewidth=1.8,
    linewidth=2.5,
    color="green",
    label="LDPC SPA BPSK",
)

# --- Uncoded BPSK theory ---
plt.semilogy(
    EbN0,
    ber_bpsk,
    linewidth=2.5,
    color="red",
    label="Uncoded BPSK (theory)",
)

# Axes
plt.ylim(1e-5, 1)
plt.xlabel("Eb/N0 [dB]", fontsize=18)
plt.ylabel("Bit Error Rate (BER)", fontsize=18)
plt.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)

# Legend
plt.legend(fontsize=14, loc="upper right", frameon=True, edgecolor="black")

# =====================================================================
#  In-graph annotation (N, wc, wr, R, iter)
# =====================================================================
text_info = (
    f"LDPC Parameters:\n"
    f"N = {N}\n"
    f"wc = {wc}, wr = {wr}\n"
    f"Rate R ≈ {R:.4f}\n"
    f"SPA iterations = {iter_spa}"
)

plt.text(
    0.03,
    0.03,
    text_info,
    fontsize=13,
    transform=plt.gca().transAxes,
    verticalalignment="bottom",
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="black"),
)

plt.tight_layout()

# =====================================================================
#  Save outputs
# =====================================================================
save_png = f"images/ldpc_ber_N{N}_wc{wc}_wr{wr}_iter{iter_spa}.png"
save_svg = f"images/ldpc_ber_N{N}_wc{wc}_wr{wr}_iter{iter_spa}.svg"

plt.savefig(save_png, dpi=300, bbox_inches="tight")
plt.savefig(save_svg, bbox_inches="tight")

print(f"Saved PNG: {save_png}")
print(f"Saved SVG: {save_svg}")

plt.show()
