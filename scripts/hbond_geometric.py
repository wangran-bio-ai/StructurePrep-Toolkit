"""Standalone geometric H-bond counter: PEB N/O vs protein N/O within 3.5 A.
Usage: python hbond_geometric.py <VARIANT>
"""
import sys, os, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = sys.argv[1] if len(sys.argv) > 1 else "Y262A"
BASE = "/root/autodl-tmp/GROMACS"
GMX = "/root/gromacs_cuda/bin/gmx"
if not os.path.exists(GMX):
    import glob
    GMX = glob.glob("/root/*/gromacs_cuda/bin/gmx")[0]
D = f"{BASE}/{V}/out_biop_{V}" if V != "WT" else f"{BASE}/WT/out_biop"
OUT = f"{BASE}/ANALYSIS/{V}"
os.makedirs(f"{OUT}/xvg", exist_ok=True)
os.makedirs(f"{OUT}/figs", exist_ok=True)
os.makedirs(f"{OUT}/logs", exist_ok=True)

TPR = f"{D}/md.tpr"
XTC = sorted([f for f in os.listdir(D) if f.endswith(".xtc") and f.startswith("md_")])
XTC = f"{D}/{XTC[0]}"


def names_from(path):
    out, sec = [], None
    for ln in open(path, errors="ignore"):
        s = ln.strip()
        if s.startswith("["):
            sec = s.strip("[] \t"); continue
        if sec == "atoms" and s and not s.startswith(";"):
            f = s.split()
            if len(f) >= 6 and f[0].isdigit():
                out.append(f[4])
    return out


pn = names_from(f"{D}/topol.top")
qn = names_from(f"{D}/PEPB.acpype/PEPB_GMX.itp")
print(f"protein atoms {len(pn)}, PEB atoms {len(qn)}")

prot_idx = [i + 1 for i, n in enumerate(pn) if n[0] in "NO"]
peb_idx = [len(pn) + i + 1 for i, n in enumerate(qn) if n[0] in "NO"]
ndx = f"{OUT}/polar.ndx"
with open(ndx, "w") as f:
    f.write("[ Prot_polar ]\n" + " ".join(map(str, prot_idx)) + "\n")
    f.write("[ PEB_polar ]\n" + " ".join(map(str, peb_idx)) + "\n")
print(f"polar atoms: protein {len(prot_idx)}, PEB {len(peb_idx)}")

cmd = [GMX, "mindist", "-s", TPR, "-f", XTC, "-n", ndx, "-d", "0.35", "-dt", "0.1",
       "-od", f"{OUT}/xvg/polar_mindist.xvg", "-on", f"{OUT}/xvg/polar_pairs.xvg"]
r = subprocess.run(cmd, input="0\n1\n", text=True, capture_output=True)
open(f"{OUT}/logs/polar.log", "w").write(r.stdout + r.stderr)
print("rc =", r.returncode)
if r.returncode != 0:
    print((r.stdout + r.stderr)[-1200:])
    sys.exit(1)

d = np.loadtxt(f"{OUT}/xvg/polar_pairs.xvg", comments=["#", "@"])
t, n = d[:, 0], d[:, 1]
print("candidate H-bonds: mean %.1f | median %.1f | max %d" % (n.mean(), np.median(n), n.max()))
fig, ax = plt.subplots(figsize=(5.5, 3.2))
ax.plot(t, n, lw=.6, color="C6")
ax.set(xlabel="time (ns)", ylabel="# polar pairs < 3.5 A",
       title=f"candidate H-bonds (geometric) - {V}")
plt.tight_layout()
plt.savefig(f"{OUT}/figs/05_hbonds.png", dpi=150)
print("saved", f"{OUT}/figs/05_hbonds.png")
