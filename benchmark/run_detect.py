import os
ROOT = os.environ.get("VIPAN_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")
import os, subprocess, sys, time, glob
from concurrent.futures import ProcessPoolExecutor

BIN = os.environ.get("OPENRDP_BIN") or os.path.join(ROOT, "tools", "OpenRDP", "openrdp", "bin")
S3, PT = f"{BIN}/3Seq/3seq.Unix", f"{BIN}/3Seq/myPvalueTable"
GC = f"{BIN}/GENECONV/geneconv.Unix"
SET = os.environ.get("SIM_SET", ""); D, O = f"{WORK}/data{SET}", f"{WORK}/out{SET}"

def run(name):
    fa = f"{D}/{name}.fa"
    t = {}
    # STALE_OUTPUTS_REMOVED: 3SEQ waits for an answer when an output file already exists and then spins forever
    for old in glob.glob(f"{O}/{name}.3s.*") + [f"{O}/{name}.frags", f"{O}/{name}.sum"]:
        try:
            os.remove(old)
        except OSError:
            pass
    t0 = time.time()
    subprocess.run([S3, "-full", fa, "-ptable", PT, "-d", "-id", f"{O}/{name}"],
                   input=b"Y\n" * 20, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=O, timeout=3600)
    t["3seq"] = time.time() - t0
    cfg = f"{O}/{name}.cfg"
    open(cfg, "w").write(f"#GCONV_CONFIG\n  -inputpath={fa}\n  -Indel_blocs\n  -Gscale=1\n"
                         "  -Minlength=1\n  -Minnpoly=2\n  -Minscore=2\n  -Maxoverlap=1\n")
    t0 = time.time()
    subprocess.run([GC, f"-Seqfile={fa}", f"-Config={cfg}", "-nolog"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=O)
    t["geneconv"] = time.time() - t0
    # geneconv writes <input>.frags next to the input; move it
    fr = fa.rsplit(".fa", 1)[0] + ".frags"
    if os.path.exists(fr):
        os.replace(fr, f"{O}/{name}.frags")
    for ext in (".sum",):
        x = fa.rsplit(".fa", 1)[0] + ext
        if os.path.exists(x):
            os.replace(x, f"{O}/{name}{ext}")
    return name, t

if __name__ == "__main__":
    os.makedirs(O, exist_ok=True)
    names = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{D}/rep*_*.fa"))
    with open(f"{O}/timings.tsv", "w") as f:
        f.write("set\tmethod\tseconds\n")
        with ProcessPoolExecutor(int(os.environ.get("WORKERS", os.cpu_count() or 2))) as ex:
            for name, t in ex.map(run, names):
                for k, v in t.items():
                    f.write(f"{name}\t{k}\t{v:.1f}\n")
                f.flush()
                print(name, {k: round(v, 1) for k, v in t.items()}, flush=True)
