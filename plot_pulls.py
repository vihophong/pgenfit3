#!/usr/bin/env python3
"""Pull distributions from a pgenfit3 pull study.

Reads the PULLDATA lines of run_pull_study.sh's <prefix>_all.txt. Every
parameter that has a "<name>true" field (l0 always; any floated species
parameter) gets pull = (fit - true)/sigma, with sigma the Minos error on the
side facing the true value (errhi if fit < true, |errlo| otherwise).
Prints mean and width of each pull distribution with their statistical
uncertainties and saves a histogram figure.

Usage: ./plot_pulls.py <prefix>_all.txt [out.png] [--strict-status] [--min-covqual N]
"""
import argparse
import math
import re
import sys

NUM = r"-?(?:nan|inf|[\d.]+(?:[eE][-+]?\d+)?)"


def fields(line):
    return {k: float(v) for k, v in re.findall(rf"(?<![\w])(\w+)=({NUM})", line)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("outfile", nargs="?")
    ap.add_argument("--strict-status", action="store_true", help="require status==0 (default: status>=0)")
    ap.add_argument("--min-covqual", type=int, default=2)
    args = ap.parse_args()

    lines = [l for l in open(args.infile) if "PULLDATA l0=" in l]
    nfail = sum(1 for l in open(args.infile) if l.strip() and "PULLDATA l0=" not in l)
    if not lines:
        sys.exit("no PULLDATA lines in " + args.infile)
    params = sorted({k[:-4] for k in fields(lines[0]) if k.endswith("true")}, key=lambda p: (p != "l0", p))

    pulls = {p: [] for p in params}
    nconv = 0
    for l in lines:
        f = fields(l)
        ok = (f["status"] == 0 if args.strict_status else f["status"] >= 0) and f["covQual"] >= args.min_covqual
        if not ok:
            continue
        nconv += 1
        for p in params:
            v, t = f.get(p), f.get(p + "true")
            hi, lo = f.get(p + "errhi"), f.get(p + "errlo")
            if None in (v, t, hi, lo):
                continue
            sigma = hi if v < t else abs(lo)
            if sigma > 0 and math.isfinite(sigma):
                pulls[p].append((v - t) / sigma)

    print(f"{len(lines)} fitted toys, {nfail} failed, {nconv} converged (status/covQual cut)")
    summary = {}
    for p in params:
        x = pulls[p]
        n = len(x)
        if n < 2:
            continue
        mean = sum(x) / n
        sd = math.sqrt(sum((a - mean) ** 2 for a in x) / (n - 1))
        summary[p] = (n, mean, sd / math.sqrt(n), sd, sd / math.sqrt(2 * (n - 1)))
        print(f"{p:>8}: n={n}  mean={mean:+.4f} +/- {sd/math.sqrt(n):.4f}  "
              f"({mean/(sd/math.sqrt(n)):+.2f} sigma)  width={sd:.4f} +/- {sd/math.sqrt(2*(n-1)):.4f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    out = args.outfile or re.sub(r"(_all)?\.txt$", "", args.infile) + "_pulls.png"
    ps = [p for p in params if p in summary]
    fig, axes = plt.subplots(1, len(ps), figsize=(4.5 * len(ps), 3.8), squeeze=False)
    for ax, p in zip(axes[0], ps):
        n, mean, emean, sd, esd = summary[p]
        ax.hist(pulls[p], bins=50, range=(-5, 5), histtype="step", color="k")
        xs = [-5 + 10 * i / 400 for i in range(401)]
        w = 10 / 50 * n
        ax.plot(xs, [w / (sd * math.sqrt(2 * math.pi)) * math.exp(-0.5 * ((u - mean) / sd) ** 2) for u in xs], "r-")
        ax.set_title(f"{p} pull")
        ax.set_xlabel("(fit - true) / sigma_MINOS")
        ax.text(0.03, 0.95, f"n={n}\nmean {mean:+.3f}±{emean:.3f}\nwidth {sd:.3f}±{esd:.3f}",
                transform=ax.transAxes, va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print("saved", out)


if __name__ == "__main__":
    main()
