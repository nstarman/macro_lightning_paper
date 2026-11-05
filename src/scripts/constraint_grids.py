"""Compute every macrolimits constraint (plus projected straight lightning) once, on the
macrolimits grid, and cache them to src/static/constraint_grids.npz for the explorer webapp.

The cache is committed. To recompute it (~1 min), clone https://github.com/zpicker/macrolimits
into the repository root and run `python src/scripts/constraint_grids.py`.
"""

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]  # the repository
REPO = ROOT / "macrolimits"  # git clone https://github.com/zpicker/macrolimits, only to recompute
CACHE = ROOT / "src" / "static" / "constraint_grids.npz"  # committed

# Talk order = Picker+2026 section order. (key, macrolimits method, label, label (m, sigma, rotation))
CONSTRAINTS = [
    ("hsc", "HSC", "HSC", (3e47, 3e8, 90)),
    ("wd", "WD", "White Dwarfs", (1e35, 2e-8, 0)),
    ("rg", "rg", "RG", (1e40, 1e5, 42)),
    ("fireballs", "fireballs", "Fireballs", (1e26, 3e-7, 30)),
    ("asteroids", "asteroid", "Asteroids", (2e28, 1e5, 24)),
    ("saturn", "saturn", "Saturn", (1e16, 7e-4, 28)),
    ("skylab", "skylab", "Skylab", (1e7, 2e-15, 40)),
    ("ohya", "ohya", "Ohya", (1e12, 5e-14, 40)),
    ("mica", "mica", "Mica", (1e27, 1e-16, 0)),
    ("chicago", "chicago", "Chicago", (3e0, 1e-25, 40)),
    ("deap", "deap", "DEAP3600", (1e19, 2e-21, 0)),
    ("dama", "dama", "DAMA", (2e3, 3e-28, 17)),
    ("cdms", "cdms", "CDMS", (1e3, 1e-34, 36)),
    ("xenon", "xenon", "Xenon1T", (1e3, 1e-39, 36)),
    ("gas", "gas", "Gas clouds", (1e37, 1e8, 36)),
]
# key, label, label position (m, sigma, rotation), arrow target (m, sigma) inside the region
LIGHTNING = ("lightning", "Straight lightning\n(projected)", (3e31, 1.5e-12, 0), (3e25, 1e-4))


def compute():
    sys.path.insert(0, str(Path(__file__).parent))
    from lightning import lightning

    cwd = os.getcwd()
    os.chdir(REPO)  # macrolimits loads data by relative path
    sys.path.insert(0, str(REPO))
    from macro_f_constraints import macro

    mc = macro()
    grids = {key: getattr(mc, meth)().astype(np.float32) for key, meth, *_ in CONSTRAINTS}
    os.chdir(cwd)
    grids["lightning"] = lightning(mc.marray, mc.sarray).astype(np.float32)
    CACHE.parent.mkdir(exist_ok=True)
    np.savez_compressed(CACHE, marray=mc.marray, sarray=mc.sarray, **grids)
    return CACHE


def load():
    """Cached grids: dict with marray, sarray and one f_DM array (n_sigma, n_mass) per key."""
    if not CACHE.exists():
        compute()
    return dict(np.load(CACHE))


if __name__ == "__main__":
    print("wrote", compute())
    g = load()
    for key, *_ in CONSTRAINTS + [LIGHTNING]:
        assert g[key].shape == (len(g["sarray"]), len(g["marray"])), key
        print(f"{key:10s} min f_DM = {np.nanmin(g[key]):.1e}")
