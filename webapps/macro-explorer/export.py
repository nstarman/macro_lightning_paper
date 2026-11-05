"""Export the macro-constraint explorer's dataset.

    python webapps/macro_explorer/export.py --out <dir>

Grids are the macrolimits f_DM maps (Picker+2026) from src/static/constraint_grids.npz, downsampled
4x in each axis. Straight lightning is recomputed in the browser from <v_rel; sigma>, so its
exposure can be varied; `check_lightning` is the Python result at the default exposure, for the test.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np

import figpress_apps

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "scripts"))
import lightning as L  # noqa: E402
from constraint_grids import CONSTRAINTS, load  # noqa: E402

DATASET = "macro-constraints"
STEP = 4  # downsampling stride on the 1000 x 999 macrolimits grid


def main(out: Path) -> Path:
    g = load()
    jm, js = np.arange(0, len(g["marray"]), STEP), np.arange(0, len(g["sarray"]), STEP)
    m, s = g["marray"][jm], g["sarray"][js]
    arrays = {f"f_{key}": g[key][np.ix_(js, jm)].astype(np.float32) for key, *_ in CONSTRAINTS}
    arrays |= {
        "log10m": np.log10(m),  # GeV
        "log10s": np.log10(s),  # cm^2
        **{f"vbar_{pl}": L.mean_speed_above(s, planet=pl) for pl in L.PLANETS},  # cm/s
        **{f"check_lightning_{pl}": L.lightning(m, s, planet=pl).astype(np.float32) for pl in L.PLANETS},
        "viridis": (matplotlib.colormaps["viridis"](np.linspace(0, 1, 256)) * 255).astype(np.uint8),
    }
    meta = {
        "title": "Macroscopic dark matter constraints",
        # label anchor (log10 m, log10 sigma) and rotation [deg], as on the slide frames
        "constraints": [{"key": k, "label": lab, "at": [float(np.log10(lm_)), float(np.log10(ls_)), rot]}
                        for k, _, lab, (lm_, ls_, rot) in CONSTRAINTS],
        "lightning": {"rho_g": L.RHO_DM, "n_obs": L.N_OBS, "sig_over_m_max": L.SIG_OVER_M_MAX,
                      "gev2g": L.GEV2G, "years_default": L.T_EXPOSURE / 3.156e7, "year_s": 3.156e7,
                      "planets": {pl: {"r": p["r"], "f_ts": p["f_ts"]} for pl, p in L.PLANETS.items()}},
        "f_range": [-15, float(np.log10(1.1))],
        "rho_default": 2.8e14,  # g/cm^3, macrolimits' nuclear density
        "source": "https://arxiv.org/abs/2609.05626",
    }
    return figpress_apps.export(
        out, DATASET, arrays, meta=meta,
        units={"log10m": "log10 GeV", "log10s": "log10 cm^2", "vbar_earth": "cm/s", "vbar_jupiter": "cm/s"},
        descriptions={f"vbar_{pl}": f"<v_rel; sigma> at {pl} for macros fast enough to seed straight lightning"
                      for pl in L.PLANETS} | {f"check_lightning_{pl}": f"Python f_DM, {pl} lightning, default 2-yr exposure"
                                               for pl in L.PLANETS},
        licence="MIT", overwrite=True,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("build/data"))
    print(main(ap.parse_args().out))
