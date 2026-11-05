"""Export the strike-count app's dataset: the counts it draws and its two confidence levels.

    python webapps/strike-count/export.py --out <dir>

The app computes everything else (Poisson probabilities, the limit, random arrivals) in the browser.
"""

import argparse
from pathlib import Path

import numpy as np

import figpress_apps

DATASET = "strike-count"


def main(out: Path) -> Path:
    return figpress_apps.export(
        out, DATASET, {"n": np.arange(0, 13, dtype=np.int32)},
        meta={
            "title": "Seeing no straight lightning",
            "levels": [0.95, 0.9],  # the paper's (Sec. IV), and Picker+2026's (the talk's)
            "lam_default": 1.0,
        },
        descriptions={"n": "strike counts drawn in the histogram"},
        licence="MIT", overwrite=True,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("build/data"))
    print(main(ap.parse_args().out))
