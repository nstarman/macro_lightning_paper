"""Projected straight-lightning constraint as f_DM(m, sigma), in macrolimits' format.

Physics from Starkman+2021 (PRD 103 063024; arXiv:2006.16272) and its code
(macro_lightning.physics.calculate_Mx_and_Sx), recast to f_DM following
Picker+2026 (arXiv:2609.05626):

    N_exp(m, sigma) = f_DM * (rho_DM / m) * pi R_E^2 f_TS * T * <v_rel ; sigma>

where <v_rel ; sigma> integrates v_rel f(v) over macros fast enough to seed
lightning, sigma >= max(6e-8 cm^2, K / v_rel^2). Non-observation of straight
lightning (N_obs = 2.3, 90% CL) gives f_DM <= 2.3 / N_exp(f_DM=1). Valid for
sigma/m <= 1e-4 cm^2/g (else the macro is stopped in the atmosphere).

This is a PROJECTION: it assumes T of complete monitoring with no straight bolt.
PLANETS holds Earth (the paper's Sec. IV) and Jupiter (its Sec. V: lightning over 10% of the
surface, a 5x larger K, same lightning physics otherwise).
"""

import numpy as np

GEV2G = 1.78e-24  # same conversion as macrolimits

# Defaults: Starkman+2021 Earth values, with rho_DM and CL matched to Picker+2026.
RHO_DM = 0.3 * GEV2G  # g/cm^3
R_EARTH = 6.371e8  # cm
F_TS = 3e-3  # fraction of Earth's surface under a thunderstorm
T_EXPOSURE = 2 * 3.156e7  # s; the paper's "two years" of non-observation
N_OBS = 2.3  # 90% CL Poisson, as in Picker+2026
SIGMA_FLOOR = 6e-8  # cm^2; t_I^0 >= 1 microsecond (paper Eq. tI0min)
K_SIGMA = 1e8  # cm^4/s^2; lambda_e^macro >= lambda_e^natural (code's sigma_factor)
SIG_OVER_M_MAX = 1e-4  # cm^2/g

# Standard halo model in the Galactic frame (as in the original code), km/s.
V_VIR, V_ESC, V_CIRC = 250.0, 550.0, 220.0
V_INFALL = np.hypot(42.1, 11.2)  # Sun escape speed at 1 AU (+) Earth escape speed

# Per-planet inputs. Jupiter's infall speed combines the Sun's escape speed at 5.2 AU (18.5 km/s)
# with Jupiter's (59.5 km/s); the original code used the Sun's escape speed at Earth instead.
PLANETS = {
    "earth": {"r": R_EARTH, "f_ts": F_TS, "k_sigma": K_SIGMA, "v_infall": V_INFALL},
    "jupiter": {"r": 6.99e9, "f_ts": 0.1, "k_sigma": 5 * K_SIGMA, "v_infall": float(np.hypot(18.5, 59.5))},
}


def speeds(dv=10.0, v_infall=V_INFALL):
    """Planet-frame speeds v_rel [cm/s] and their SHM weights, summing to 1."""
    g = np.arange(-V_ESC, V_ESC + dv, dv)
    vx, vy, vz = np.meshgrid(g, g, g, indexing="ij", sparse=True)
    vgal = np.sqrt(vx**2 + vy**2 + vz**2)
    w = np.exp(-((vgal / V_VIR) ** 2)) * (vgal <= V_ESC)
    vrel = np.sqrt(v_infall**2 + vx**2 + (vy - V_CIRC) ** 2 + vz**2)
    vrel, w = np.broadcast_arrays(vrel, w)
    keep = w > 0
    return vrel[keep] * 1e5, w[keep] / w[keep].sum()


def mean_speed_above(sarray, dv=10.0, planet="earth"):
    """<v_rel ; sigma>: integral of v_rel f(v) over speeds that seed lightning."""
    p = PLANETS[planet]
    v, w = speeds(dv, p["v_infall"])
    order = np.argsort(v)
    v, vw = v[order], (v * w)[order]
    tail = np.cumsum(vw[::-1])[::-1]  # integral of v f over speeds >= v[i]
    # sigma >= K/v^2  <=>  v >= sqrt(K/sigma)
    vmin = np.sqrt(p["k_sigma"] / np.asarray(sarray, float))
    idx = np.searchsorted(v, vmin)
    out = np.where(idx < len(v), tail[np.minimum(idx, len(v) - 1)], 0.0)
    return np.where(np.asarray(sarray) >= SIGMA_FLOOR, out, 0.0)


def lightning(marray, sarray, T=T_EXPOSURE, dv=10.0, planet="earth"):
    """f_DM grid, shape (len(sarray), len(marray)), masses in GeV, sigma in cm^2.

    Unconstrained cells are set to 99, matching macrolimits' convention.
    """
    m = np.asarray(marray) * GEV2G  # g
    s = np.asarray(sarray)[:, None]
    p = PLANETS[planet]
    exposure = RHO_DM * np.pi * p["r"] ** 2 * p["f_ts"] * T  # g s / cm
    n_exp = exposure * mean_speed_above(s, dv, planet) / m[None, :]  # events at f_DM = 1
    with np.errstate(divide="ignore"):
        f = np.where(n_exp > 0, N_OBS / n_exp, 99.0)
    f[(s / m[None, :]) > SIG_OVER_M_MAX] = 99.0
    return np.minimum(f, 99.0)


if __name__ == "__main__":
    # Check against Starkman+2021: ~1e6 g reach at f_DM = 1 for large sigma.
    m = np.logspace(20, 35, 3001)
    f = lightning(m, [1.0])[0]
    m_reach_g = m[f <= 1].max() * GEV2G
    print(f"mass reach at f_DM=1, sigma=1 cm^2: {m_reach_g:.2e} g")
    assert 3e5 < m_reach_g < 3e6, m_reach_g
    # Below the floor nothing is constrained; f_DM scales as m (flux-limited).
    assert (lightning(m, [1e-9]) == 99).all()
    f1, f2 = lightning([1e26, 1e27], [1e-3])[0]
    assert np.isclose(f2 / f1, 10)
    mj = m[lightning(m, [1.0], planet="jupiter")[0] <= 1].max() * GEV2G
    print(f"Jupiter mass reach at f_DM=1, sigma=1 cm^2: {mj:.2e} g")
    assert mj > 1e3 * m_reach_g  # ~120x area and ~33x storm fraction
    print("ok")
