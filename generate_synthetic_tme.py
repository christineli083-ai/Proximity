"""
Synthetic melanoma TME spatial data generator.

WHY SYNTHETIC: Real HTAN/TCGA-SKCM multiplexed IF data requires registered
access + multi-GB downloads + QuPath cell segmentation, none of which fit a
4.5-hour build window. This generator instead encodes the SPECIFIC mechanism
described in the accompanying literature review: cancer-associated fibroblasts
(CAFs) remodel the ECM and create physical barriers that push CD8+ T-cells
away from the tumor margin, and this exclusion effect intensifies as tumors
become thicker/more aggressive (Falcone et al. 2020; Sikorski et al. 2025).

We simulate N virtual "patients," each with:
  - a melanoma "tumor mass" as a 2D point cloud (malignant melanocytes)
  - a Breslow depth (mm), drawn from a realistic clinical distribution
  - CD8+ T-cells scattered around the tumor, with their MEAN distance from
    the tumor margin increasing with Breslow depth (the exclusion effect)
  - a survival outcome probabilistically linked to that same exclusion
    distance, independent noise added so the correlation isn't trivial/perfect

This is a controlled proof-of-concept dataset: the ground-truth relationship
is designed in, so the analysis pipeline can be validated end-to-end. It is
NOT real patient data and every output built from it should say so.
"""

import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(42)  # fixed seed -> reproducible, citable in appendix

OUT_DIR = Path(__file__).parent
N_PATIENTS = 60          # comparable to real published melanoma spatial cohorts (e.g. n=54, n=60 range)
N_TUMOR_CELLS = 250
N_TCELLS = 120
TUMOR_RADIUS_UM = 300    # rough tumor nest radius in micrometers


def sample_breslow_depth():
    # Clinically-informed: right-skewed, most cases thin/intermediate, some thick.
    return float(np.clip(RNG.lognormal(mean=0.3, sigma=0.6), 0.2, 12.0))


def simulate_patient(patient_id: int) -> dict:
    breslow = sample_breslow_depth()

    # Tumor mass: melanocytes clustered near origin (roughly circular nest,
    # with some irregularity to mimic real tumor margins)
    angles = RNG.uniform(0, 2 * np.pi, N_TUMOR_CELLS)
    radii = TUMOR_RADIUS_UM * np.sqrt(RNG.uniform(0, 1, N_TUMOR_CELLS)) * (
        1 + 0.15 * RNG.standard_normal(N_TUMOR_CELLS)
    )
    tumor_x = radii * np.cos(angles)
    tumor_y = radii * np.sin(angles)

    # CD8+ T-cells: base ring around the tumor margin, but pushed further out
    # as Breslow depth increases (stromal exclusion effect). Effect size and
    # noise are explicit constants below so they can be cited/justified.
    EXCLUSION_SLOPE_UM_PER_MM = 22.0   # um of extra exclusion per mm Breslow depth (population average)
    BASE_INFILTRATION_UM = 15.0        # um, thin tumors: T-cells reach near margin
    NOISE_SD_UM = 95.0                 # biological + measurement noise (per-cell)

    # Patient-level heterogeneity: not every patient follows the population
    # trend equally (real tumors vary in stromal density, immune fitness,
    # etc. beyond what Breslow depth alone captures). This is what keeps the
    # correlation realistic (moderate) instead of a trivially perfect line.
    patient_slope = max(0.0, RNG.normal(EXCLUSION_SLOPE_UM_PER_MM, 14.0))
    patient_baseline_shift = RNG.normal(0, 25.0)

    mean_offset = BASE_INFILTRATION_UM + patient_baseline_shift + patient_slope * breslow
    tcell_angles = RNG.uniform(0, 2 * np.pi, N_TCELLS)
    tcell_radial_offset = np.clip(
        RNG.normal(mean_offset, NOISE_SD_UM, N_TCELLS), -TUMOR_RADIUS_UM * 0.9, None
    )
    tcell_radii = TUMOR_RADIUS_UM + tcell_radial_offset
    tcell_x = tcell_radii * np.cos(tcell_angles)
    tcell_y = tcell_radii * np.sin(tcell_angles)

    return {
        "patient_id": patient_id,
        "breslow_depth_mm": breslow,
        "tumor_cells": np.column_stack([tumor_x, tumor_y]),
        "tcells": np.column_stack([tcell_x, tcell_y]),
    }


def simulate_survival(nn_mean_um: np.ndarray, breslow: np.ndarray) -> np.ndarray:
    """
    5-year survival probability as a logistic function of NN distance and
    Breslow depth (both independently reported as prognostic in the lit
    review's cited sources), plus noise -> binary outcome per patient.
    Coefficients are illustrative, not fitted to real outcome data — this
    must be stated as a modeling assumption, not a validated clinical result.
    """
    z = 2.2 - 0.012 * nn_mean_um - 0.15 * breslow + RNG.normal(0, 0.6, len(nn_mean_um))
    p_survival = 1 / (1 + np.exp(-z))
    return (RNG.uniform(0, 1, len(nn_mean_um)) < p_survival).astype(int)


def main():
    patients = [simulate_patient(i) for i in range(N_PATIENTS)]

    # Save per-patient raw coordinates (appendix-ready, one CSV per patient
    # would be excessive -> pack into a single long-format file)
    rows = []
    for p in patients:
        for x, y in p["tumor_cells"]:
            rows.append({"patient_id": p["patient_id"], "cell_type": "melanocyte", "x_um": x, "y_um": y})
        for x, y in p["tcells"]:
            rows.append({"patient_id": p["patient_id"], "cell_type": "CD8_Tcell", "x_um": x, "y_um": y})
    coords_df = pd.DataFrame(rows)
    coords_df.to_csv(OUT_DIR / "synthetic_cell_coordinates.csv", index=False)

    clinical_df = pd.DataFrame(
        {"patient_id": [p["patient_id"] for p in patients],
         "breslow_depth_mm": [p["breslow_depth_mm"] for p in patients]}
    )
    clinical_df.to_csv(OUT_DIR / "synthetic_clinical.csv", index=False)

    print(f"Generated {N_PATIENTS} synthetic patients")
    print(f"  -> {OUT_DIR / 'synthetic_cell_coordinates.csv'}")
    print(f"  -> {OUT_DIR / 'synthetic_clinical.csv'}")


if __name__ == "__main__":
    main()
