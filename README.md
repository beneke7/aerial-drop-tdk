# aerial-drop-tdk

Companion code to the TDK thesis *A légi tűzoltás fizikai modellezése és hatékonyságának vizsgálata* (ELTE TTK, 2026).

## Contents

| Path | What |
|---|---|
| `cases/rouaix_coarse`, `rouaix_fineA`, `rouaix_fineB` | OpenFOAM v2512 case setup of the Rouaix et al. (2023) Case 1 jet in crossflow (Section 4.1) |
| `cases/bambi_full`, `bambi_f60`, `bambi_rotor_f60` | Bambi bucket release, full and 60 % fill, with and without rotor (Section 4.2) |
| `cases/velocity_matched_jets`, `velocity_matched_curtain` | Velocity-matched releases (Chapter 5) |
| `scripts/handover_crossing.py`, `handover_drops.py`, `handover_flux_check.py` | VOF-to-parcel conversion at the surface S and the parcel flights to the ground (drag, breakup, evaporation, wind, turbulence) |
| `scripts/evaporation.py` | Single-drop evaporation (Denner 2026) |
| `scripts/band_interp.py`, `metric_measured_cases.py` | Reconstruction of measured cup patterns from published contours and their evaluation |
| `src/aerial_drop/hold.py` | Evaluation of ground patterns: discrete and continuous efficiency |

Each case has the standard `0.orig/`, `constant/` and `system/` folders; meshes, fields and logs are not included. The settings are summarised in Appendix C of the thesis.

## Use

```
pip install -e .
python scripts/handover_crossing.py <liquidBand.vtm.series> --U 21.607 ...
```

The scripts read OpenFOAM exports and digitised measurements that are not part of this repository (size).

## Not published

HPC job and run scripts, plotting and auxiliary analysis scripts are not published.

## Licence

MIT
