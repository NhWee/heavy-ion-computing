"""Reference particle properties (world averages) with their source.

Two kinds of "PDG value" appear in this project and must not be mixed up:

* generator INPUTS -- what the toy generator samples from, i.e. the closure
  truth.  They live in hicpy.toy_generator and are frozen there on purpose:
  changing them changes every committed event sample.
* external REFERENCES -- the current world average an analysis compares a
  measurement with.  They live here.

For the Z they differ.  The toy generator uses 91.1876 / 2.4952 GeV, the
long-standing LEP-era values (labelled "PDG 2024" until 2026-11-07, which was
wrong).  The values below were read on 2026-11-07 from the PDG mass-width
tables shipped with the `particle` package 1.1.0:

    mass_width_2025.txt  (RPP 2025):  m_Z = 91.1880 +- 0.0020,  Gamma_Z = 2.4955 +- 0.0023
    mass_width_2026.txt  (RPP 2026):  m_Z = 91.1879 +- 0.0020,  Gamma_Z = 2.4955 +- 0.0023

The 0.3 MeV difference from the toy input is far below anything Module 1 can
resolve (its Z-mass error is ~100 MeV); it matters only for labelling.
"""
from __future__ import annotations

REFERENCE_EDITION = "RPP 2026 (particle 1.1.0, data/mass_width_2026.txt)"

Z_MASS = 91.1879          # GeV, +- 0.0020
Z_MASS_ERR = 0.0020
Z_WIDTH = 2.4955          # GeV, +- 0.0023
Z_WIDTH_ERR = 0.0023

__all__ = ["REFERENCE_EDITION", "Z_MASS", "Z_MASS_ERR", "Z_WIDTH", "Z_WIDTH_ERR"]
