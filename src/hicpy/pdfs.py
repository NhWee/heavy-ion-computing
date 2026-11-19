"""Parton distribution functions for Module 5, from PYTHIA's own PDF grids.

WHY NOT LHAPDF
--------------
LHAPDF is the standard library, but it is not on PyPI and its PDF-set servers
were unreachable from the cloud session.  PYTHIA ships LHAPDF-format grids for
the PDFs its tunes use (NNPDF2.3/3.1, CTEQ6) and the EPS09 LO nuclear
modification for lead -- enough for every question in this module.  The same
code runs on the conda `pythia8` in the WSL2 environment.

CONVENTIONS
-----------
xf(id, x, Q2) returns x times the number density of parton id at momentum
fraction x and scale Q^2 (GeV^2).  PDG codes: 1 d, 2 u, 3 s, 4 c, 5 b, 21 g,
negative = antiquark.  Valence = q - qbar.

NUCLEAR PDFs
------------
EPS09 gives, per nucleon, f_i^A(x, Q2) = R_i^A(x, Q2) f_i^nucleon(x, Q2) with
isospin: f^nucleon = (Z f^p + (A - Z) f^n) / A, f^n from f^p by u <-> d.  The
ratio R_i^A < 1 at small x is SHADOWING, > 1 around x ~ 0.1 ANTISHADOWING,
< 1 at x ~ 0.3-0.7 the EMC effect, and > 1 near x -> 1 Fermi motion.  EPS09
uncertainties are Hessian (15 eigenvector pairs, members 1-30, member 0 is
the central fit) at 90 % confidence level.
"""
from __future__ import annotations

import os
import pathlib

import numpy as np

PB208 = 1000822080
PROTON_SET = "NNPDF23_lo_as_0130_qed_0000.dat"     # the PDF the Monash tune uses


def pythia_module():
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    return pythia8


def pdfdata_path() -> str:
    """Directory of PYTHIA's PDF grids, found the way PYTHIA finds xmldoc."""
    if os.environ.get("PYTHIA8DATA"):
        cand = pathlib.Path(os.environ["PYTHIA8DATA"]).parent / "pdfdata"
        if cand.exists():
            return str(cand)
    py = pythia_module().Pythia("", False)
    xml = pathlib.Path(py.settings.word("xmlPath"))
    cand = xml.parent / "pdfdata"
    if cand.exists():
        return str(cand)
    raise FileNotFoundError(f"PYTHIA pdfdata not found next to {xml}")


class ProtonPDF:
    def __init__(self, grid: str = PROTON_SET, beam: int = 2212):
        self.p8 = pythia_module()
        self.grid = grid
        self.pdf = self.p8.LHAGrid1(beam, grid, pdfdata_path())

    def xf(self, pid, x, q2):
        x = np.atleast_1d(np.asarray(x, float))
        return np.array([self.pdf.xf(int(pid), float(xi), float(q2)) for xi in x])

    def valence(self, flavour, x, q2):
        return self.xf(flavour, x, q2) - self.xf(-flavour, x, q2)


class NuclearRatio:
    """EPS09 LO for Pb-208: R_i(x, Q2) = f_i^Pb / f_i^(isospin-averaged nucleon)."""
    N_MEMBERS = 31

    def __init__(self, proton: ProtonPDF | None = None):
        self.p8 = pythia_module()
        self.proton = proton or ProtonPDF()
        self.neutron = ProtonPDF(self.proton.grid, beam=2112)
        self._members = {}
        self.A, self.Z = 208, 82

    def _member(self, m):
        if m not in self._members:
            self._members[m] = self.p8.EPS09(PB208, 1, m, pdfdata_path(), self.proton.pdf)
        return self._members[m]

    def nucleon_xf(self, pid, x, q2):
        return (self.Z * self.proton.xf(pid, x, q2) + (self.A - self.Z) * self.neutron.xf(pid, x, q2)) / self.A

    def ratio(self, pid, x, q2, member=0):
        x = np.atleast_1d(np.asarray(x, float))
        e = self._member(member)
        num = np.array([e.xf(int(pid), float(xi), float(q2)) for xi in x])
        return num / self.nucleon_xf(pid, x, q2)

    def ratio_band(self, pid, x, q2):
        """Central value and symmetric Hessian uncertainty (90 % CL)."""
        r0 = self.ratio(pid, x, q2, 0)
        d2 = np.zeros_like(r0)
        for k in range(1, 16):
            d2 += (self.ratio(pid, x, q2, 2 * k - 1) - self.ratio(pid, x, q2, 2 * k)) ** 2
        return r0, 0.5 * np.sqrt(d2)


def integrate_logx(f, x_min=1e-6, x_max=1.0, n=4000):
    """Integral of f(x) dx on a log grid (trapezoid in x)."""
    x = np.geomspace(x_min, x_max * (1 - 1e-9), n)
    return float(np.trapezoid(f(x), x))


__all__ = ["ProtonPDF", "NuclearRatio", "pdfdata_path", "integrate_logx", "PROTON_SET", "PB208"]
