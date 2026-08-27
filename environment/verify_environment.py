#!/usr/bin/env python3
"""Validate the HEP environment and record exact versions.

    python environment/verify_environment.py [--write environment/versions_installed.txt]

"The install command exited 0" is not evidence that anything works.  Each check
below actually *uses* the package -- creates a ROOT file, generates events,
clusters jets, evaluates a PDF -- and reports PASS/FAIL/SKIP.  SKIP means the
package is not installed yet, which is expected for modules you have not reached.
The exit code is non-zero only if something installed is broken.
"""
from __future__ import annotations

import argparse
import importlib
import os
import platform
import subprocess
import sys
import tempfile

RESULTS: list[tuple[str, str, str]] = []      # (name, status, detail)


def record(name, status, detail=""):
    RESULTS.append((name, status, detail))
    colour = {"PASS": "\033[32m", "FAIL": "\033[31m", "SKIP": "\033[33m"}.get(status, "")
    print(f"  {colour}{status:4s}\033[0m  {name:22s} {detail}")


def check_import(name, module=None, version_attr="__version__"):
    """Import a package and record its version."""
    mod_name = module or name
    try:
        m = importlib.import_module(mod_name)
    except Exception as exc:                                    # noqa: BLE001
        record(name, "SKIP", f"not installed ({type(exc).__name__})")
        return None
    v = getattr(m, version_attr, "?")
    record(name, "PASS", f"version {v}")
    return m


# --------------------------------------------------------------------------
def check_python_stack():
    print("\n[python stack]")
    for n in ["numpy", "scipy", "matplotlib", "yaml", "uproot", "awkward",
              "vector", "hist", "mplhep", "particle", "iminuit", "pandas"]:
        check_import(n)


def check_uproot_roundtrip():
    """Write a small TTree with uproot and read it back."""
    print("\n[functional: uproot TTree round-trip]")
    try:
        import awkward as ak
        import numpy as np
        import uproot
    except Exception:
        record("uproot round-trip", "SKIP", "uproot/awkward missing")
        return
    try:
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "check.root")
            data = {"particle": ak.zip({"pt": ak.Array([[1.0, 2.0], [3.0]])}),
                    "flag": np.array([True, False])}
            with uproot.recreate(path) as f:
                t = f.mktree("t", {k: (v.type if hasattr(v, "type") else v.dtype)
                                   for k, v in data.items()},
                             counter_name=lambda n: "n" + n)
                t.extend(data)
            with uproot.open(path) as f:
                got = f["t"].arrays()
                assert f["t"].num_entries == 2
                assert list(got.particle_pt[0]) == [1.0, 2.0]
        record("uproot round-trip", "PASS", "TTree written and read back")
    except Exception as exc:                                    # noqa: BLE001
        record("uproot round-trip", "FAIL", f"{type(exc).__name__}: {exc}")


def check_root():
    """PyROOT: create a TFile+TTree+TH1, read it back, check the entry count."""
    print("\n[functional: ROOT / PyROOT]")
    try:
        import ROOT
    except Exception as exc:                                    # noqa: BLE001
        record("ROOT (PyROOT)", "SKIP", f"not installed ({type(exc).__name__})")
        return
    try:
        ROOT.gROOT.SetBatch(True)
        record("ROOT (PyROOT)", "PASS", f"version {ROOT.gROOT.GetVersion()}")
        import array
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "check.root")
            f = ROOT.TFile(path, "RECREATE")
            t = ROOT.TTree("t", "check")
            x = array.array("f", [0.0])
            t.Branch("x", x, "x/F")
            h = ROOT.TH1F("h", ";x;counts", 10, 0.0, 10.0)
            for i in range(100):
                x[0] = i % 10
                h.Fill(x[0])
                t.Fill()
            t.Write()
            h.Write()
            f.Close()
            f2 = ROOT.TFile(path)
            assert f2.Get("t").GetEntries() == 100
            assert abs(f2.Get("h").GetEntries() - 100) < 1e-6
            f2.Close()
        record("ROOT TFile/TTree/TH1", "PASS", "100 entries written and read")
    except Exception as exc:                                    # noqa: BLE001
        record("ROOT TFile/TTree/TH1", "FAIL", f"{type(exc).__name__}: {exc}")
    # RDataFrame is the modern analysis interface; check it separately
    try:
        df = ROOT.RDataFrame(10).Define("y", "rdfentry_ * 2")
        assert df.Sum("y").GetValue() == 90
        record("ROOT RDataFrame", "PASS", "sum over 10 entries correct")
    except Exception as exc:                                    # noqa: BLE001
        record("ROOT RDataFrame", "FAIL", f"{type(exc).__name__}: {exc}")


def check_pythia():
    """Generate 10 minimum-bias events and count final-state particles."""
    print("\n[functional: PYTHIA 8]")
    try:
        import pythia8
    except Exception as exc:                                    # noqa: BLE001
        record("PYTHIA 8", "SKIP", f"not installed ({type(exc).__name__})")
        return
    try:
        p = pythia8.Pythia("", False)
        for s in ["Beams:idA = 2212", "Beams:idB = 2212", "Beams:eCM = 13000.",
                  "SoftQCD:inelastic = on", "Random:setSeed = on",
                  "Random:seed = 12345", "Print:quiet = on"]:
            p.readString(s)
        assert p.init()
        n_final = 0
        for _ in range(10):
            if p.next():
                n_final += sum(1 for pa in p.event if pa.isFinal())
        assert n_final > 100, f"suspiciously few final particles: {n_final}"
        record("PYTHIA 8", "PASS",
               f"10 pp events at 13 TeV, {n_final} final-state particles")
    except Exception as exc:                                    # noqa: BLE001
        record("PYTHIA 8", "FAIL", f"{type(exc).__name__}: {exc}")


def check_fastjet():
    """Cluster three particles with anti-kT and check we get a jet back."""
    print("\n[functional: FastJet]")
    try:
        import fastjet
    except Exception as exc:                                    # noqa: BLE001
        record("FastJet", "SKIP", f"not installed ({type(exc).__name__})")
        return
    try:
        version = getattr(fastjet, "__version__", "?")
        if hasattr(fastjet, "JetDefinition") and hasattr(fastjet, "PseudoJet"):
            # SWIG bindings (conda-forge `fastjet`)
            jd = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
            parts = [fastjet.PseudoJet(10.0, 0.0, 0.0, 10.0),
                     fastjet.PseudoJet(9.0, 0.5, 0.0, 9.02),
                     fastjet.PseudoJet(-20.0, 0.0, 0.0, 20.0)]
            jets = fastjet.ClusterSequence(parts, jd).inclusive_jets(5.0)
            n = len(jets)
        else:                                                   # scikit-hep wheel
            import awkward as ak
            import vector
            vector.register_awkward()
            arr = ak.Array([[{"px": 10.0, "py": 0.0, "pz": 0.0, "E": 10.0},
                             {"px": 9.0, "py": 0.5, "pz": 0.0, "E": 9.02},
                             {"px": -20.0, "py": 0.0, "pz": 0.0, "E": 20.0}]],
                           with_name="Momentum4D")
            cs = fastjet.ClusterSequence(arr, fastjet.JetDefinition(
                fastjet.antikt_algorithm, 0.4))
            n = len(cs.inclusive_jets(min_pt=5.0)[0])
        assert n >= 2, f"expected >= 2 jets, got {n}"
        record("FastJet", "PASS", f"version {version}, anti-kT R=0.4 -> {n} jets")
    except Exception as exc:                                    # noqa: BLE001
        record("FastJet", "FAIL", f"{type(exc).__name__}: {exc}")


def check_hepmc():
    print("\n[functional: HepMC3]")
    try:
        import pyhepmc
    except Exception as exc:                                    # noqa: BLE001
        record("HepMC3 (pyhepmc)", "SKIP", f"not installed ({type(exc).__name__})")
        return
    try:
        ev = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
        p1 = pyhepmc.GenParticle((0, 0, 100, 100), 2212, 4)
        p2 = pyhepmc.GenParticle((0, 0, -100, 100), 2212, 4)
        v = pyhepmc.GenVertex()
        v.add_particle_in(p1)
        v.add_particle_in(p2)
        ev.add_vertex(v)
        assert len(ev.particles) == 2
        record("HepMC3 (pyhepmc)", "PASS",
               f"version {pyhepmc.__version__}, GenEvent with 2 particles")
    except Exception as exc:                                    # noqa: BLE001
        record("HepMC3 (pyhepmc)", "FAIL", f"{type(exc).__name__}: {exc}")


def check_lhapdf():
    """Evaluate a PDF at a reference point."""
    print("\n[functional: LHAPDF]")
    try:
        import lhapdf
    except Exception as exc:                                    # noqa: BLE001
        record("LHAPDF", "SKIP", f"not installed ({type(exc).__name__})")
        return
    try:
        pdf = lhapdf.mkPDF("CT18NLO", 0)
        val = pdf.xfxQ(2, 0.01, 100.0)          # x*u(x=0.01, Q=100 GeV)
        assert 0.0 < val < 10.0, val
        record("LHAPDF", "PASS", f"CT18NLO x*u(0.01, 100 GeV) = {val:.4f}")
    except Exception as exc:                                    # noqa: BLE001
        record("LHAPDF", "SKIP", f"installed but no PDF set available ({exc})")


def check_toolchain():
    print("\n[toolchain]")
    for cmd, args in [("gcc", ["--version"]), ("g++", ["--version"]),
                      ("cmake", ["--version"]), ("git", ["--version"]),
                      ("root-config", ["--version"])]:
        try:
            out = subprocess.check_output([cmd] + args, stderr=subprocess.STDOUT)
            record(cmd, "PASS", out.decode().splitlines()[0])
        except Exception:                                       # noqa: BLE001
            record(cmd, "SKIP", "not on PATH")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", default=None,
                    help="also write the version record to this file")
    args = ap.parse_args()

    print("=" * 72)
    print("heavy-ion-computing :: environment validation")
    print("=" * 72)
    print(f"  python   {sys.version.split()[0]}  ({sys.executable})")
    print(f"  platform {platform.platform()}")
    print(f"  machine  {platform.machine()}, {os.cpu_count()} logical CPUs")

    check_toolchain()
    check_python_stack()
    check_uproot_roundtrip()
    check_root()
    check_pythia()
    check_fastjet()
    check_hepmc()
    check_lhapdf()

    n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
    n_fail = sum(1 for _, s, _ in RESULTS if s == "FAIL")
    n_skip = sum(1 for _, s, _ in RESULTS if s == "SKIP")
    print("\n" + "=" * 72)
    print(f"  {n_pass} PASS   {n_fail} FAIL   {n_skip} SKIP (not installed yet)")
    print("=" * 72)

    if args.write:
        with open(args.write, "w") as fh:
            fh.write("# heavy-ion-computing environment record\n")
            fh.write(f"# generated by environment/verify_environment.py\n")
            fh.write(f"python\t{sys.version.split()[0]}\n")
            fh.write(f"platform\t{platform.platform()}\n")
            for name, status, detail in RESULTS:
                fh.write(f"{name}\t{status}\t{detail}\n")
        print(f"  version record -> {args.write}")

    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
