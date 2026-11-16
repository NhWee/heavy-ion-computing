"""scripts/check_results_drift.py: which values are compared how strictly.

Central values must reproduce to 1e-6 across machines; HESSE-derived
uncertainties only to 1e-2 (they move by ~1e-3 between machines); the
ill-conditioned Barlow shortcut is not compared at all."""
import sys

from conftest import REPO

sys.path.insert(0, str(REPO / "scripts"))
import check_results_drift as cd  # noqa: E402


def test_classification_of_paths():
    assert not cd.is_uncertainty("/fits[0]/params[1]")
    assert not cd.is_uncertainty("/inclusive/summary/T/syst")
    assert not cd.is_uncertainty("/inclusive/variations[4]/T_delta")
    for p in ("/fits[0]/errors[1]", "/fits[0]/covariance[0][1]", "/checks[3]/pull",
              "/inclusive/summary/T/stat", "/z_fit/m_err", "/tsallis/pi+-/poisson/T/pull_sd"):
        assert cd.is_uncertainty(p), p
    assert cd.is_informational("/inclusive/variations[4]/T_barlow_sigma")
    assert cd.is_informational("/inclusive/variations[4]/n_significant_by_barlow")
    assert not cd.is_informational("/inclusive/variations[4]/T_significant")


def test_walk_applies_the_right_tolerance():
    problems = []
    cd.walk({"params": [1.0], "errors": [0.1]}, {"params": [1.0 + 5e-6], "errors": [0.1005]},
            "", problems)
    assert len(problems) == 1 and "params" in problems[0] and "VALUE" in problems[0]
    problems = []
    cd.walk({"errors": [0.1], "x_barlow_sigma": 1.0}, {"errors": [0.12], "x_barlow_sigma": 5.0},
            "", problems)
    assert len(problems) == 1 and "errors" in problems[0]        # barlow ignored
    problems = []
    cd.walk({"n": 3, "flag": True}, {"n": 4, "flag": True}, "", problems)
    assert problems == ["/n: 3 vs 4"]
