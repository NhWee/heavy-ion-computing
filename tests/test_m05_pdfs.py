"""hicpy.pdfs against facts that do not depend on any fit: sum rules, isospin,
and the sign pattern of nuclear modifications."""
import numpy as np
import pytest

try:
    from hicpy.pdfs import pythia_module
    pythia_module()
except ImportError:
    pytest.skip("no PYTHIA Python module with PDF classes (pip install pythia8mc)", allow_module_level=True)

from hicpy.pdfs import ProtonPDF, NuclearRatio, integrate_logx   # noqa: E402


@pytest.fixture(scope="module")
def proton():
    return ProtonPDF()


@pytest.mark.parametrize("q2", [10.0, 1e4])
def test_valence_and_momentum_sum_rules(proton, q2):
    assert integrate_logx(lambda x: proton.valence(2, x, q2) / x) == pytest.approx(2.0, abs=0.01)
    assert integrate_logx(lambda x: proton.valence(1, x, q2) / x) == pytest.approx(1.0, abs=0.01)
    assert integrate_logx(lambda x: proton.valence(3, x, q2) / x) == pytest.approx(0.0, abs=1e-3)
    mom = integrate_logx(lambda x: sum(proton.xf(i, x, q2) for i in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21)))
    assert mom == pytest.approx(1.0, abs=0.01)


def test_gluon_grows_at_small_x_with_q2(proton):
    assert proton.xf(21, 1e-4, 1e4)[0] > 3 * proton.xf(21, 1e-4, 10.0)[0]


def test_isospin_neutron_u_is_proton_d(proton):
    neutron = ProtonPDF(beam=2112)
    assert neutron.xf(2, 0.1, 10.0)[0] == pytest.approx(proton.xf(1, 0.1, 10.0)[0], rel=1e-6)


def test_eps09_shadowing_and_antishadowing(proton):
    nr = NuclearRatio(proton)
    r, d = nr.ratio_band(21, [1e-3, 0.1], 10.0)
    assert r[0] < 0.9 and r[1] > 1.05          # shadowing at small x, antishadowing at 0.1
    assert np.all(d > 0)
