"""hicpy.eventrecord on hand-made event records (no PYTHIA needed).

The synthetic record mimics the two features of a real PYTHIA history that
broke the first converter: a mother stored AFTER its daughter (backward ISR
evolution) and a vertex with two incoming particles (a string).
"""
import numpy as np
import pytest

pyhepmc = pytest.importorskip("pyhepmc")

from hicpy.eventrecord import (Rec, to_hepmc, topological_order, vertex_balance,  # noqa: E402
                               global_checks, ancestry)


def boosted_pair(m, pz, frac=0.5):
    """Two massless daughters sharing a parent (E, pz) along z."""
    e = np.hypot(m, pz)
    a = (0.0, 0.0, frac * pz, frac * e)
    b = (0.0, 0.0, (1 - frac) * pz, (1 - frac) * e)
    return a, b


def toy_record():
    """beams 1, 2 -> (via a gluon 5 stored AFTER its daughter 3) -> string of
    gluons 3 and 4 -> neutral hadron 6 -> pi+ pi- (final).  Each beam loses a
    gluon only, so the two remnants keep the proton charge: total +2."""
    E = 10.0
    b1 = (0, 0, E, E); b2 = (0, 0, -E, E)
    q = (0, 0, 2.0, 2.0)           # parton entering the 'hard process'
    g = (0, 0, -1.0, 1.0)
    h = (0, 0, 1.0, 3.0)           # hadron from the string (q + g)
    d1, d2 = boosted_pair(np.sqrt(3.0 ** 2 - 1.0 ** 2), 1.0, 0.3)
    rem1 = tuple(np.subtract(b1, (0, 0, 2.0, 2.0)))
    rem2 = tuple(np.subtract(b2, (0, 0, -1.0, 1.0)))
    return [
        Rec(1, 2212, -12, 4, (), b1, 3),
        Rec(2, 2212, -12, 4, (), b2, 3),
        Rec(3, 21, -21, 21, (5,), q, 0),         # mother 5 has a HIGHER index
        Rec(4, 21, -21, 21, (2,), g, 0),
        Rec(5, 21, -41, 41, (1,), q, 0),
        Rec(6, 113, -83, 2, (3, 4), h, 0),       # two incoming: a string
        Rec(7, 211, 91, 1, (6,), d1, 3),
        Rec(8, -211, 91, 1, (6,), d2, -3),
        Rec(9, 2212, 63, 1, (1,), rem1, 3),
        Rec(10, 2212, 63, 1, (2,), rem2, 3),
    ]


def test_topological_order_puts_mothers_first():
    order = [r.index for r in topological_order(toy_record())]
    pos = {i: k for k, i in enumerate(order)}
    for r in toy_record():
        for m in r.mothers:
            assert pos[m] < pos[r.index], (m, r.index)


def test_cycle_is_rejected():
    recs = [Rec(1, 1, 0, 0, (2,), (0, 0, 0, 1), 0), Rec(2, 1, 0, 0, (1,), (0, 0, 0, 1), 0)]
    with pytest.raises(ValueError):
        topological_order(recs)


def test_hepmc_round_trip(tmp_path):
    ev = to_hepmc(toy_record(), event_number=7)
    path = tmp_path / "toy.hepmc"
    with pyhepmc.open(str(path), "w") as f:
        f.write(ev)
    with pyhepmc.open(str(path)) as f:
        back = f.read()                      # keep the event alive (troubleshooting #18)
    assert back is not None and back.event_number == 7
    assert len(back.particles) == len(ev.particles) == 10
    assert len(back.vertices) == len(ev.vertices)
    string = [v for v in back.vertices if len(v.particles_in) == 2]
    assert len(string) == 1 and {p.pid for p in string[0].particles_out} == {113}


def test_global_checks_and_decay_vertex_balance():
    ev = to_hepmc(toy_record())
    g = global_checks(ev)
    assert g["n_final"] == 4                       # pi+, pi-, two remnants
    assert g["charge_final"] == pytest.approx(2.0)
    assert g["p4_rel_diff"] < 1e-12
    vb = vertex_balance(ev)
    assert vb["9"]["n"] == 1 and vb["9"]["conserved"] == 1      # the decay
    assert vb["8"]["conserved"] == 1                              # the string


def test_ancestry_reaches_a_beam():
    ev = to_hepmc(toy_record())
    pion = next(p for p in ev.particles if p.pid == 211)
    chain = ancestry(pion)
    assert chain[0] == (211, 91)
    assert chain[-1][0] == 2212
