"""The PYTHIA card must state every physics choice explicitly (generators/README.md).
This test needs no PYTHIA: it only reads the card."""
import pathlib
import re
import sys

import pytest

from conftest import REPO

sys.path.insert(0, str(REPO / "scripts"))
from m02_run_pythia import card_value   # noqa: E402

CARD = REPO / "generators" / "pythia_mb_13TeV.cmnd"
REQUIRED = {
    "Beams:idA": "2212", "Beams:idB": "2212", "Beams:eCM": "13000.",
    "SoftQCD:inelastic": "on", "Tune:pp": "14", "PDF:pSet": "13",
    "ParticleDecays:limitTau0": "on", "ParticleDecays:tau0Max": "10.",
    "Random:setSeed": "on", "Random:seed": "20260827",
    "Main:numberOfEvents": "20000",
}


@pytest.fixture(scope="module")
def card_text():
    return CARD.read_text()


@pytest.mark.parametrize("key,value", sorted(REQUIRED.items()))
def test_card_states_every_physics_choice(card_text, key, value):
    assert card_value(card_text, key) == value


def test_card_has_no_duplicate_settings(card_text):
    keys = [re.match(r"^\s*([\w:]+)\s*=", ln.split("!")[0]).group(1)
            for ln in card_text.splitlines()
            if re.match(r"^\s*[\w:]+\s*=", ln.split("!")[0])]
    assert len(keys) == len(set(keys)), "a setting appears twice; the last one wins silently"


def test_card_value_ignores_comments():
    assert card_value("! Random:seed = 1\nRandom:seed = 2 ! x", "Random:seed") == "2"
    assert card_value("Foo = 1", "Bar") is None
