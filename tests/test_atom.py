import pytest
import numpy as np
import lightshifts.atom as atom


def test_branching_ratio_LS(atom_file):
    yb = atom.from_json(atom_file)
    state_i = ('6s8s','3S1')
    state_f = ('6s6p','3P2')
    br = yb.branching_ratio_LS(state_i, state_f)
    assert pytest.approx(0.49202592529614686)==float(br)
