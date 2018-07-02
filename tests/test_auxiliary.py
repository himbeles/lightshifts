import pytest
import numpy as np
from lightshifts.auxiliary import atom, smart_gen_array, laser_intensity, frequency_from_wavenumber


def test_branching_ratio_LS(atom_file):
    yb = atom.from_json(atom_file)
    state_i = ('6s8s','3S1')
    state_f = ('6s6p','3P2')
    br = yb.branching_ratio_LS(state_i, state_f)
    assert pytest.approx(0.49202592529614686)==float(br)

def test_smart_gen_array():
    xx, yy = smart_gen_array(lambda x: x, 1, 2, sing=[1.5], n=5, eps=1e-2)
    print(xx,yy)
    assert (len(xx)==12)

def test_laser_intensity():
    li = laser_intensity(1, 1e-2)
    assert li == pytest.approx(2/np.pi)

def test_frequency_from_wavenumber():
    f = frequency_from_wavenumber(1)
    print(f)
    assert f == pytest.approx(29979245800)
