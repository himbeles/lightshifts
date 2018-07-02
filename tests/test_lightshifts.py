import pytest
from lightshifts import lightshift_solver


def test_lightshift_solver_scalar_lightshift(atom_file, trans_file):
    ls = lightshift_solver(atom_file, trans_file)
    ls_scalar = ls.scalar_lightshift(lam=670e-9, laser_intensity=2)
    assert pytest.approx(-19.2007990788261)==float(ls_scalar)


def test_lightshift_solver_lightshifts(atom_file, trans_file):
    ls = lightshift_solver(atom_file, trans_file)
    ls_all = ls.lightshifts(lam=555.7e-9, q=1, mFi=-5/2, laser_intensity=2)
    a = [float(x) for x in ls_all]
    assert pytest.approx([583.819297319377, 16.7205709924352,
        -1.18364597616775])==a


def test_lightshift_solver_total_lightshift(atom_file, trans_file):
    ls = lightshift_solver(atom_file, trans_file)
    ls_tot = ls.total_lightshift(lam=555.7e-9, q=1, mFi=-5/2, laser_intensity=2)
    assert pytest.approx(599.356222335644)==float(ls_tot)
