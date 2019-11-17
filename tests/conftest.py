from pytest import fixture
import sys
import os

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))


@fixture
def fixtures_path():
    return os.path.abspath(os.path.join(here, 'fixtures'))


@fixture
def atom_file(fixtures_path):
    return os.path.join(fixtures_path, 'atom_yb173.json')


@fixture
def trans_file(fixtures_path):
    return os.path.join(fixtures_path, 'transitions_1S0.json')
