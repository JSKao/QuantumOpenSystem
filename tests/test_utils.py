import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import numpy as np
from utils import dagger, commutator, jump_probability

def test_dagger():
    A = np.array([[1+2j, 2], [3, 4-1j]])
    A_dag = dagger(A)
    expected = np.array([[1-2j, 3], [2, 4+1j]])
    assert np.allclose(A_dag, expected)

def test_commutator():
    A = np.array([[0, 1], [1, 0]], dtype=complex)
    B = np.array([[0, -1j], [1j, 0]], dtype=complex)
    comm = commutator(A, B, 0)
    expected = np.array([[0.+2.j, 0.+0.j], [0.+0.j, 0.-2.j]])
    assert np.allclose(comm, expected)

def test_jump_probability():
    psi = np.array([[1.0], [0.0]], dtype=complex)
    L = np.array([[0, 1], [0, 0]], dtype=complex)
    dt = 0.1
    prob = jump_probability(psi, L, dt)
    assert np.isclose(prob, 0.0)