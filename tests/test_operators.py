import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import numpy as np
from operators import S_x, S_minus, S_plus, phi_basis_vecs, rho_basis_ops, id

def test_S_x_1qubit():
    sx = S_x(1, 1)
    expected = np.array([[0, 1], [1, 0]], dtype=complex)
    assert np.allclose(sx, expected)

def test_S_minus_1qubit():
    sm = S_minus(1, 1)
    # print(sm)
    expected = np.array([[0, 1], [0, 0]], dtype=complex)
    assert np.allclose(sm, expected)

def test_S_plus_1qubit():
    sp = S_plus(1, 1)
    # print(sp)
    expected = np.array([[0, 0], [1, 0]], dtype=complex)
    assert np.allclose(sp, expected)

def test_phi_basis_vecs_2qubit():
    vecs = phi_basis_vecs(2)
    assert len(vecs) == 4
    assert np.allclose(vecs[0], np.array([[1],[0],[0],[0]], dtype=complex))
    assert np.allclose(vecs[3], np.array([[0],[0],[0],[1]], dtype=complex))

def test_rho_basis_ops_1qubit():
    ops = rho_basis_ops(1)
    assert len(ops) == 4
    # |0><0|
    assert np.allclose(ops[0], np.array([[1,0],[0,0]], dtype=complex))
    # |0><1|
    assert np.allclose(ops[1], np.array([[0,1],[0,0]], dtype=complex))
    # |1><0|
    assert np.allclose(ops[2], np.array([[0,0],[1,0]], dtype=complex))
    # |1><1|
    assert np.allclose(ops[3], np.array([[0,0],[0,1]], dtype=complex))

def test_id_2qubit():
    id2 = id(2)
    assert np.allclose(id2, np.eye(4, dtype=complex))