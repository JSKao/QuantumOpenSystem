import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import numpy as np
from simulation import (
    N_trajectory_twoLevelSystem,
    N_trajectory_twoAtomSystem,
    exact_twoLevelSystem,
    exact_twoAtomSystem,
)
from operators import S_x, S_minus, phi_basis_vecs
from quantum_system import QuantumSystem

def test_exact_twoLevelSystem_runs():
    omega = 1.0
    gamma = 0.5
    psi0 = phi_basis_vecs(1)[0]
    tlist = np.linspace(0, 1, 10)
    ne, ne_inf = exact_twoLevelSystem(omega, gamma, psi0, tlist)
    assert isinstance(ne, (np.ndarray, list))
    assert len(ne) == len(tlist)

def test_N_trajectory_twoLevelSystem_runs():
    omega = 1.0
    gamma = 0.5
    H = omega * S_x(1,1)
    L = S_minus(1,1)
    dt = 0.01
    psi0 = phi_basis_vecs(1)[0]
    system = QuantumSystem(H, gamma, L, dt)
    states, ne, ne_avg, tlist = N_trajectory_twoLevelSystem(
        3, psi0, 0.1, dt, 0.01, system
    )
    assert len(tlist) > 0
    assert len(ne_avg) == len(tlist)

def test_exact_twoAtomSystem_runs():
    omega = 1.0
    gamma = 0.5
    psi0 = phi_basis_vecs(2)[0]
    tlist = np.linspace(0, 1, 5)
    # 需要一個 dummy QuantumSystem 物件
    H = omega * (S_x(2,1) + S_x(2,2))
    L = (S_minus(2,1) + S_minus(2,2)) / np.sqrt(2)
    dt = 0.01
    system = QuantumSystem(H, gamma, L, dt)
    ne, ne_inf = exact_twoAtomSystem(omega, gamma, psi0, tlist, system)
    assert isinstance(ne, (np.ndarray, list))
    assert len(ne[0]) == len(tlist)

def test_N_trajectory_twoAtomSystem_runs():
    omega = 1.0
    gamma = 0.5
    H = omega * (S_x(2,1) + S_x(2,2))
    L = (S_minus(2,1) + S_minus(2,2)) / np.sqrt(2)
    dt = 0.01
    psi0 = phi_basis_vecs(2)[0]
    system = QuantumSystem(H, gamma, L, dt)
    states, ne, ne_avg, tlist = N_trajectory_twoAtomSystem(
        2, psi0, 0.1, dt, 0.01, system
    )
    assert len(tlist) > 0
    assert len(ne_avg[0]) == len(tlist)