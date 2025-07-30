import numpy as np
from scipy.linalg import expm
from utils import dagger

class QuantumSystem:
    def __init__(self, H, gamma, L, dt):
        self.H = H                                     # Hamiltonian
        self.gamma = gamma                             # Decay rate
        self.L = L                                     # Jump operator
        self.H_eff = H - 1j/2 * gamma * dagger(L) @ L  # Effective Hamiltonian
        self.U_eff = expm(-1j * self.H_eff * dt)       # Propagator