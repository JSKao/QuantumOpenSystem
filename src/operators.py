import numpy as np
from functools import reduce

def phi_basis_vecs(N):
    """Generate computational basis vectors for N qubits."""
    dim = 2**N
    return [np.eye(dim, dtype=complex)[:, [i]] for i in range(dim)]

def rho_basis_ops(N):
    """Generate computational basis operators for N qubits."""
    dim = 2**N
    return [np.eye(dim, dtype=complex)[:, [i]] @ np.eye(dim, dtype=complex)[:, [j]].T for i in range(dim) for j in range(dim)]

def id(N):
    """Return the identity matrix for N qubits."""
    dim = 2**N
    return np.eye(dim, dtype=complex)

def S_x(N, n):
    """Construct the Pauli X operator for the n-th atom in an N-atom system."""
    s_x = np.array([[0, 1], [1, 0]], dtype=complex)
    identity = np.eye(2, dtype=complex)
    ops = [s_x if i == n-1 else identity for i in range(N)]
    return reduce(np.kron, ops)

def S_plus(N, n):
    """Construct the Pauli raising operator (σ⁺) for the n-th atom in an N-atom system."""
    s_plus = np.array([[0, 0], [1, 0]], dtype=complex)
    identity = np.eye(2, dtype=complex)
    ops = [s_plus if i == n-1 else identity for i in range(N)]
    return reduce(np.kron, ops)

def S_minus(N, n):
    """Construct the Pauli lowering operator (σ⁻) for the n-th atom in an N-atom system."""
    s_minus = np.array([[0, 1], [0, 0]], dtype=complex)
    identity = np.eye(2, dtype=complex)
    ops = [s_minus if i == n-1 else identity for i in range(N)]
    return reduce(np.kron, ops)