import numpy as np
import matplotlib.pyplot as plt

def dagger(matrix):
    """Return the Hermitian conjugate (dagger) of a matrix."""
    return np.transpose(np.conj(matrix))

def commutator(matrix_1, matrix_2, sign):
    """Compute the commutator or anti-commutator of two matrices.
    sign=0: [A,B]; sign=1: {A,B}
    """
    if sign == 0:
        return matrix_1 @ matrix_2 - matrix_2 @ matrix_1
    else:
        return matrix_1 @ matrix_2 + matrix_2 @ matrix_1

def jump_probability(psi, L, dt):
    """Calculate the probability of a quantum jump for a given state."""
    prob = float(np.real(dt * np.vdot(psi, L.conj().T @ L @ psi)))
    return prob

def plot_excitation_density(time, n_e, title_text, compare):
    plt.figure(figsize=(8,6))
    if len(n_e) > 100:
        plt.plot(time, n_e, label=r'$n_{e}(t)$')
    elif compare == 1:
        plt.plot(time, n_e[0], label=r'$n_{e,1}(t)$ - exact')
        plt.plot(time, n_e[1], label=r'$n_{e,1}(t)$ - numerics')
        plt.plot(time, n_e[2], label=r'Stationary value $n_{e}(\infty)$')
    else:
        plt.plot(time, n_e[0], label=r'$n_{e}(t)$ for trajectory 1')
        plt.plot(time, n_e[1], label=r'$n_{e}(t)$ for trajectory 2')
    plt.ylim(0, 1)
    plt.title(title_text)
    plt.xlabel('Time')
    plt.ylabel(r'$n_{e}(t)$')
    plt.legend()
    plt.show()