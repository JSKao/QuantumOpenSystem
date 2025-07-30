import numpy as np
from tqdm import tqdm
from sympy import Matrix
from utils import jump_probability, dagger, commutator
from operators import S_x, S_minus, phi_basis_vecs, rho_basis_ops, id, S_plus

def N_trajectory_twoLevelSystem(N, init_state, total_time, dt, jump_dt, QuantumSystem):
    t_length = int(total_time/dt)
    L = QuantumSystem.L
    U_eff = QuantumSystem.U_eff
    gamma = QuantumSystem.gamma

    N_ne_values = []
    N_state_vectors = []
    for n in tqdm(range(N)):
        time_values = []
        state_vectors = []
        ne_values = []
        current_state = init_state
        current_time = 0
        while current_time < total_time:
            current_state_1 = U_eff @ current_state
            dp = jump_probability(current_state_1, L, dt)
            r = np.random.rand()
            if r < dp:
                time_values.append(current_time)
                current_state = np.sqrt(gamma) * (L @ current_state)
                current_state = current_state / np.sqrt(dp/dt)
                state_vectors.append(current_state)
                ne = np.abs(current_state[1][0])**2
                ne_values.append(ne)
                current_time += jump_dt * dt
            else:
                time_values.append(current_time)
                current_state = current_state_1 / np.sqrt(1-dp)
                state_vectors.append(current_state)
                ne = np.abs(current_state[1][0])**2
                ne_values.append(ne)
                current_time += dt
        N_state_vectors.append(np.array(state_vectors))
        N_ne_values.append(np.array(ne_values))
    print("Simulation complete")
    time_values = time_values[:t_length]
    for n in range(N):
        N_ne_values[n] = N_ne_values[n][:t_length]
    N_ne_sum = np.array([0+0j]*t_length, dtype=complex)
    for n in range(N):
        N_ne_sum += N_ne_values[n]
    N_ne_avg = N_ne_sum / N
    return N_state_vectors, N_ne_values, N_ne_avg, time_values


def N_trajectory_twoAtomSystem(N, init_state, total_time, dt, jump_dt, QuantumSystem):
    
    t_length = int(total_time/dt) 
    gamma = QuantumSystem.gamma
    L = QuantumSystem.L
    U_eff = QuantumSystem.U_eff

    N_ne_values = [[],[]]
    N_ne_avg = [[],[]]
    N_state_vectors = []
    for n in tqdm(range(N)):
        
        # start a new trajectory simulation        
        time_values = []
        state_vectors = []
        ne_values = [[],[]]
    
        current_state = init_state
        current_time = 0
        while current_time < total_time:
                  
            # propagation
            current_state_1 = U_eff @ current_state                      # get the candidate propagated state
    
            dp = jump_probability(current_state_1, L, dt)
            r = np.random.rand()
            if r < dp:                                                   # jump occurs
                time_values.append(current_time)
                current_state = np.sqrt(gamma) * (L @ current_state)     # propagate by acting jump
                current_state = current_state / np.sqrt(dp/dt)           # normalize the state
                state_vectors.append(current_state)       
                ne1 = np.abs(current_state[2][0])**2 + np.abs(current_state[3][0])**2
                ne_values[0].append(ne1) 
                ne2 = np.abs(current_state[1][0])**2 + np.abs(current_state[3][0])**2
                ne_values[1].append(ne2)
                current_time += jump_dt * dt                             # update time with smaller dt
                           
            else:
                time_values.append(current_time)
                current_state = current_state_1/ np.sqrt(1-dp)           # choose the candidate propagated state
                state_vectors.append(current_state)
                ne1 = np.abs(current_state[2][0])**2 + np.abs(current_state[3][0])**2
                ne_values[0].append(ne1) 
                ne2 = np.abs(current_state[1][0])**2 + np.abs(current_state[3][0])**2
                ne_values[1].append(ne2)
                current_time += dt                                       # update time with regular dt
                     
        N_state_vectors.append(np.array(state_vectors))
        N_ne_values[0].append(np.array(ne_values[0]))
        N_ne_values[1].append(np.array(ne_values[1]))
    
    print("Simulation complete")
    
    time_values = time_values[:t_length]
    for n in range(N):
        N_ne_values[0][n] = N_ne_values[0][n][:t_length]
        N_ne_values[1][n] = N_ne_values[1][n][:t_length]

    N_ne1_sum = np.array([0+0j]*len(N_ne_values[0][0]), dtype=complex)
    N_ne2_sum = np.array([0+0j]*len(N_ne_values[0][0]), dtype=complex)
    for n in range(N):
        N_ne1_sum += N_ne_values[0][n]
        N_ne2_sum += N_ne_values[1][n]
    
    
    N_ne_avg[0] = N_ne1_sum / N
    N_ne_avg[1] = N_ne2_sum / N
    
    
    return N_state_vectors, N_ne_values, N_ne_avg, time_values


def exact_twoLevelSystem(omega, gamma, init_state, time_values):
    ne_infty = 1/ (2+gamma**2/(4*omega**2))
    u0 = np.outer(init_state, init_state).reshape(4)
    # Define the eigenvectors we've solved
    u1 = np.array([1+gamma**2/(4*omega**2), 1j*gamma/(2*omega), -1j*gamma/(2*omega), 1], dtype=complex)
    u2 = np.array([0, 1, 1, 0], dtype=complex)
    u3 = np.array([-1, 1j*8*omega/(gamma+np.sqrt(gamma**2-64*omega**2,dtype=complex)), -1j*8*omega/(gamma+np.sqrt(gamma**2-64*omega**2,dtype=complex)), 1], dtype=complex)
    u4 = np.array([-1, -1j*8*omega/(-gamma+np.sqrt(gamma**2-64*omega**2,dtype=complex)), 1j*8*omega/(-gamma+np.sqrt(gamma**2-64*omega**2,dtype=complex)), 1], dtype=complex)
    u = [u1, u2, u3, u4]
    
    # Define the eigenvalues we've solved
    l1 = 0
    l2 = -gamma/2
    l3 = 1/4 * (-3*gamma - np.sqrt(gamma**2 - 64*omega**2, dtype=complex))
    l4 = 1/4 * (-3*gamma + np.sqrt(gamma**2 - 64*omega**2, dtype=complex))
    l = [l1, l2, l3, l4]
    
    # Solve c1, c2, c3, c4 by Uc = d
    U = np.zeros((4,4), dtype=complex)
    d = np.zeros((4), dtype=complex)
    for m in range(4):
        d[m] = dagger(u[m]) @ u0
        for n in range(4):
            U[m][n] = dagger(u[m]) @ u[n]
    
    c = np.linalg.solve(U,d)
    
    parameters = [c[0], c[1], c[2], c[3], l1, l2, l3, l4, u1, u2, u3, u4]
    # Construct the solution rho(t) of the master equation
    def Rho(parameters, t):
        state = np.array([0, 0, 0, 0], dtype=complex)
        for i in range(4):
            state += parameters[i] * np.exp(parameters[i+4]*t) * parameters[i+8] 
        return state
    
    
    
    # Construct the excitation density ne(t)
    ne_exact = []
    ne_infty_values = []
    for t in time_values:
        norm = Rho(parameters, t)[0] + Rho(parameters, t)[3]
        ne_exact.append(Rho(parameters, t)[3]/norm)
        ne_infty_values.append(ne_infty)
        
    return ne_exact, ne_infty_values


def exact_twoAtomSystem(omega, gamma, init_state, time_values, QuantumSystem):
    # Build the matrix representation of the Lindbladian superoperator
    gamma = omega/3
    L = QuantumSystem.L
    
    def Lind(rho):
        L1 = commutator(S_x(2,1)+S_x(2,2), rho, 0)
        L2 = L @ rho @ dagger(L) - 1/2 * commutator(dagger(L) @ L, rho, 1)
        return -1j * omega * L1 + 2 * gamma * L2

    M2 = np.zeros((16,16),dtype=complex)
    for m in range(16):
        for n in range(16):
            M2[m][n] = np.trace( dagger(rho_basis_ops(2)[m]) @ Lind(rho_basis_ops(2)[n]))
            
    # Calculate eigenvalues and eigenvectors by sympy module
    M2 = Matrix(M2)
    eigenvalues = M2.eigenvals()
    eigenvectors = M2.eigenvects()

    # Store in numpy list
    vectors = []
    values = []
    for i in range(16):
        value = np.squeeze(np.array(eigenvectors[i][0]).astype(complex))
        vector = np.squeeze(np.array(eigenvectors[i][2]).astype(complex))
        values.append(value)
        vectors.append(vector)
    values = np.squeeze(values)
    
    # Set up initial state rho(0)
    phiphi = np.outer(init_state, init_state)
    u0 = phiphi.reshape(16)

    # Solve c1, c2, ..., c16 by Uc = d
    U = np.zeros((16,16), dtype=complex)
    d = np.zeros((16), dtype=complex)
    for m in range(16):
        d[m] = dagger(vectors[m]) @ u0
        for n in range(16):
            U[m][n] = dagger(vectors[m]) @ vectors[n]
    
    c = np.linalg.solve(U,d)

    # Construct the full solution rho(t) and rho_ss
    parameters = []
    parameters.append(c)
    parameters.append(values)
    parameters.append(vectors)
    
    # rho_ss and n_e(infinity)
    rho_stationary = parameters[0][13]*parameters[2][13] + parameters[0][14]*parameters[2][14]
    norm = rho_stationary[0]+rho_stationary[5]+rho_stationary[10]+rho_stationary[15]
    n_e_1_infty = (rho_stationary[10]+rho_stationary[15])/norm
    n_e_2_infty = (rho_stationary[5]+rho_stationary[15])/norm
    
    # rho(t)
    def Rho(parameters, t):
        state = np.zeros((16), dtype=complex)
        for i in range(16):
            state += parameters[0][i] * np.exp(parameters[1][i]*t) * parameters[2][i] 
        return state

    # Construct the excitation density n_e_1(t), n_e_2(t)
    # n_e_1 = tr(|1)(1|_1 . rho) = tr(E3E3.rho) + tr(E4E4.rho) = rho[10] + rho[15]
    # n_e_2 = tr(|1)(1|_2 . rho) = tr(E2E2.rho) + tr(E4E4.rho) = rho[5] + rho[15]
    ne_exact = [[],[]]
    ne_infty_values = [[],[]]
    
    
    for t in time_values:
        ne_exact[0].append( (Rho(parameters, t)[10] + Rho(parameters, t)[15]) )
        ne_exact[1].append( (Rho(parameters, t)[5] + Rho(parameters, t)[15]) )
        ne_infty_values[0].append(n_e_1_infty)
        ne_infty_values[1].append(n_e_2_infty)


    return ne_exact, ne_infty_values