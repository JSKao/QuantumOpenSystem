![CI](https://github.com/JSKao/QuantumOpenSystem/actions/workflows/pytest.yml/badge.svg)
# QuantumOpenSystem

A modular Python project for simulating open quantum systems using the Lindblad master equation and Quantum Jump Monte Carlo (QJMC) method.  
Includes both analytical and numerical solutions for single and two-atom systems, with clear code structure and Jupyter notebook demonstrations.

---

## Features

- **Lindblad & QJMC simulation:** Supports one- and two-qubit open quantum systems.
- **Modular codebase:** All core logic is in `src/` modules (`operators.py`, `simulation.py`, `quantum_system.py`, `utils.py`).
- **Jupyter notebook demo:** Step-by-step workflow in `notebooks/QJMC.ipynb`.
- **Unit tests:** All major modules are covered by tests in `tests/`.
- **Easy to extend:** Add more atoms, observables, or dissipation models.

---

## Installation

1. Clone this repository:

    ```bash
    git clone https://github.com/JSKao/QuantumOpenSystem.git
    cd QuantumOpenSystem
    ```

2. (Recommended) Create a virtual environment:

    ```bash
    python -m venv venv
    source venv/bin/activate   # On Windows: venv\Scripts\activate
    ```

3. Install dependencies:

    ```bash
    pip install -r requirements.txt
    ```

---

## Usage

- **Run the main notebook:**

    Open `notebooks/QJMC.ipynb` in JupyterLab or VS Code and execute the cells in order.  
    The notebook demonstrates both the theory and the code for QJMC and Lindblad simulations.

- **Run unit tests:**

    ```bash
    pytest tests/
    ```

---

## Project Structure

```
QuantumOpenSystem/
├── src/
│   ├── operators.py         # Quantum operators (Pauli, raising/lowering, etc.)
│   ├── simulation.py        # QJMC and analytical solution functions
│   ├── quantum_system.py    # QuantumSystem class (Hamiltonian, Lindblad, propagator)
│   └── utils.py             # Plotting, commutator, dagger, etc.
├── notebooks/
│   └── QJMC.ipynb           # Main demonstration notebook
├── tests/
│   ├── test_operators.py
│   ├── test_simulation.py
│   └── test_utils.py
├── requirements.txt
└── README.md
```

---

## Example Results

You can find example plots and results in the notebook, such as:

- Excitation density evolution for single and two-atom systems
- Comparison between QJMC numerics and analytical Lindblad solutions
- Visualization of dark states and collective effects

---

## References

- Mølmer, K., Dalibard, J., & Castin, Y. Monte Carlo wave-function method in quantum optics. Journal of the Optical Society of America B Vol. 10, Issue 3, pp. 524-538 (1993).
- Andrew J Daley, Quantum trajectories and open many-body quantum systems. Advances in Physics. Vol. 63, 2, 77-149 (2014)

---

## How to Extend

- To simulate more atoms: generalize the operators and initial states in `operators.py` and `simulation.py`.
- To measure other observables: add new measurement functions in `utils.py`.
- To implement more complex dissipation: modify the Lindblad operators in `quantum_system.py`.

---

## License

MIT License (or your preferred license)
