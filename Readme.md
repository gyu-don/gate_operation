Library for practice and exersize of quantum information.

Example:

```py
# Quantum teleportation

import np_gate # NumPy backend.
# or
# import sp_gate # SymPy backend. (Performance is worse than np_gate)

alice = np_gate.Qubit(1).h(0).t(0)

# This is simulation. The cloning is possible:)
# |circuit> = |alice>
circuit = alice.clone()

# |circuit> = |alice>|00>
circuit.append_qubit(2)

# |circuit> = (|alice>|00> + |alice>|11>) / sqrt(2)
circuit.h(1).cx(1, 2)

# --*--H-----*--
# --|-----*--|--
# --X-----X--Z--
circuit.cx(0, 1).h(0).cx(1, 2).cz(0, 2)

# Qubit 1 and 0 are measured and ignored.
circuit.drop_qubit(1).drop_qubit(0)

print(alice.fidelity(circuit)) # => Almost 1
```
