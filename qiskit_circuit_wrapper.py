import sys
import sympy
from qiskit import QuantumCircuit

class CircuitWrapper:
    def __init__(self, qc):
        if not isinstance(qc, QuantumCircuit):
            print('WARNING: qc is not QuantumCircuit. It may not work.', file=sys.stderr)
        self.qc = qc

    def apply_gate(self, gate):
        if isinstance(gate, sympy.physics.quantum.gate.IdentityGate):
            pass
        elif isinstance(gate, sympy.physics.quantum.gate.HadamardGate):
            self.qc.h(int(gate.args[0]))
        elif isinstance(gate, sympy.physics.quantum.gate.XGate):
            self.qc.x(int(gate.args[0]))
        elif isinstance(gate, sympy.physics.quantum.gate.YGate):
            self.qc.y(int(gate.args[0]))
        elif isinstance(gate, sympy.physics.quantum.gate.ZGate):
            self.qc.z(int(gate.args[0]))
        elif isinstance(gate, sympy.physics.quantum.gate.PhaseGate):
            self.qc.s(int(gate.args[0]))
        elif isinstance(gate, sympy.physics.quantum.gate.TGate):
            self.qc.t(int(gate.args[0]))
        elif isinstance(gate, sympy.physics.quantum.gate.CNotGate):
            self.qc.cx(int(gate.args[0]), int(gate.args[1]))
        elif isinstance(gate, sympy.physics.quantum.gate.SwapGate):
            self.qc.cx(int(gate.args[0]), int(gate.args[1]))
            self.qc.cx(int(gate.args[1]), int(gate.args[0]))
            self.qc.cx(int(gate.args[0]), int(gate.args[1]))
        else:
            assert False, '{}はゲートじゃないか、対応してないゲートです'.format(repr(gate))

    def apply(self, sympy_expr):
        if isinstance(sympy_expr, sympy.mul.Mul):
            for expr in reversed(sympy_expr.args):
                self.apply(expr)
        elif isinstance(sympy_expr, sympy.physics.quantum.gate.Gate):
            self.apply_gate(sympy_expr)

    def __rmul__(self, sympy_expr):
        self.apply(sympy_expr)
