import math
import sympy
from sympy.physics.quantum import TensorProduct

from _gate_common import *

class SymPyGateOperation(IGateOperation):
    _gate = build_gate_class(
            'SymPyGate',
            sympy.Matrix,
            sympy.sqrt(2) / 2,
            lambda a, b: a + b * sympy.I
    )
    _identity = staticmethod(sympy.eye)
    _kron = staticmethod(TensorProduct)
    _matmul = staticmethod(sympy.Matrix.multiply)

    @staticmethod
    def _multiply(mat1, mat2):
        return mat1.multiply_elementwise(mat2)

    def __init__(self, n_bits, data):
        IGateOperation.__init__(self, n_bits, data)


class Qubit(SymPyGateOperation, IQubitOperation):
    @staticmethod
    def _zeros(n):
        return sympy.zeros(n, 1)

    @staticmethod
    def _abssq(v):
        return v.conjugate() * v

    @staticmethod
    def _innerproduct(v1, v2):
        return (v1.H * v2)[0, 0]

    _sqrt = staticmethod(sympy.sqrt)

    def _del_idx(self, i):
        self.data.row_del(i)

    def fidelity(self, other):
        return sympy.simplify(IQubitOperation.fidelity(self, other))

    def __init__(self, n_bits, arr=None, measured=None, rng=None):
        self.formal_repr = False
        self.measured = measured
        if arr is None:
            data = self._generate_data(n_bits)
        else:
            if arr.shape == (2**n_bits, 1):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        SymPyGateOperation.__init__(self, n_bits, data)
        IQubitOperation.__init__(self, measured, rng)

    def __repr__(self):
        if not self.formal_repr:
            return str(self)
        return 'Qubit(' + str(self.n_bits) + ', arr=\n' + str(self.data) + ', measured=' + bin(self.measured) + ')'

    def __str__(self):
        def qubit_format(i):
            if isinstance(self.data[i], sympy.Add):
                fmt = '({})|{:0%db}>' % self.n_bits
            else:
                fmt = '{}|{:0%db}>' % self.n_bits
            return fmt.format(self.data[i], i)

        fmtc = '{:0%db}' % self.n_bits
        return ' + '.join(qubit_format(i) for i in range(2**self.n_bits) if self.data[i] != 0) + ' Measured: ' + fmtc.format(self.measured)


class Unitary(SymPyGateOperation):
    def __init__(self, n_bits, arr=None):
        self.n_bits = n_bits
        if arr is None:
            data = sympy.eye(2**n_bits)
            data[0, 0] = 1
        else:
            if arr.shape == (2**n_bits, 2**n_bits):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        SymPyGateOperation.__init__(self, n_bits, data)

    def __repr__(self):
        return 'Unitary(' + str(self.n_bits) + ', arr=' + repr(self.data) + ')'

    def __str__(self):
        return str(self.data)
