import math
import sympy
from sympy.physics.quantum import TensorProduct

from _gate_common import *

class SymPyGateOperation(GateOperation):
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
    def _inplace_multiply(mat1, mat2):
        mat1 = mat1.multiply_elementwise(mat2)

    def __init__(self, n_bits, data, rng=None):
        self.n_bits = n_bits
        self.data = data
        if rng:
            self.rng = rng
        else:
            self.rng = random.Random()


class Qubit(SymPyGateOperation):
    def __init__(self, n_bits, arr=None, measured=0, rng=None):
        self.n_bits = n_bits
        self.measured = measured
        if arr is None:
            data = sympy.zeros(2**n_bits, 1)
            data[0, 0] = 1
        else:
            if arr.shape == (2**n_bits, 1):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        SymPyGateOperation.__init__(self, n_bits, data, rng)

    def __repr__(self):
        return 'Qubit(' + str(self.n_bits) + ', arr=\n' + str(self.data) + ', measured=' + bin(self.measured) + ')'

    def __str__(self):
        def qubit_format(i):
            if isinstance(self.data[i], sympy.Add):
                fmt = '({})|{:0%db}>' % self.n_bits
            else:
                fmt = '{}|{:0%db}>' % self.n_bits
            return fmt.format(self.data[i], i)

        fmtc = '{:0%db}' % self.n_bits
        return ' + '.join(qubit_format(i) for i in range(2**self.n_bits) if abs(self.data[i]) > 0.00001) + ' Measured: ' + fmtc.format(self.measured)

    def measure(self, i):
        normsq = 0
        d = self.data
        for j in self._bit_indices(i, 0):
            normsq += d[j].conjugate() * d[j]
        r = self.rng.random()
        if r < normsq:
            norm = sympy.sqrt(normsq)
            for j in self._bit_indices(i, 0):
                self.data[j] /= norm
            for j in self._bit_indices(i, 1):
                self.data[j] = 0
            self.measured ^= self.measured & (1 << (self.n_bits - i - 1))
            return self
        else:
            norm = sympy.sqrt(1 - normsq)
            for j in self._bit_indices(i, 1):
                self.data[j] /= norm
            for j in self._bit_indices(i, 0):
                self.data[j] = 0
            self.measured |= 1 << (self.n_bits - i - 1)
            return self

    m = measure


class Unitary(SymPyGateOperation):
    def __init__(self, n_bits, arr=None, rng=None):
        self.n_bits = n_bits
        if arr is None:
            data = sympy.eye(2**n_bits)
            data[0, 0] = 1
        else:
            if arr.shape == (2**n_bits, 2**n_bits):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        SymPyGateOperation.__init__(self, n_bits, data, rng)

    def __repr__(self):
        return 'Unitary(' + str(self.n_bits) + ', arr=' + repr(self.data) + ')'

    def __str__(self):
        return str(self.data)
