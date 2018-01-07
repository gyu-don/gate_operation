import math
import numpy as np

from _gate_common import *

class NumPyGateOperation(IGateOperation):
    _gate = build_gate_class('NumPyGate', np.array, math.sqrt(2.0) / 2.0, complex)
    _identity = staticmethod(np.identity)
    _kron = staticmethod(np.kron)
    _matmul = staticmethod(np.dot)
    _inplace_multiply = staticmethod(np.ndarray.__imul__)
    _multiply = staticmethod(np.multiply)

    def __init__(self, n_bits, data):
        IGateOperation.__init__(self, n_bits, data)


class Qubit(NumPyGateOperation, IQubitOperation):
    @staticmethod
    def _zeros(n):
        return np.zeros((n, 1), dtype=complex)

    @staticmethod
    def _abssq(v):
        return np.abs(v)**2

    @staticmethod
    def _innerproduct(v1, v2):
        return np.sum(v1.conj() * v2)

    _sqrt = staticmethod(math.sqrt)

    def _del_idx(self, i):
        self.data = np.delete(self.data, i, 0)

    def __init__(self, n_bits, arr=None, measured=None, rng=None):
        if arr is None:
            data = self._generate_data(n_bits)
        else:
            if arr.shape == (2**n_bits, 1):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        NumPyGateOperation.__init__(self, n_bits, data)
        IQubitOperation.__init__(self, measured, rng)

    def __repr__(self):
        return 'Qubit(' + str(self.n_bits) + ', arr=\n' + str(self.data) + ', measured=' + bin(self.measured) + ')'

    def __str__(self):
        fmt = '{}|{:0%db}>' % self.n_bits
        fmtc = '{:0%db}' % self.n_bits
        return ' + '.join(fmt.format(self.data[i], i) for i in range(2**self.n_bits) if abs(self.data[i]) > 0.00001) + ' Measured: ' + fmtc.format(self.measured)


class Unitary(NumPyGateOperation):
    def __init__(self, n_bits, arr=None, rng=None):
        if arr is None:
            data = np.identity(2**n_bits, dtype=complex)
            data[0, 0] = 1
        else:
            if arr.shape == (2**n_bits, 2**n_bits):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        NumPyGateOperation.__init__(self, n_bits, data)

    def __repr__(self):
        return 'Unitary(' + str(self.n_bits) + ', arr=\n' + str(self.data) + ')'

    def __str__(self):
        return str(self.data)
