import math
import numpy as np

from _gate_common import *

class NumPyGateOperation(GateOperation):
    _gate = build_gate_class('NumPyGate', np.array, math.sqrt(2.0) / 2.0, complex)
    _identity = staticmethod(np.identity)
    _kron = staticmethod(np.kron)
    _matmul = staticmethod(np.dot)
    _inplace_multiply = staticmethod(np.ndarray.__imul__)

    def __init__(self, n_bits, data, rng=None):
        GateOperation.__init__(self, n_bits, data, rng)


class Qubit(NumPyGateOperation):
    def __init__(self, n_bits, arr=None, measured=0, rng=None):
        self.n_bits = n_bits
        self.measured = measured
        if arr is None:
            data = np.zeros((2**n_bits, 1), dtype=complex)
            data[0, 0] = 1
        else:
            if arr.shape == (2**n_bits, 1):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        NumPyGateOperation.__init__(self, n_bits, data, rng)

    def __repr__(self):
        return 'Qubit(' + str(self.n_bits) + ', arr=\n' + str(self.data) + ', measured=' + bin(self.measured) + ')'

    def __str__(self):
        fmt = '{}|{:0%db}>' % self.n_bits
        fmtc = '{:0%db}' % self.n_bits
        return ' + '.join(fmt.format(self.data[i], i) for i in range(2**self.n_bits) if abs(self.data[i]) > 0.00001) + ' Measured: ' + fmtc.format(self.measured)

    def measure(self, i):
        normsq = 0
        d = self.data
        for j in self._bit_indices(i, 0):
            normsq += np.abs(d[j])**2
        r = self.rng.random()
        if r < normsq:
            norm = math.sqrt(normsq)
            for j in self._bit_indices(i, 0):
                self.data[j] /= norm
            for j in self._bit_indices(i, 1):
                self.data[j] = 0
            self.measured ^= self.measured & (1 << (self.n_bits - i - 1))
            return self
        else:
            norm = math.sqrt(1 - normsq)
            for j in self._bit_indices(i, 1):
                self.data[j] /= norm
            for j in self._bit_indices(i, 0):
                self.data[j] = 0
            self.measured |= 1 << (self.n_bits - i - 1)
            return self

    m = measure


class Unitary(NumPyGateOperation):
    def __init__(self, n_bits, arr=None, rng=None):
        self.n_bits = n_bits
        if arr is None:
            data = np.identity(2**n_bits, dtype=complex)
            data[0, 0] = 1
        else:
            if arr.shape == (2**n_bits, 2**n_bits):
                data = arr
            else:
                raise ValueError('Unexpected length of array is given.')
        NumPyGateOperation.__init__(self, n_bits, data, rng)

    def __repr__(self):
        return 'Unitary(' + str(self.n_bits) + ', arr=\n' + str(self.data) + ')'

    def __str__(self):
        return str(self.data)
