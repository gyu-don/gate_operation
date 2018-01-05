import itertools
import math
import random
import sympy
from sympy.physics.quantum import TensorProduct

from _gate_common import *


class SymPyGateOperation(GateOperation):
    def __init__(self, n_bits, data, rng=None):
        self.n_bits = n_bits
        self.data = data
        if rng:
            self.rng = rng
        else:
            self.rng = random.Random()

    def _make_single_gate_mat(self, gate, i):
        if i > 0:
            mat = TensorProduct(sympy.eye(2**i), gate)
        else:
            mat = gate
        j = self.n_bits - i - 1
        if j > 0:
            mat = TensorProduct(mat, sympy.eye(2**j))
        return mat

    def apply_gate(self, gate, i):
        self.data = self._make_single_gate_mat(gate, i) * self.data
        return self

    def u(self, gate_operation):
        self.data = gate_operation.data * self.data
        return self

    def i(self, i):
        return self

    def h(self, i):
        return self.apply_gate(_gate.h, i)

    def x(self, i):
        return self.apply_gate(_gate.x, i)

    def y(self, i):
        return self.apply_gate(_gate.y, i)

    def z(self, i):
        return self.apply_gate(_gate.z, i)

    def s(self, i):
        return self.apply_gate(_gate.s, i)

    def t(self, i):
        return self.apply_gate(_gate.t, i)

    def s_dag(self, i):
        return self.apply_gate(_gate.s_dag, i)

    def t_dag(self, i):
        return self.apply_gate(_gate.t_dag, i)

    def apply_cgate(self, gate, c, i):
        if c == i:
            raise ValueError('Control bit and operation bit shall be different')
        mat = self._make_single_gate_mat(gate, i)
        mat = mat * self._make_single_gate_mat(_gate._one, c)
        mat += self._make_single_gate_mat(_gate._zero, c)
        self.data = mat * self.data
        return self

    def cu(self, c, gate_operation):
        if c == i:
            raise ValueError('Control bit and operation bit shall be different')
        mat = self._make_single_gate_mat(_gate._one, c)
        mat = gate_operation.data * mat
        mat += self._make_single_gate_mat(_gate._zero, c)
        self.data = mat * self.data
        return self

    def ci(self, c, i):
        return self

    def cx(self, c, i):
        return self.apply_cgate(_gate.x, c, i)

    def cy(self, c, i):
        return self.apply_cgate(_gate.y, c, i)

    def cz(self, c, i):
        return self.apply_cgate(_gate.z, c, i)

    def cs(self, c, i):
        return self.apply_cgate(_gate.s, c, i)

    def ct(self, c, i):
        return self.apply_cgate(_gate.t, c, i)

    def ch(self, c, i):
        return self.apply_cgate(_gate.h, c, i)

    def cs_dag(self, c, i):
        return self.apply_cgate(_gate.s_dag, c, i)

    def ct_dag(self, c, i):
        return self.apply_cgate(_gate.t_dag, c, i)

    cnot = cx

    def swap(self, i, j):
        return self.cx(i, j).cx(j, i).cx(i, j)

    def apply_ccgate(self, gate, c1, c2, i):
        if c1 == i or c2 == i:
            raise ValueError('Control bit and operation bit shall be different')
        if c1 == c2:
            return self.apply_cgate(gate, c1, i)
        mat = self._make_single_gate_mat(gate, i)
        matc = self._make_single_gate_mat(_gate._one, c1)
        matc *= self._make_single_gate_mat(_gate._one, c2)
        mat = mat * matc
        mat += sympy.eye(2**self.n_bits) - matc
        self.data = mat * self.data
        return self

    def ccu(self, c1, c2, gate_operation):
        if c1 == i or c2 == i:
            raise ValueError('Control bit and operation bit shall be different')
        if c1 == c2:
            return self.cu(c1, gate_operation)
        matc = self._make_single_gate_mat(_gate._one, c1)
        matc *= self._make_single_gate_mat(_gate._one, c2)
        mat = gate_operation.data * matc
        mat += sympy.eye(2**self.n_bits) - matc
        self.data = mat * self.data
        return self

    def ccnot(self, c1, c2, i):
        return self.apply_ccgate(_gate.x, c1, c2, i)

    toffoli = ccnot
    ccx = ccnot


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
