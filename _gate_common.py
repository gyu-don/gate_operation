import itertools
import random

class FakeRandom:
    def __init__(self, random_sequence, cycle=False):
        if cycle:
            self.seq = itertools.cycle(iter(random_sequence))
        else:
            self.seq = iter(random_sequence)

    def random(self):
        try:
            return next(self.seq)
        except StopIteration:
            raise ValueError('Given random sequence is over.')

class GateOperation:
    # Shall be set by subclass.
    _gate = None

    @staticmethod
    def _identity(n):
        'Returns an identity matrix of n rows and n colums.'
        # Shall be overrided by subclass.
        raise NotImplemented

    @staticmethod
    def _kron(cls, mat1, mat2):
        'Returns a Kronecker product of mat1 and mat2.'
        # Shall be overrided by subclass.
        raise NotImplemented

    @staticmethod
    def _matmul(cls, mat1, mat2):
        'Returns a matrix product of mat1 and mat2'
        # Shall be overrided by subclass.
        raise NotImplemented

    @staticmethod
    def _inplace_matmul(cls, mat1, mat2):
        'Do an inplace matmul operation, `mat1 = matmul(mat1, mat2)`.'
        mat1 = self._matmul(mat1, mat2)

    def __init__(self, n_bits, data, rng=None):
        '''Initialize self.
        Args:
            n_bits (int): The number of bits.
            data:         The matrix.
            rng:          The random number generator.
        '''
        self.n_bits = n_bits
        self.data = data
        if rng:
            self.rng = rng
        else:
            self.rng = random.Random()

    # Gates
    def u(self, gate_operation):
        self.data = self._matmul(gate_operation.data, self.data)
        return self

    def i(self, i):
        return self

    def h(self, i):
        return self.apply_gate(self._gate.h, i)

    def x(self, i):
        return self.apply_gate(self._gate.x, i)

    def y(self, i):
        return self.apply_gate(self._gate.y, i)

    def z(self, i):
        return self.apply_gate(self._gate.z, i)

    def s(self, i):
        return self.apply_gate(self._gate.s, i)

    def t(self, i):
        return self.apply_gate(self._gate.t, i)

    def s_dag(self, i):
        return self.apply_gate(self._gate.s_dag, i)

    def t_dag(self, i):
        return self.apply_gate(self._gate.t_dag, i)

    def cu(self, c, gate_operation):
        _gate = self._gate
        if c == i:
            raise ValueError('Control bit and operation bit shall be different')
        mat = self._make_single_gate_mat(_gate._one, c)
        mat = self._matmul(gate_operation.data, mat)
        mat += self._make_single_gate_mat(_gate._zero, c)
        self.data = self._matmul(mat, self.data)
        return self

    def ci(self, c, i):
        return self

    def cx(self, c, i):
        return self.apply_cgate(self._gate.x, c, i)

    def cy(self, c, i):
        return self.apply_cgate(self._gate.y, c, i)

    def cz(self, c, i):
        return self.apply_cgate(self._gate.z, c, i)

    def cs(self, c, i):
        return self.apply_cgate(self._gate.s, c, i)

    def ct(self, c, i):
        return self.apply_cgate(self._gate.t, c, i)

    def ch(self, c, i):
        return self.apply_cgate(self._gate.h, c, i)

    def cs_dag(self, c, i):
        return self.apply_cgate(self._gate.s_dag, c, i)

    def ct_dag(self, c, i):
        return self.apply_cgate(self._gate.t_dag, c, i)

    cnot = cx

    def swap(self, i, j):
        return self.cx(i, j).cx(j, i).cx(i, j)

    def ccu(self, c1, c2, gate_operation):
        _gate = self._gate
        if c1 == c2:
            return self.cu(c1, gate_operation)
        matc = self._make_single_gate_mat(_gate._one, c1)
        self._inplace_multiply(matc, self._make_single_gate_mat(_gate._one, c2))
        mat = self._matmul(gate_operation.data, matc)
        mat += self._identity(2**self.n_bits) - matc
        self.data = self._matmul(mat, self.data)
        return self

    def ccnot(self, c1, c2, i):
        return self.apply_ccgate(self._gate.x, c1, c2, i)

    toffoli = ccnot
    ccx = ccnot

    # Utility functions
    def _bit_indices(self, bit_no, bit_value=0):
        stop = 2**self.n_bits
        longstep = 2**(self.n_bits - bit_no - 1)
        shortstep = 2**bit_no
        i = bit_value * longstep
        for _ in range(shortstep):
            for _ in range(longstep):
                yield i
                i += 1
            i += longstep

    def _make_single_gate_mat(self, gate, i):
        kron = self._kron
        identity = self._identity
        if i > 0:
            mat = kron(identity(2**i), gate)
        else:
            mat = gate.copy()
        j = self.n_bits - i - 1
        if j > 0:
            mat = kron(mat, identity(2**j))
        return mat

    def apply_gate(self, gate, i):
        self.data = self._matmul(self._make_single_gate_mat(gate, i), self.data)
        return self

    def apply_cgate(self, gate, c, i):
        _gate = self._gate
        if c == i:
            raise ValueError('Control bit and operation bit shall be different')
        mat = self._make_single_gate_mat(gate, i)
        mat = self._matmul(mat, self._make_single_gate_mat(_gate._one, c))
        mat += self._make_single_gate_mat(_gate._zero, c)
        self.data = self._matmul(mat, self.data)
        return self

    def apply_ccgate(self, gate, c1, c2, i):
        _gate = self._gate
        if c1 == i or c2 == i:
            raise ValueError('Control bit and operation bit shall be different')
        if c1 == c2:
            return self.apply_cgate(gate, c1, i)
        mat = self._make_single_gate_mat(gate, i)
        matc = self._make_single_gate_mat(_gate._one, c1)
        self._inplace_multiply(matc, self._make_single_gate_mat(_gate._one, c2))
        mat = self._matmul(mat, matc)
        mat += self._identity(2**self.n_bits) - matc
        self.data = self._matmul(mat, self.data)
        return self

    def tensorproduct(self, other, rng=None):
        return type(self)(self.n_bits + other.n_bits, self._kron(self.data, other.data), rng)


class IQubitOperation:
    @staticmethod
    def _zeros(n):
        raise NotImplemented

    @staticmethod
    def _abssq(val):
        raise NotImplemented

    @staticmethod
    def _sqrt(val):
        raise NotImplemented

    def _del_idx(self, i):
        raise NotImplemented

    @classmethod
    def _generate_data(cls, n_bits):
        data = cls._zeros(2**n_bits)
        data[0, 0] = 1
        return data

    def _measure_bit(self, i):
        abssq = self._abssq
        bi = self._bit_indices
        normsq = 0
        d = self.data
        for j in bi(i, 0):
            normsq += abssq(d[j])
        r = self.rng.random()
        if r < normsq:
            norm = self._sqrt(normsq)
            for j in bi(i, 0):
                d[j] /= norm
            for j in bi(i, 1):
                d[j] = 0
            return 0
        else:
            norm = self._sqrt(1 - normsq)
            for j in bi(i, 1):
                d[j] /= norm
            for j in bi(i, 0):
                d[j] = 0
            return 1

    def measure(self, i):
        b = self._measure_bit(i)
        j = self.n_bits - i - 1
        self.measured = (self.measured ^ (0 << j)) | (1 << j)
        return self

    m = measure

    def append_qubit(self, n):
        self.data = self._kron(self.data, self._generate_data(n))
        self.n_bits += n
        self.measured <<= n
        return self

    def merge_qubit(self, other):
        self.data = self._kron(self.data, other.data)
        self.n_bits += other.n_bits
        self.measured = (self.measured << other.n_bits) | other.measured
        return self

    def drop_qubit(self, i):
        data = self.data
        del_idx = self._del_idx
        b = self._measure_bit(i)
        if b == 1:
            b = 0
        else:
            b = 1
        indices = list(self._bit_indices(i, b))
        indices.sort(reverse=True)
        for j in indices:
            del_idx(j)
        return self


def build_gate_class(name, f_matrix, sqrt2_inv, make_complex):
    d = {}
    d['__sqrt2_inv'] = sqrt2_inv
    d['h'] = f_matrix([[sqrt2_inv, sqrt2_inv], [sqrt2_inv, -sqrt2_inv]])
    d['i'] = f_matrix([[1, 0], [0, 1]])
    d['x'] = f_matrix([[0, 1], [1, 0]])
    d['y'] = f_matrix([[0, make_complex(0, -1)], [make_complex(0, 1), 0]])
    d['z'] = f_matrix([[1, 0], [0, -1]])
    d['s'] = f_matrix([[1, 0], [0, make_complex(0, 1)]])
    d['t'] = f_matrix([[1, 0], [0, make_complex(sqrt2_inv, sqrt2_inv)]])
    d['s_dag'] = f_matrix([[1, 0], [0, make_complex(0, -1)]])
    d['t_dag'] = f_matrix([[1, 0], [0, make_complex(sqrt2_inv, -sqrt2_inv)]])
    # _zero and _one are not unitary but _zero + _one is identity. It is used for Control-U gate.
    d['_zero'] = f_matrix([[1, 0], [0, 0]])
    d['_one'] = f_matrix([[0, 0], [0, 1]])
    return type(name, (object,), d)
