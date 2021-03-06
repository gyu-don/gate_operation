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

def single_bit(f):
    def g(self, *args):
        if not args:
            raise TypeError(f.__name__ + '() takes at least one argument (0 given)')
        last_i = -1
        last_ellipsis = False
        for i in args:
            if i is Ellipsis:
                if last_ellipsis:
                    raise ValueError('Consecutive Ellipsis is not allowed.')
                last_ellipsis = True
            else:
                if last_ellipsis:
                    if i <= last_i:
                        raise ValueError('i < j is required when arguments contain i, ..., j')
                    for j in range(last_i + 1, i):
                        f(self, j)
                    last_ellipsis = False
                f(self, i)
                last_i = i
        if last_ellipsis:
            for j in range(last_i + 1, self.n_bits):
                f(self, j)
        return self
    return g

class IGateOperation:
    # Shall be set by subclass.
    _gate = None

    @staticmethod
    def _identity(n):
        'Returns an identity matrix of n rows and n colums.'
        # Shall be overrided by subclass.
        raise NotImplementedError

    @staticmethod
    def _kron(cls, mat1, mat2):
        'Returns a Kronecker product of mat1 and mat2.'
        # Shall be overrided by subclass.
        raise NotImplementedError

    @staticmethod
    def _matmul(cls, mat1, mat2):
        'Returns a matrix product of mat1 and mat2'
        # Shall be overrided by subclass.
        raise NotImplementedError

    @staticmethod
    def _inplace_multiply(mat1, mat2):
        'Do an inplace element-wise multiply operation.'
        raise NotImplementedError

    @staticmethod
    def _multiply(mat1, mat2):
        'Do an element-wise multiply operation.'
        # Shall be overrided by subclass.
        raise NotImplementedError

    def __init__(self, n_bits, data):
        '''Initialize self.
        Args:
            n_bits (int): The number of bits.
            data:         The matrix.
        '''
        self.n_bits = n_bits
        self.data = data

    # Gates
    def u(self, gate_operation):
        self.data = self._matmul(gate_operation.data, self.data)
        return self

    def i(self, i):
        return self

    @single_bit
    def h(self, i):
        return self.apply_gate(self._gate.h, i)

    @single_bit
    def x(self, i):
        return self.apply_gate(self._gate.x, i)

    @single_bit
    def y(self, i):
        return self.apply_gate(self._gate.y, i)

    @single_bit
    def z(self, i):
        return self.apply_gate(self._gate.z, i)

    def rx(self, i, theta):
        return self.apply_gate(self._gate.rx(theta), i)

    def ry(self, i, theta):
        return self.apply_gate(self._gate.ry(theta), i)

    def rz(self, i, theta):
        return self.apply_gate(self._gate.rz(theta), i)

    @single_bit
    def s(self, i):
        return self.apply_gate(self._gate.s, i)

    @single_bit
    def t(self, i):
        return self.apply_gate(self._gate.t, i)

    @single_bit
    def s_dag(self, i):
        return self.apply_gate(self._gate.s_dag, i)

    @single_bit
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
        matc2 = self._make_single_gate_mat(_gate._one, c2)
        try:
            self._inplace_multiply(matc, matc2)
        except NotImplementedError:
            matc = self._multiply(matc, matc2)
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

    def apply(self, operation, **placeholder):
        return operation.apply_to(self, **placeholder)

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
        matc2 = self._make_single_gate_mat(_gate._one, c2)
        try:
            self._inplace_multiply(matc, matc2)
        except NotImplementedError:
            matc = self._multiply(matc, matc2)
        mat = self._matmul(mat, matc)
        mat += self._identity(2**self.n_bits) - matc
        self.data = self._matmul(mat, self.data)
        return self

    def tensorproduct(self, other, *args, **kwargs):
        return type(self)(self.n_bits + other.n_bits, self._kron(self.data, other.data), *args, **kwargs)

    def clone(self, *args, **kwargs):
        return type(self)(self.n_bits, self.data.copy(), *args, **kwargs)


class IQubitOperation:
    @staticmethod
    def _zeros(n):
        raise NotImplementedError

    @staticmethod
    def _abssq(val):
        raise NotImplementedError

    @staticmethod
    def _sqrt(val):
        raise NotImplementedError

    @staticmethod
    def _innerproduct(vec1, vec2):
        raise NotImplementedError

    @classmethod
    def _generate_data(cls, n_bits):
        data = cls._zeros(2**n_bits)
        data[0, 0] = 1
        return data

    def _del_idx(self, i):
        raise NotImplementedError

    def __init__(self, measured=None, rng=None):
        if rng is None:
            rng = random.Random()
        self.rng = rng
        if measured is None:
            measured = 0
        self.measured = measured

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

    @single_bit
    def measure(self, i):
        b = self._measure_bit(i)
        j = self.n_bits - i - 1
        self.measured = (self.measured ^ (0 << j)) | (b << j)
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
        self.n_bits -= 1
        return self

    def ignore_global(self):
        d = self.data
        for i in range(2**self.n_bits - 1):
            abssq = self._abssq(d[i])
            if abssq < 0.00001:
                continue
            absval = self._sqrt(abssq)
            angle = absval / d[i]
            d[i] = absval
            for j in range(i + 1, 2**self.n_bits):
                d[j] *= angle
            break
        return self

    @single_bit
    def reset(self, i):
        if self._measure_bit(i):
            self.x(i)
        return self

    def fidelity(self, other):
        v = self._innerproduct(self.data, other.data)
        return abs(v)


def build_gate_class(name, f_matrix, sqrt2_inv, make_complex, sin, cos, exp):
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
    # Rotations
    d['rx'] = lambda theta: f_matrix(
            [[cos(theta / 2), make_complex(0, -sin(theta / 2))],
             [make_complex(0, -sin(theta / 2)), cos(theta / 2)]])
    d['ry'] = lambda theta: f_matrix(
            [[cos(theta / 2), -sin(theta / 2)],
             [sin(theta / 2), cos(theta / 2)]])
    d['rz'] = lambda theta: f_matrix(
            [[exp(make_complex(0, -theta / 2)), 0],
             [0, exp(make_complex(0, theta / 2))]])
    # _zero and _one are not unitary but _zero + _one is identity. It is used for Control-U gate.
    d['_zero'] = f_matrix([[1, 0], [0, 0]])
    d['_one'] = f_matrix([[0, 0], [0, 1]])
    return type(name, (object,), d)
