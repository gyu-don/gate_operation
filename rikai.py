import numpy as np
import math


class _gate:
    __sqrt2_inv = math.sqrt(2.0) / 2.0
    h = np.array([[1, 1], [1, -1]], dtype=float) * __sqrt2_inv
    i = np.array([[1, 0], [0, 1]], dtype=float)
    x = np.array([[0, 1], [1, 0]], dtype=float)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=float)
    s = np.array([[1, 0], [0, 1j]], dtype=complex)
    t = np.array([[1, 0], [1, complex(__sqrt2_inv, __sqrt2_inv)]], dtype=complex)


# とりあえず速度は度外視でやってみる。
class NpQubit:
    def __init__(self, n_bits):
        self.n_bits = n_bits
        self.v = np.zeros((2**n_bits, 1), dtype=complex)
        self.v[0, 0] = 1

    def _make_single_gate_mat(self, gate, i):
        j = self.n_bits - i - 1
        if j > 0:
            mat = np.kron(np.identity(2**j), gate)
        else:
            mat = gate
        if i > 0:
            mat = np.kron(mat, np.identity(2**i))
        return mat

    def h(self, i):
        self.v = self._make_single_gate_mat(_gate.h, i) @ self.v
