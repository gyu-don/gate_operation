class Operation:
    def __init__(self):
        self._operations = []

    def apply_to(self, obj, **placeholder):
        if placeholder:
            return self.apply(self.placeholder(**placeholder))
        else:
            for name, _args, _kwargs in self._operations:
                obj = getattr(obj, name)(*_args, **_kwargs)
            return obj

    def placeholder(self, **placeholder):
        op = Operation()
        for name, _args, _kwargs in self._operations:
            _args = (placeholder.get(x, x) for x in _args)
            _kwargs = {k: placeholder.get(v, v) for k,v in _kwargs.items()}
            getattr(op, name)(*_args, **_kwargs)
        return op

    def __repr__(self):
        chain = ['Operation()']
        for name, _args, _kwargs in self._operations:
            s = name + '('
            if _args:
                s += ', '.join(repr(x) for x in _args)
                if _kwargs:
                    s += ', '
            if _kwargs:
                s += ', '.join(k + '=' + repr(v) for k, v in _kwargs.items())
            s += ')'
            chain.append(s)
        return '.'.join(chain)

    def _append_to_operations(self, name, _args, _kwargs):
        self._operations.append((name, _args, _kwargs))


    def __wellknown_gate(g):
        def f(_self, *args, **kwargs):
            _self._append_to_operations(g.__name__, args, kwargs)
            return _self
        return f

    def __getattr__(self, name):
        def f(*args, **kwargs):
            self._append_to_operations(name, args, kwargs)
            return self
        return f

    # Gates
    @__wellknown_gate
    def u(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def i(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def h(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def x(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def y(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def z(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def s(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def t(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def s_dag(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def t_dag(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cu(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ci(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cx(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cy(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cz(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cs(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ct(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ch(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cs_dag(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ct_dag(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def cnot(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def swap(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ccu(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ccnot(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def toffoli(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def ccx(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def clone(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def measure(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def m(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def append_qubit(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def merge_qubit(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def drop_qubit(self, *args, **kwargs):
        pass

    @__wellknown_gate
    def fidelity(self, *args, **kwargs):
        pass
