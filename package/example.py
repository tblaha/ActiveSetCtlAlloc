# https://github.com/chcomin/ctypes-numpy-example/tree/master/simplest
import numpy as np
import ctypes as ct
funcs = ct.CDLL("./bin/solveActiveSet.so")

CA_N_V = 6
CA_N_U = 20
CA_N_C = CA_N_U + CA_N_V

n_v = 4
n_u = 4
theta = ct.c_double(0.002)
cond_bound = ct.c_double(1e8)

JG = np.random.random((CA_N_V, CA_N_U))
Wv = np.random.random(CA_N_V)
Wu = np.random.random(CA_N_U)
up = np.random.random(CA_N_U)
dv = np.random.random(CA_N_V)

A = np.zeros((CA_N_C, CA_N_U))
b = np.zeros((CA_N_C))
gamma = ct.c_double()

c_doublep = ct.POINTER(ct.c_double)
c_intp = ct.POINTER(ct.c_int)
c_int8p = ct.POINTER(ct.c_int8)

funcs.setup_wls.restype = None
funcs.setup_wls(
    n_v, n_u,
    JG.ctypes.data_as(c_doublep),
    Wv.ctypes.data_as(c_doublep),
    Wu.ctypes.data_as(c_doublep),
    up.ctypes.data_as(c_doublep),
    dv.ctypes.data_as(c_doublep),
    theta,
    cond_bound,
    A.ctypes.data_as(c_doublep),
    b.ctypes.data_as(c_doublep),
    ct.byref(gamma),
    )

umin = np.zeros((CA_N_U))
umax = np.ones((CA_N_U))
us = np.zeros((CA_N_U))
Ws = np.zeros((CA_N_U), dtype=np.byte)
updating = ct.c_bool(True)
imax = 10
iter = ct.c_int()
n_free = ct.c_int()
costs = np.zeros((15), )

funcs.setup_wls.restype = ct.c_int8
res = funcs.solveActiveSet_chol(
    A.ctypes.data_as(c_doublep),
    b.ctypes.data_as(c_doublep),
    umin.ctypes.data_as(c_doublep),
    umax.ctypes.data_as(c_doublep),
    us.ctypes.data_as(c_doublep),
    Ws.ctypes.data_as(c_int8p),
    updating,
    imax,
    n_u, n_v,
    ct.byref(iter), ct.byref(n_free),
    costs.ctypes.data_as(c_doublep),
)

print(f'Return: {res}')
print(f'us    : {us[:n_u]}')