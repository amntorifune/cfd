import numpy as np
import matplotlib.pyplot as plt
import math

cfl = 1.0
gamma = 1.4
dom = (0.0, 1.0)
nc = 1000
nx = nc + 1
dx = (dom[1] - dom[0]) / nc
x = np.linspace(dom[0], dom[1], nx)
a = 1.0
dt = (cfl * dx) / a
ter = 0.3

def lmtr(r):
    return max(0, min(2 * r, 1.0), min(r, 2.0))

def calcr(ua, ub, uc):
    # r = (ub - ua) / (uc - ub)
    e = 1e-9
    if ub - ua < 0:
        nom = ub - ua - e
    else:
        nom = ub - ua + e
    if uc - ub < 0:
        denom = uc - ub - e
    else:
        denom = uc - ub + e
    return nom / denom

def interpol(ua, ub, uc):
    r = calcr(ua, ub, uc)
    return ub + 0.5 * lmtr(r) * (uc - ub)

def flux(ua, ub, uc, ud):
    r = calcr(ua, ub, uc)
    uL = ub + 0.5 * lmtr(r) * (uc - ub)
    r = calcr(ub, uc, ud)
    uR = uc - 0.5 * lmtr(r) * (ud - uc)
    F = (a * (uR + uL) - np.abs(a) * (uR - uL)) / 2.0
    return F

def ffunc(ut):
    ft = np.zeros(nx)
    for i in range(2, nx - 2):
        dF = flux(ut[i-1], ut[i], ut[i+1], ut[i+2]) - flux(ut[i-2], ut[i-1], ut[i], ut[i+1])
        ft[i] = -dF / dx
    dF = flux(ut[-2], ut[0], ut[1], ut[2]) - flux(ut[-3], ut[-2], ut[0], ut[1])
    ft[0] = -dF / dx
    dF = flux(ut[0], ut[1], ut[2], ut[3]) - flux(ut[-2], ut[0], ut[1], ut[2])
    ft[1] = -dF / dx
    dF = flux(ut[-3], ut[-2], ut[0], ut[1]) - flux(ut[-4], ut[-3], ut[2], ut[0])
    ft[-2] = -dF / dx
    ft[-1] = ft[0]
    return ft

def rk4(un):
    f0 = ffunc(un)
    ut = un + 0.5 * dt * f0
    f1 = ffunc(ut)
    ut = un + 0.5 * dt * f1
    f2 = ffunc(ut)
    ut = un + dt * f2
    f3 = ffunc(ut)
    ut = un + dt * (f0 + 2 * f1 + 2 * f2 + f3) / 6.0
    return ut

u = np.zeros(nx)
# for i in range(int(0.25 / dx), int(0.75 / dx) + 1):
#     u[i] = np.cos((i * dx - 0.5) * 4 * np.pi) + 1
# for i in range(nx):
#     u[i] = np.cos(4 * np.pi * i / nc)
u[int(0.25 / dx): int(0.5 / dx) + 1] = 1
u0 = u.copy()

t = 0
while t < ter:
    u = rk4(u)
    t += dt

def index(i):
    j = int((i * dx - ter * a) / dx)
    if j < 0:
        j = j + nx - 1
    return j

ua = np.zeros(nx)
for i in range(0, nx - 1):
    ua[i] = u0[index(i)]
ua[-1] = ua[0]

plt.plot(x, u0, 'b-', linewidth = 1)
plt.plot(x, ua, 'g-', linewidth = 1)
plt.plot(x, u, 'r-', linewidth = 1)
plt.ylabel(r'$\phi$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
# plt.grid(True)
plt.show() 