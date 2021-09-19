import numpy as np
import matplotlib.pyplot as plt

cfl = 0.5
gamma = 1.4
dom = (0.0, 1.0)
nc = 1000
nx = nc + 1
dx = (dom[1] - dom[0]) / nc
x = np.linspace(dom[0], dom[1], nx)
mid = int(nc / 2)

def flux(r, u, E, p, gamma):
    f0 = r * u
    f1 = r * u ** 2 + p
    f2 = u * (r * E + p)
    return np.array([f0, f1, f2])

def q2f(q, gamma):
    r = q[0]
    u = q[1] / r
    E = q[2] / r
    p = (gamma - 1.0) * r * (E - 0.5 * u ** 2)
    return flux(r, u, E, p, gamma)

def lwflux(q, dx, gamma, a, nx, dt):
    r = q[0]
    u = q[1] / r
    E = q[2] / r
    p = (gamma - 1.0) * r * (E - 0.5 * u ** 2)
    F = flux(r, u, E, p, gamma)

    qhalf = 0.5 * (q[:, 1:] + q[:, :-1]) - dt * (F[:, 1:] - F[:, :-1]) / (2.0 * dx)
    Fhalf = q2f(qhalf, gamma)
    return (Fhalf[:, 1:] - Fhalf[:, :-1])

r0 = np.zeros(nx)
u0 = np.zeros(nx)
p0 = np.zeros(nx)

r0[:mid] = 1.0
r0[mid:] = 0.125
u0[:mid] = 0.0
u0[mid:] = 0.0
p0[:mid] = 1.0
p0[mid:] = 0.1

E0 = p0 / ((gamma - 1.0) * r0) + 0.5 * u0 ** 2
a0 = np.sqrt(gamma * p0 / r0)
q = np.array([r0, r0 * u0, r0 * E0])
ter = 0.05
t = 0
it = 0
dt = cfl * dx / max(np.abs(u0) + a0)
a = a0

while t < ter:
    q0 = q.copy()
    dF = lwflux(q0, dx, gamma, a, nx, dt)
    q[:, 1:-1] = q0[:, 1:-1] - dF * dt / dx
    q[:, 0] = q0[:, 0]
    q[:, -1] = q0[:, -1]

    r = q[0]
    u = q[1] / r
    E = q[2] / r
    p = (gamma - 1.0) * r * (E - 0.5 * u ** 2)
    a = np.sqrt(gamma * p / r)
    if min(p) < 0:
        print('Negative pressure!')
    
    dt = cfl * dx / max(np.abs(u) + a)
    t += dt
    it += 1

r = q[0]
u = q[1] / r
E = q[2] / r
p = (gamma - 1.0) * r * (E - 0.5 * u ** 2)

fig, ax = plt.subplots(nrows = 4, ncols = 1)
plt.subplot(4, 1, 1)
plt.plot(x, r, 'k-', linewidth = 1)
plt.ylabel(r'$\rho$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(x, u, 'r-', linewidth = 1)
plt.ylabel('$u$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(x, p, 'b-', linewidth = 1)
plt.ylabel('$p$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(x, E, 'g-', linewidth = 1)
plt.ylabel('$E$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
plt.grid(True)

plt.show()