import numpy as np
import matplotlib.pyplot as plt

cfl = 1.0
gamma = 1.4
dom = (0.0, 1.0)
nc = 2000
nx = nc + 1
dx = (dom[1] - dom[0]) / nc
x = np.linspace(dom[0], dom[1], nx)
u = 2.0
dt = (cfl * dx) / u
ter = 0.4

def flux(a, b):
    F = u * (a + b - dt * u * (b - a) / dx) / 2
    return F

q = np.zeros(nx)
# for i in range(int(0.25 / dx), int(0.75 / dx) + 1):
#     q[i] = np.cos((i * dx - 0.5) * 4 * np.pi) + 1
for i in range(nx):
    q[i] = np.cos(4 * np.pi * i / nc)

q0 = q.copy()

t = 0
while t < ter:
    qn = q.copy()
    for i in range(1, nx - 1):
        q[i] = qn[i] - dt * (flux(qn[i], qn[i + 1]) - flux(qn[i - 1], qn[i])) / dx
    q[0] = qn[0] - dt * (flux(qn[0], qn[1]) - flux(qn[-2], qn[0])) / dx
    q[-1] = q[0]
    t += dt

plt.plot(x, q0, 'b-', linewidth = 1)
plt.plot(x, q, 'r-', linewidth = 1)
plt.ylabel(r'$\phi$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
plt.grid(True)
plt.show() 