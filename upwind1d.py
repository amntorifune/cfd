import numpy as np
import matplotlib.pyplot as plt

cfl = 0.5
gamma = 1.4
dom = (0.0, 1.0)
nc = 2000
nx = nc + 1
dx = (dom[1] - dom[0]) / nc
x = np.linspace(dom[0], dom[1], nx)
a = 1.0
dt = (cfl * dx) / a
ter = 1.25

def flux(ua, ub):
    F = a * (ub + ua) - np.abs(a) * (ub - ua)
    F = F / 2.0
    return F

u = np.zeros(nx)
# for i in range(int(0.25 / dx), int(0.75 / dx) + 1):
#     u[i] = np.cos((i * dx - 0.5) * 4 * np.pi) + 1
for i in range(nx):
    u[i] = np.cos(4 * np.pi * i / nc)
# u[500:1000] = 1
u0 = u.copy()

t = 0
while t < ter:
    un = u.copy()
    for i in range(1, nx - 1):
        u[i] = un[i] - dt * (flux(un[i], un[i+1]) - flux(un[i-1], un[i])) / dx
    u[0] = un[0] - dt * (flux(un[0], un[1]) - flux(un[-2], un[0])) / dx
    u[-1] = u[0]
    t += dt

plt.plot(x, u0, 'b-', linewidth = 1)
plt.plot(x, u, 'r-', linewidth = 1)
plt.ylabel(r'$\phi$')
plt.tick_params(axis = 'x', bottom = False, labelbottom = False)
plt.grid(True)
plt.show() 