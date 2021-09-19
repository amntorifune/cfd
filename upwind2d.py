import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

cfl = 1.0
domx = domy = (0.0, 1.0)
ncx = ncy = 100
nx = ncx + 1
ny = ncy + 1
x = np.linspace(domx[0], domx[1], nx)
y = np.linspace(domy[0], domy[1], ny)
A = 2.0
B = 0.5
dx = (domx[1] - domx[0]) / ncx
dy = (domy[1] - domx[0]) / ncy
dt = cfl * min(dx, dy) / (np.sqrt(A ** 2 + B ** 2) * np.sqrt(2))
ter = 1

u = np.zeros((nx, ny))
# for i in range(int(0.25 / dx), int(0.5 / dx) + 1):
#     for j in range(int(0.25 / dy), int(0.5 / dy) + 1):
#         u[i, j] = (np.cos((i * dx - 0.375) * 8 * np.pi) + 1) * (np.cos((j * dy - 0.375) * 8 * np.pi) + 1) / 4
for i in range(0, nx):
    for j in range(0, ny):
        u[j, i] = np.cos(4 * np.pi * i / ncx) + np.cos(4 * np.pi * j / ncy)
# u[25:50, 25:50] = 1

X, Y = np.meshgrid(x, y)

def flux(d, a, l, r):
    F = a * (l + r) - np.abs(a) * (r - l)
    F = F / 2.0
    return F
u0 = u.copy()

t = 0.0
while t < ter:
    un = u.copy()
    for j in range(0, ny - 1):
        for i in range(1, nx - 1):
            dF = flux(dx, A, un[j, i], un[j, i + 1]) - flux(dx, A, un[j, i - 1], un[j, i])
            u[j, i] = un[j, i] - dt * dF / dx
        dF = flux(dx, A, un[j, 0], un[j, 1]) - flux(dx, A, un[j, -2], un[j, 0])
        u[j, 0] = un[j, 0] - dt * dF / dx
    u[:, -1] = u[:, 0]
    u[-1, :] = u[0, :]
    un = u.copy()
    for i in range(0, nx - 1):
        for j in range(1, ny - 1):
            dF = flux(dy, B, un[j, i], un[j + 1, i]) - flux(dy, B, un[j - 1, i], un[j, i])
            u[j, i] = un[j, i] - dt * dF / dy
        dF = flux(dy, B, un[0, i], un[1, i]) - flux(dy, B, un[-2, i], un[0, i])
        u[0, i] = un[0, i] - dt * dF / dy
    u[:, -1] = u[:, 0]
    u[-1, :] = u[0, :]
    t += dt

print(t)

fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (13, 5))
plt.subplot(1, 2, 1)
plt.contourf(X, Y, u0, 16, alpha = 0.7, cmap=cm.viridis)  
plt.colorbar()
plt.contour(X, Y, u0, 16, cmap=cm.viridis)  

plt.subplot(1, 2, 2)
plt.contourf(X, Y, u, 16, alpha = 0.7, cmap=cm.viridis)  
plt.colorbar()
plt.contour(X, Y, u, 16, cmap=cm.viridis)

ax[0].set_title('t = 0')
ax[1].set_title('t = %f'%(ter))
plt.show()