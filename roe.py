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

def roeflux(q, dx, gamma, a, nx):
    r = q[0]
    u = q[1] / r
    E = q[2] / r
    p = (gamma - 1.0) * r * (E - 0.5 * u ** 2)
    H = gamma / (gamma - 1) * p / r + 0.5 * u ** 2
    F = flux(r, u, E, p, gamma)

    Ph = np.zeros([3, nx - 1])
    for i in range(nx - 1):
        R = np.sqrt(r[i + 1] / r[i - 1])
        rh = R * r[i]
        uh = (R * u[i + 1] + u[i]) / (R + 1)
        Hh = (R * H[i + 1] + H[i]) / (R + 1)
        ah = np.sqrt((gamma - 1.0) * (Hh - 0.5 * uh ** 2))

        alpha1 = (gamma - 1) * uh ** 2 / (2 * ah ** 2)
        alpha2 = (gamma - 1) / ah ** 2

        d = q[:, i + 1] - q[:, i]

        L = np.array([
            [np.abs(uh - ah),               0,                                  0               ],
            [0,                             np.abs(uh),                         0               ],
            [0,                             0,                                  np.abs(uh + ah) ]
        ])
        P = np.array([
            [1,                             1,                                  1               ],
            [uh - ah,                       uh,                                 uh + ah         ],
            [Hh - ah * uh,                  0.5 * uh ** 2,                      Hh + ah * uh    ]
        ])
        Pinv = np.array([
            [0.5 * (alpha1 + uh / ah),      -0.5 * (alpha2 * uh + 1 / ah),       alpha2 / 2      ],
            [1 - alpha1,                    alpha2 * uh,                         -alpha2         ],
            [0.5 * (alpha1 - uh / ah),      -0.5 * (alpha2 * uh - 1 / ah),       alpha2 / 2      ]
        ])

        A = P @ L @ Pinv
        Ph[:, i] = A @ d

    Ph = 0.5 * (F[:, 0:nx -1] + F[:, 1:nx] - Ph)
    return (Ph[:, 1:] - Ph[:, :-1])

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
    dF = roeflux(q0, dx, gamma, a, nx)
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