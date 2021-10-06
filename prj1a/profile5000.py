import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

def draw(Re, ls = 'b-'):
    x = []
    y = []
    u = []
    v = []
    with open("UVP_Re%d_t100.ns2dmuscl1.csv"%Re) as f:
        l = f.readline()
        for l in f:
            para = np.array(l.split(","), dtype=float)
            if para[0] == 0.5:
                y.append(para[1])
                u.append(para[3])
            if para[1] == 0.5:
                x.append(para[0])
                v.append(para[4])
        f.close()
    print(len(x), len(y), len(u), len(v))
    plt.subplot(1, 2, 1)
    plt.plot(x, v, ls, linewidth = 0.75, label = 'Upwind(256x256)')
    
    plt.subplot(1, 2, 2)
    plt.plot(u, y, ls, linewidth = 0.75, label = 'Upwind(256x256)')

def drawcomp():
    y = np.array([
        0.0,
        0.0547,
        0.0625,
        0.0703,
        0.1016,
        0.1719,
        0.2813,
        0.4531,
        0.5000,
        0.6172,
        0.7344,
        0.8516,
        0.9531,
        0.9609,
        0.9688,
        0.9766,
        1.0
    ])
    u = np.array([
        0.0,
        -0.41165,
        -0.42901,
        -0.43643,
        -0.40435,
        -0.33050,
        -0.22855,
        -0.07404,
        -0.03039,
        0.08183,
        0.20087,
        0.33556,
        0.46036,
        0.45992,
        0.46120,
        0.48223,
        1.0
    ])
    plt.subplot(1, 2, 2)
    plt.scatter(u, y, marker = 'x', color = 'g', label = 'Ghia et al.(256x256)')
    x = np.array([
        0.0,
        0.0625,
        0.0703,
        0.0781,
        0.0938,
        0.1563,
        0.2266,
        0.2344,
        0.5000,
        0.8047,
        0.8594,
        0.9063,
        0.9453,
        0.9531,
        0.9609,
        0.9688,
        1.0
    ])
    v = np.array([
        0.0,
        0.42447,
        0.43329,
        0.43648,
        0.42951,
        0.35368,
        0.28066,
        0.27280,
        0.00945,
        -0.30018,
        -0.36214,
        -0.41442,
        -0.52876,
        -0.55408,
        -0.55069,
        -0.49774,
        0.0
    ])
    plt.subplot(1, 2, 1)
    plt.scatter(x, v, marker = 'x', color = 'g', label = 'Ghia et al.(256x256)')

if __name__ == "__main__":
    fig, ax = plt.subplots(nrows = 1, ncols = 2)
    plt.suptitle('Re = 5000')
    drawcomp()
    draw(5000)

    plt.subplot(1, 2, 1)
    plt.xlabel('x')
    plt.ylabel('v')
    plt.grid(True)
    plt.legend(loc = 'lower left')

    plt.subplot(1, 2, 2)
    plt.xlabel('u')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend(loc = 'lower right')

    plt.show()