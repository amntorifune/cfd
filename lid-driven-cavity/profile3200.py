import csv
import matplotlib.pyplot as plt
import sys

def draw(Re, resol, linsol, ls):
    x = []
    y = []
    u = []
    v = []
    fname = 'Re%d@fx.%s%d.csv'%(Re, linsol, resol)
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['x']) == 0.5 and 0 <= float(row['y']) <= 1:
                y.append(float(row['y']))
                u.append(float(row['uf']))
    
    fname = 'Re%d@fy.%s%d.csv'%(Re, linsol, resol)
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['y']) == 0.5 and 0 <= float(row['x']) <= 1:
                x.append(float(row['x']))
                v.append(float(row['vf']))
    
    plt.subplot(1, 2, 1)
    plt.plot(x, v, ls, linewidth = 0.75, label = '3rd-order upwind (%dx%d)'%(resol, resol))
    plt.xlabel('x')
    plt.ylabel('v')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(u, y, ls, linewidth = 0.75, label = '3rd-order upwind (%dx%d)'%(resol, resol))
    plt.xlabel('u')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()

def drawGhia():
    x = [
        1.0000,
        0.9688,
        0.9609,
        0.9531,
        0.9453,
        0.9063,
        0.8594,
        0.8074,
        0.5000,
        0.2344,
        0.2266,
        0.1563,
        0.0938,
        0.0781,
        0.0703,
        0.0625,
        0.0000
    ]
    v = [
         0.00000,
        -0.39017,
        -0.47425,
        -0.52357,
        -0.54053,
        -0.44307,
        -0.37401,
        -0.31184,
         0.00999,
         0.28188,
         0.29030,
         0.37119,
         0.42768,
         0.41906,
         0.40917,
         0.39560,
         0.00000,
    ]
    y = [
        1.0000,
        0.9766,
        0.9688,
        0.9609,
        0.9531,
        0.8516,
        0.7344,
        0.6172,
        0.5000,
        0.4531,
        0.2813,
        0.1719,
        0.1016,
        0.0703,
        0.0625,
        0.0547,
        0.0000
    ]
    u = [
         1.00000,
         0.53236,
         0.48296,
         0.46547,
         0.46101,
         0.34682,
         0.19791,
         0.07156,
        -0.04272,
        -0.86636,
        -0.24427,
        -0.34323,
        -0.41933,
        -0.37827,
        -0.35344,
        -0.32407,
         0.00000
    ]

    plt.subplot(1, 2, 1)
    plt.scatter(x, v, c = 'g', marker = 'x', label = 'Ghia et al. (128x128)')
    plt.xlabel('x')
    plt.ylabel('v')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.scatter(u, y, c = 'g', marker = 'x', label = 'Ghia et al. (128x128)')
    plt.xlabel('u')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()

if __name__ == '__main__':
    Re = 3200
    fig, ax = plt.subplots(nrows = 1, ncols = 2)
    draw(Re, 128, sys.argv[1], 'b-' )
    draw(Re, 256, sys.argv[1], 'r--')
    drawGhia()
    plt.suptitle('Velocity at geometry center (Re = %d)'%Re)
    plt.show()
