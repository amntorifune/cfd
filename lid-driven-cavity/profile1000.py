import csv
import matplotlib.pyplot as plt

def draw128(Re):
    x = []
    y = []
    u = []
    v = []
    fname = 'Re%d@fx.128.csv'%Re
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['x']) == 0.5 and 0 <= float(row['y']) <= 1:
                y.append(float(row['y']))
                u.append(float(row['uf']))
    
    fname = 'Re%d@fy.128.csv'%Re
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['y']) == 0.5 and 0 <= float(row['x']) <= 1:
                x.append(float(row['x']))
                v.append(float(row['vf']))
    
    plt.subplot(1, 2, 1)
    plt.plot(x, v, 'b-', linewidth = 0.75, label = '3rd-order upwind (128x128)')
    plt.xlabel('x')
    plt.ylabel('v')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(u, y, 'b-', linewidth = 0.75, label = '3rd-order upwind (128x128)')
    plt.xlabel('u')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()

def draw256(Re):
    x = []
    y = []
    u = []
    v = []
    fname = 'Re%d@f.256.csv'%Re
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['x']) == 0.5 and 0 <= float(row['y']) <= 1:
                y.append(float(row['y']))
                u.append(float(row['uf']))
    
    fname = 'Re%d@f.256.csv'%Re
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['y']) == 0.5 and 0 <= float(row['x']) <= 1:
                x.append(float(row['x']))
                v.append(float(row['vf']))

    plt.subplot(1, 2, 1)
    plt.plot(x, v, 'r--', linewidth = 0.75, label = '3rd-order upwind (256x256)')
    plt.xlabel('x')
    plt.ylabel('v')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(u, y, 'r--', linewidth = 0.75, label = '3rd-order upwind (256x256)')
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
        -0.21388,
        -0.27669,
        -0.33714,
        -0.39188,
        -0.51550,
        -0.42665,
        -0.31966,
         0.02526,
         0.32235,
         0.33075,
         0.37095,
         0.32627,
         0.30353,
         0.29012,
         0.27485,
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
         0.65928,
         0.57492,
         0.51117,
         0.46604,
         0.33304,
         0.18719,
         0.05702,
        -0.06080,
        -0.10648,
        -0.27805,
        -0.38289,
        -0.29730,
        -0.22220,
        -0.20196,
        -0.18109,
         0.00000
    ]

    plt.subplot(1, 2, 1)
    plt.scatter(x, v, c = 'g', marker = 'x', label = 'Ghia et al. (256x256)')
    plt.xlabel('x')
    plt.ylabel('v')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.scatter(u, y, c = 'g', marker = 'x', label = 'Ghia et al. (256x256)')
    plt.xlabel('u')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()

if __name__ == '__main__':
    fig, ax = plt.subplots(nrows = 1, ncols = 2)
    draw128(1000)
    # draw256(3200)
    drawGhia()
    plt.suptitle('Velocity at geometry center (Re = 1000)')
    plt.show()