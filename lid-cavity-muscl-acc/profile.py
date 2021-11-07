import csv
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def draw128(Re, ls = '-'):
    fname = 'UVP_Re%d.ns2dccg.csv'%Re
    x = []
    y = []
    u = []
    v = []
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['x']) == 0.5:
                y.append(float(row['y']))
                u.append(float(row['u']))
            if float(row['y']) == 0.5:
                x.append(float(row['x']))
                v.append(float(row['v']))
    
    x = np.array(x)
    y = np.array(y)
    x = x * 2 - 1
    y = y * 2 - 1

    plt.plot(x, v, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (128x128)'%Re)
    plt.plot(u, y, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (128x128)'%Re)

def draw256(Re, c = 'b', ls = '-'):
    fname = 'UVP_Re%d.ns2dccg256.csv'%Re
    x = []
    y = []
    u = []
    v = []
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['x']) == 0.5:
                y.append(float(row['y']))
                u.append(float(row['u']))
            if float(row['y']) == 0.5:
                x.append(float(row['x']))
                v.append(float(row['v']))

    x = np.array(x)
    y = np.array(y)
    x = x * 2 - 1
    y = y * 2 - 1

    plt.plot(x, v, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (256x256)'%Re)
    plt.plot(u, y, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (256x256)'%Re)

if __name__ == '__main__':
    draw128(100)
    draw128(400)
    draw128(1000)
    draw256(3200)
    draw256(5000)
    draw256(7500)
    draw256(10000)
    xt = np.linspace(- 1, 1, 11)
    yt = np.linspace(- 1, 1, 11)
    plt.xticks(xt)
    plt.yticks(yt)
    plt.grid(True)
    plt.legend()
    plt.show()