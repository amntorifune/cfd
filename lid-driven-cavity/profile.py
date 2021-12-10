import csv
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def draw128(Re, ls = '-'):
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
    
    x = np.array(x)
    y = np.array(y)
    x = x * 2 - 1
    y = y * 2 - 1

    plt.plot(x, v, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75)
    plt.plot(u, y, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (128x128)'%Re)

def draw256(Re, ls = '-'):
    x = []
    y = []
    u = []
    v = []
    fname = 'Re%d@fx.256.csv'%Re
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['x']) == 0.5 and 0 <= float(row['y']) <= 1:
                y.append(float(row['y']))
                u.append(float(row['uf']))
    
    fname = 'Re%d@fy.256.csv'%Re
    with open(fname) as f:
        csvf = csv.DictReader(f)
        for row in csvf:
            if float(row['y']) == 0.5 and 0 <= float(row['x']) <= 1:
                x.append(float(row['x']))
                v.append(float(row['vf']))

    x = np.array(x)
    y = np.array(y)
    x = x * 2 - 1
    y = y * 2 - 1

    plt.plot(x, v, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75)
    plt.plot(u, y, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (256x256)'%Re)

if __name__ == '__main__':
    draw128(100  )
    draw128(400  , ls = '--')
    draw128(1000 )
    draw256(3200 , ls = '--')
    draw256(5000 )
    draw256(7500 , ls = '--')
    draw256(10000)
    xt = np.linspace(- 1, 1, 11)
    yt = np.linspace(- 1, 1, 11)
    plt.xticks(xt)
    plt.yticks(yt)
    plt.grid(True)
    plt.legend()
    plt.show()