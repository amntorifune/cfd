import csv
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys

def draw(Re, resol, linsol, ls = '-'):
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
    
    x = np.array(x)
    y = np.array(y)
    x = x * 2 - 1
    y = y * 2 - 1

    plt.plot(x, v, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75)
    plt.plot(u, y, c = cm.viridis(Re / 10000), ls = ls, linewidth = 0.75, label = '%d (128x128)'%Re)

if __name__ == '__main__':
    linsol = sys.argv[1]
    draw(100  , 128, linsol)
    draw(400  , 128, linsol, ls = '--')
    draw(1000 , 128, linsol)
    draw(3200 , 256, linsol, ls = '--')
    draw(5000 , 256, linsol)
    draw(7500 , 256, linsol, ls = '--')
    draw(10000, 256, linsol)
    # draw256(3200 , ls = '--')
    # draw256(5000 )
    # draw256(7500 , ls = '--')
    # draw256(10000)
    xt = np.linspace(- 1, 1, 11)
    yt = np.linspace(- 1, 1, 11)
    plt.xticks(xt)
    plt.yticks(yt)
    plt.grid(True)
    plt.legend()
    plt.show()