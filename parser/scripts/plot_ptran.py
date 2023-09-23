import matplotlib.pyplot as plt
import numpy as np
from pylab import mpl
from sklearn.preprocessing import StandardScaler

mpl.rcParams["font.sans-serif"] = ["SimHei"]
mpl.rcParams["axes.unicode_minus"] = False


def readfile(path):
    timepoints = []
    voltages = []
    with open(path, 'r') as file:
        for line in file:
            fields = line.strip().split(',')
            t = float(fields[1].strip().split('=')[1])
            b = int(fields[2].strip().split('=')[1])
            x = fields[-1].strip().split('=')[1]
            if b == 0:
                timepoints.append(t)
                voltages.append([float(num) for num in x.split()])
    
    timepoints = np.array(timepoints).reshape(-1, 1)
    voltages = np.array(voltages)

    # tscaler = StandardScaler()
    # vscaler = StandardScaler()

    # timepoints = tscaler.fit_transform(timepoints)
    # voltages = vscaler.fit_transform(voltages)

    return timepoints.reshape(-1), voltages


if __name__ == '__main__':
    path = "../bin/ptran3.txt"
    t, v = readfile(path)
    for iter in range(len(v[0])):
        plt.figure()
        plt.scatter(t, v[:, iter], color='blue', marker='o')
        plt.xlabel('timepoints')
        plt.ylabel('voltages')
        plt.title('x' + str(iter + 1) + ' process')
        plt.show()
        plt.close()