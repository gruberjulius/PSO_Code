#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np


def plot_2D(N, R, name,xl,yl):
    plt.cla()
    plt.plot(N, R)
    #plt.plot(N, np.sqrt(N)) #auskommentiert für S
    plt.xlabel(xl)
    plt.ylabel(yl, rotation = 0)
    #plt.legend(["MD","sqrt(N)"]) #auskom für S
    plt.title(name)
    plt.savefig(name + '.pdf')


def main():
    r0 = 0.00154
    N, T, R, S = np.loadtxt('endtoend.out', delimiter = ' ', unpack = True)
    plot_2D(S, R, "end2end_R2S","Steps","R(S)")


if __name__ == '__main__':
    main()
