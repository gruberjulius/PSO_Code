import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import csv
import string


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
    I, T = np.loadtxt('times_ps.txt', delimiter = ';', unpack = True)
    plot_2D(I, T, "time_ps","Iterations","Time")


if __name__ == '__main__':
    main()


"""    x = []
    y = []
    with open('results/times_ps.txt','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
           y.append(float(row[1])) 
"""
""""    plt.plot(x,y, label='Loaded from file!')
    plt.xlabel('Iteration')
    plt.ylabel('Time')
    plt.title('Interesting Graph\nCheck it out')
    plt.legend()
    plt.show()"""