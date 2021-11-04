import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np


def plot_2D(N, R, name,xl,yl):
    plt.cla()
    plt.plot(N, R)
    #plt.plot(N, np.sqrt(N)) #auskommentiert für S
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(xl)
    plt.ylabel(yl, rotation = 0)
    #plt.legend(["MD","sqrt(N)"]) #auskom für S
    plt.title("log log plot")
    plt.savefig(name + '.pdf')
'''
def plot_bar(values,name):
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    params = ['C', 'C++', 'Java', 'Python', 'PHP','C', 'C++', 'Java', 'Python', 'PHP']
    ax.bar(params,values)
    plt.savefig(name + '.pdf')
'''
def main():
    G, X, R = np.loadtxt('bestfit.dat', delimiter = ' ', unpack = True)
    plot_2D(G, -R, "bestfit1_3","Generation","-(Best Value)")
    plot_2D(G,X,"2ndcol_to_G","Generation","second columns")
#after bash ./../testing --load gen_200000 -o outpufile200000 -n 100
#    values = np.loadtxt('outputfile200000', delimiter = ' ', unpack = True)
#    plot_bar(values,"best_ind_parameters")


if __name__ == '__main__':
    main()
