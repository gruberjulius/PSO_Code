import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import csv
import string
import statistics


def plot_2D(x_1, y_1, x_2, y_2, name,xl,yl):
    plt.cla()
    plt.plot(x_1, y_1, 'C1', label='s')
    plt.plot(x_2, y_2, 'C2', label='m')
    #plt.plot(x_1, np.sqrt(N)) #auskommentiert für S
    plt.xlabel(xl)
    plt.ylabel(yl, rotation = 0)
    #plt.legend(["MD","sqrt(N)"]) #auskom für S
    plt.title(name)
    plt.savefig(name + '.pdf')

def mean_median_min_max(path, i1_s, t1_s, i1_m, t1_m, name1, i2_s, t2_s, i2_m, t2_m, name2):

    w_file = path + name1 + "_" + name2 + "_stat.txt"
    f = open(w_file, "a")
    #f.write("Now the file has more content!")

    print(name1 , "    vs.    ", name2)
    mean1_s = statistics.mean(t1_s)
    mean1_m = statistics.mean(t1_m)
    mean2_s = statistics.mean(t2_s)
    mean2_m = statistics.mean(t2_m)

    print("Mean s: %s, mean m: %s | prec_v s: %s, prec_v m: %s" % (mean1_s, mean1_m, mean2_s, mean2_m)) 
    f.write("Mean s: %s, mean m: %s | prec_v s: %s, prec_v m: %s \n" % (mean1_s, mean1_m, mean2_s, mean2_m)) 

    med1_s = statistics.median(t1_s)
    med1_m = statistics.median(t1_m)
    med2_s = statistics.median(t2_s)
    med2_m = statistics.median(t2_m)

    print("Median s: %s, median m: %s | prec_v s: %s, prec_v m: %s" % (med1_s, med1_m, med2_s, med2_m)) 
    f.write("Median s: %s, median m: %s | prec_v s: %s, prec_v m: %s \n" % (med1_s, med1_m, med2_s, med2_m)) 

    min1_s = min(t1_s)
    min1_m = min(t1_m)
    min2_s = min(t2_s)
    min2_m = min(t2_m)

    print("Min s: %s, min m: %s | prec_v s: %s, prec_v m: %s" % (min1_s, min1_m, min2_s, min2_m)) 
    f.write("Min s: %s, min m: %s | prec_v s: %s, prec_v m: %s \n" % (min1_s, min1_m, min2_s, min2_m)) 

    max1_s = max(t1_s)
    max1_m = max(t1_m)
    max2_s = max(t2_s)
    max2_m = max(t2_m)

    print("Max s: %s, max m: %s | prec_v s: %s, prec_v m: %s" % ( max1_s, max1_m, max2_s, max2_m))
    f.write("Max s: %s, max m: %s | prec_v s: %s, prec_v m: %s \n" % ( max1_s, max1_m, max2_s, max2_m))
    
    f.close()

def all(path, name1_, name2_):
    name1_s = path + name1_ + "_s.dat"
    name1_m = path + name1_ + "_m.dat"

    i1_s, t1_s = np.loadtxt(name1_s, delimiter = ';', unpack = True)
    i1_m, t1_m = np.loadtxt(name1_m, delimiter = ';', unpack = True)

    name2_s = path + name2_ + "_s.dat"
    name2_m = path + name2_ + "_m.dat"

    i2_s, t2_s = np.loadtxt(name2_s, delimiter = ';', unpack = True)
    i2_m, t2_m = np.loadtxt(name2_m, delimiter = ';', unpack = True)
    
    
    plot_2D(i1_s, t1_s, i1_m, t1_m, name1_,"Iterations","Time")    
    plot_2D(i2_s, t2_s, i2_m, t2_m, name2_,"Iterations","Time")
    mean_median_min_max(path, i1_s, t1_s, i1_m, t1_m, name1_, i2_s, t2_s, i2_m, t2_m, name2_)

def main():
    path = "results_benchmark/"
    it = 10000
    thr = 4
    name_1 = "benchmark_times_paps_" + str(it) + "_" + str(thr) #"benchmark_times_paps_10000_4"
    name_2 = "benchmark_times_paps_prec_vel_" + str(it) + "_" + str(thr) #"benchmark_times_paps_prec_vel_10000_4"

    all(path, name_1, name_2)

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