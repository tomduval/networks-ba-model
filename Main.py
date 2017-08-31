import BAModelCPP as bam
import BAAnalysis as baa
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp

def VaryingN():
    """produces data for the BAModel using the CPP module with increasing N
    Uses: the model BAModelCPP"""
    data = []
    maxdeg = []
    for i in range(2, 5):
        print i
        data.append(bam.BAModel(N=10**i, n0=4, m=3, L=1, initiate=2, repeats=1000, graph = "BA"))
        data[-1].Iterate()
        maxdeg.append(data[-1].maxdeg)
    return data, maxdeg
    
def VaryingNMany():
    """produces data for the BAModel using the CPP module with increasing N
    Uses: the model BAModelCPP"""
    data = []
    i_range = np.linspace(2,5,10)
    maxdeg = []
    for i in i_range:
        print i
        data.append(bam.BAModel(N=int(10**i), n0=4, m=3, L=1, initiate=2, repeats=200, graph = "BA"))
        data[-1].Iterate()
        maxdeg.append(data[-1].maxdeg)
    return data, maxdeg

def VaryingM():
    """produces data for the BAModel using the CPP module with increasing m
    Uses: the model BAModelCPP"""
    data = []
    maxdeg = []
    for i in range(1,5):
        print i
        data.append(bam.BAModel(N=10**4, n0=3**i+1, m=3**i, L=1, initiate=2, repeats=1000, graph = "BA"))
        data[-1].Iterate()
        maxdeg.append(data[-1].maxdeg)
    return data, maxdeg
        
def Graphs(data, maxdeg):
    """produces graphs
    Uses the BAAnalysis module"""
    analysis = []
    for i in range(len(data)):
        analysis.append(baa.BAModelAnalysis(degrees=data[i].degrees, n0=data[i].n0, m=data[i].m, N=data[i].N, graphtype=data[i].graphtype))
        analysis[i].LogBin(a=1.2)
        analysis[i].DegreeProbabilityRawVsLBPlot()
        analysis[i].GradientDegreeProbability()
        analysis[i].Statistics()
        analysis[i].DegreeProbabilityCollapsedPlot(maxdeg[i])
        analysis[i].PValue()
        analysis[i].DegreeProbabilityPlot()
        analysis[i].DegreeProbabilityTheoreticalPlot()
    
def MaxDegreesPlot(data, maxdeg):
    """plots the maximum degree against the theoretical maximum degree"""
    maxdeg_mean = np.mean(np.array(maxdeg), axis=1)
    maxdeg_std = np.std(np.array(maxdeg), axis=1)
    N = [data[i].N for i in range(len(data))]
    
    if data[0].graphtype == "Pure Random":
        theoretical = [data[i].m-np.log(data[i].N)/(np.log(data[i].m)-np.log(data[i].m+1)) for i in range(len(data))]
    else:
        theoretical = [(-1+np.sqrt(1+4*data[i].N*data[i].m*(data[i].m+1)))/2 for i in range(len(data))]
        
    maxdeg_grad, maxdeg_intercept = sp.linregress(np.log10(N), np.log10(maxdeg_mean))[:2]
    plt.figure('MaximumDegreeAgainstN')
    plt.errorbar(x=N, y=maxdeg_mean, yerr=maxdeg_std, label="Numerical")
    plt.plot(N, theoretical, '--', label="Theoretical")
    plt.title('Maximum degree against the total number of nodes')
    plt.xlabel(r'$N$')
    plt.ylabel(r'$k_1$')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc='best')
    print "Maximum degree log-log gradient: ", maxdeg_grad
    
#G = nx.barabasi_albert_graph(n=1000, m=4, seed=500)
#analysis = BAModelAnalysis(model.edges)
#analysis.Plot()
#analysis.ColourPlot(ncols=1000)
#analysis.DegreeFrequencyPlot()