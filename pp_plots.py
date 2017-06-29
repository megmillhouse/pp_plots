import numpy as np
import matplotlib.pyplot as plt
import math


def pp_cdf(megachain, tru, Nbins=11):
    """
    For each run, finds the smallest credible interval that contains the true value of the parameter, and returns a cumulative sum of those
    
    Parameters
    ----------
    megachain : array_like[Nchains, Nsample]
        The all the chains of the parameter you're looking at
                
    tru : array[Nchains]
        Array of the injected values of your parameter (drawn from your prior). Must be in same order as megachain
        
    Nbins : int
        Number of bins
        
    Returns
    ------
    bins : array_like[Nbins]
        Center of the bins (i.e. x-axis of plot)
        
    CDF : array_like[Nbins]
        The normalized cumulative number of events in the ith credible interval
        
    """
    
    Nruns = len(megachain) # found out how many runs we did
    
    bins = np.linspace(0.0,1.0,num=Nbins)
    
    ci_count = np.zeros(Nbins)
    ci = np.zeros(2)
    
    for i in range(0,Nruns):
        for j in range(0,Nbins):
            ci[0] = np.percentile(megachain[i],(1.-bins[j])/2.*100.) # Get upper and lower bound on CIs
            ci[1] = np.percentile(megachain[i],(1.+bins[j])/2.*100.)
            if ci[0] <= tru[i] <= ci[1]: # find lowest CI that contains true value
                ci_count[j] += 1
                break

    # return bins and cumumulative summed CI counts
    return (bins,np.cumsum(ci_count)/Nruns)

def CDF_stds(Nb, Nsample, sigma):
    """
        Find the expected sigma-th standard deviation at each credible interval. Used for illustrative plotting.
        
        Parameters
        ----------
        Nbins : int
            Number of bins
        
        Nsample : int
            Number of samples in your distribution
        
        sigma : float
            Standard deviation
        
        Returns
        ------
        err : array_like[Nbins]
            The sigma-th standard deviation
        
    """
    n = float(Nsample)/float(Nb)
    NN = float(Nsample)
    err = []
    for ii in range(0,int(Nb)+1):
        i = float(ii)
        err.append(i*n/NN+sigma*math.sqrt(i*n*(1.-i*n/NN))/NN)
    return(err)


def make_plots(bins, pp_cdf, Nruns, param_name, sigmas=True, directdraw=False, diag=False):
    """
        Make your p-p plots! With options for how to present them
        
        Parameters
        ----------
        pp_cdf : array_like
            CDF values of your distribution
        
        Nbins : int
            Number of bins
        
        Nruns : int
            Number of chains
        
        param_name : string
            Name of particular parameter in this plot
            
        sigmas : bool
            Show or hide standard deviation
            
        directdraw : bool
            Show or hide CDFs of direct draws from U[0,1] for comparison
            
        diag : bool
            Show or hide diagonal like of slope 1
        
        Returns
        ------
        fig : figure
            Your p-p plot!
        
    """
    
    Nbins = len(bins)
    
    fig = plt.figure()
    
    plt.grid(linestyle='-',alpha=0.1)
    plt.rcParams['axes.axisbelow'] = True
    
    if sigmas:
        for i in range(1,4):
            err1 = CDF_stds(Nbins, Nruns, float(i))
            err2 = CDF_stds(Nbins, Nruns, -float(i))
            x = np.linspace(0,1.,int(Nbins)+1)
            plt.fill_betweenx(x,err1,err2,alpha=1.-(i-1)*.3,color='slategray')


    if directdraw:
         for i in range(0,100):
             h, e = np.histogram(np.random.uniform(0,1,Nruns),bins=Nbins)
             h = np.cumsum(h)/float(Nruns)
             c = (e[:-1] + e[1:]) / 2
             plt.plot(c,h,color='lightblue',alpha=0.5)

    if diag:
        plt.plot([0,1],[0,1],'r--',alpha=0.5)

    plt.plot(bins,pp_cdf,color='indigo')


    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xticks(np.arange(0,1.1,.1))
    plt.yticks(np.arange(0,1.1,.1))
    plt.xlabel('p',fontsize=16)
    plt.ylabel('P(p)',fontsize=16)
    plt.title(param_name,fontsize=24)

    return fig



