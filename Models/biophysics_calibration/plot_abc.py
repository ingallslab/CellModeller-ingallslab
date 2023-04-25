import pyabc
import matplotlib.pyplot as plt
import numpy as np
import pandas

plt.rcParams['font.size'] = 24
colors = ['b','grey','r','g','chocolate','gold','m']

def ABCsimulation(params):
    return None

def plot_1d(history, limits, exact, key, **kwargs):
    """
    Plot posterior distributions of an ABC run
    @param history  ABC database
    @param limits   dict containing upper/lower bounds of priors
    @param exact    Exact value of the observed data
    @param key      Key in exact dict
    
    The figure is saved in the same directory as this script
    """

    xmin = limits[key][0]
    xmax = limits[key][1]
    
    fig = plt.figure(figsize=(9,11))
   
    #Plot prior distributions
    y = 1/(xmax-xmin) # uniform prior is equal to 1/max(param) for each parameter 
    plt.plot(np.linspace(xmin,xmax,100), y*np.ones(100,), color='k', lw=3, label='Prior')
    
    #Plot exact value
    plt.plot(exact[key], 0, ls=':', color='black', label='Observed Value')
    
    ax = plt.gca()
    c = 0
    for t in range(0,history.n_populations,1):
        df, w = history.get_distribution(t=t)
        pyabc.visualization.plot_kde_1d(
            df, w, 
            x=key, ax=ax,
            xmin=xmin,
            xmax=xmax,              
            numx=10000,
            label="Stage={}".format(t+1),
            refval=exact, refval_color="black", lw=3, color=colors[c])    
        plt.xlabel(key)
        ax.set_ylabel("") 
        c+=1
    ax.set_ylabel("")
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', ncol=1, prop={'size': 18})
    plt.savefig("%s_1D_KDE-t=%s.png" % (key,t), bbox_inches='tight')

if __name__ == '__main__':
    # Create "Null" ABCSMC object that has the corresponding parameter priors used in ABC run 
    lower_bound = 0
    scale = 1
    prior = pyabc.Distribution(mu=pyabc.RV("uniform", lower_bound, scale))
    abc = pyabc.ABCSMC(ABCsimulation, prior)
    
    # Load database
    db_path = ("sqlite:///" + "results.db")
    run_id = 1
    history = abc.load(db_path, run_id)
    
    ### Plots ###
    exact = {"gamma": 10, "reg_param": 0.1}
    limits = {"gamma": [10, 1000], "reg_param": [1/1000, 1/10]}
    for key in exact.keys():
        plot_1d(history, limits, exact, key)
    
