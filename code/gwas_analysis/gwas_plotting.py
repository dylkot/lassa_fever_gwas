import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta


def qqplot(P, num_top_points=1000, numsample = 100000,
           fillcolor='b', dotcolor='k', fill_dens=0.1,
           title='title', ax=None,
           dot_size=5, rasterized=False):

    N = len(P)
    observed = -1*np.log10(np.sort(P))[::-1]
    expected = -1*np.log10(np.arange(N, 0, -1) / (N+1))
    
    if num_top_points > N:
        num_top_points = N
    
    ## Use only the last K theoretical quantiles to compute the confidence band
    bandK = np.arange(num_top_points, 0, -1)
    bandResid = N+1-bandK    
    lower = -1*np.log10(beta.ppf(0.975, bandK, bandResid))
    upper = -1*np.log10(beta.ppf(0.025, bandK, bandResid))
    
    if ax is None:
        (fig,ax) = plt.subplots(1,1, figsize=(5,5))

    #ax.fill_between(expected[(-1*num_top_points):], lower,upper, color=fillcolor, alpha=fill_dens)
    #ax.plot(expected[:(-1*numbandvals)], observed[:(-1*numbandvals)], color=dotcolor)
    #ax.scatter(expected[(-1*numbandvals):], observed[(-1*numbandvals):], s=dot_size, edgecolor='None', color=dotcolor, rasterized=rasterized)
        
    expected_upper = expected[(-1*num_top_points):]
    observed_upper = observed[(-1*num_top_points):]
    
    sample = sorted(np.random.choice(N-num_top_points, size=numsample))
    expected_lower = expected[sample]
    observed_lower = observed[sample]
    expected_plot = np.concatenate((expected_lower, expected_upper))
    observed_plot = np.concatenate((observed_lower, observed_upper))   
    
    ax.scatter(expected_plot, observed_plot, s=dot_size,
               edgecolor='None', color=dotcolor, rasterized=rasterized)
    
    maxval = np.max([upper[-1], observed[-1], expected[-1]])
    plotmax = maxval*1.02
    linemax = maxval*2
    
    ax.set_xlabel('Expected Quantile')
    ax.set_ylabel('Observed Quantile')
    ax.plot([0, linemax], [0, linemax],'--r')
    ax.set_xlim([0, plotmax])
    ax.set_ylim([0, plotmax])
    ax.set_title(title)
    return(ax)


def manhattan(BP, CHR, P=None, minusLogP=None, ax=None, colors=['k', 'r', 'b'],
             chrlist=None, chrlabels=None, plotheight=None, cut=0, spacer=500, dotsize=20,
             ylabel='-Log P-value', xlabel='Coordinate', title='Manhattan Plot',
             only_altticks=False, rasterized=False, linewidth=0.5):
    
    if ax is None:
        (fig,ax) = plt.subplots(1,1)
        
    if (minusLogP is None) and (P is None):
        sys.exit('Requires minusLogP or P')
    elif minusLogP is None:
        Y = -np.log10(P)
    elif P is None:
        Y = minusLogP
        
    if chrlist is None:
        chrlist = sorted(CHR.unique())

    if chrlabels is None:
        chrlabels = [str(x) for x in chrlist]     

    if plotheight is None:
        plotheight = np.max(Y)*1.1
    
    numchroms = len(chrlist)
    numcolors = len(colors)

    chrstart = [0]
    Ycut = Y>cut
    tickpositions = []
        
    for (i,c) in enumerate(chrlist):
        ind = (CHR==c) & Ycut
        color = colors[i%numcolors]
        X = np.arange(ind.sum())+chrstart[-1]
        ax.scatter(X, Y[ind], color=color, edgecolor='None', s=dotsize, alpha=.7, rasterized=rasterized)
        chrstart.append(X.max()+spacer)
        if only_altticks:
            if (i % 2) == 0:
                tickpositions.append(X[int(len(X)/2)])
        else:
            tickpositions.append(X[int(len(X)/2)])            
        
        if i != (numchroms-1):
            ax.plot([chrstart[-1] - spacer/2, chrstart[-1] - spacer/2], [cut, plotheight*1.5], '--', lw=linewidth, color='lightgray')

    if only_altticks:
        chrlabels = [x for x in chrlabels if x != '']

            
    ax.set_xticks(tickpositions)
    ax.set_xticklabels(chrlabels)
    ax.set_xlim([0-spacer, X[-1]+spacer])
    ax.set_ylim([cut, plotheight])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    
    
def manhattan_fixwidth(BP, CHR, P=None, minusLogP=None, ax=None, colors=['k', 'r', 'b'],
             chrlist=None, chrlabels=None, plotheight=None, cut=0, spacer=.05, dotsize=20,
             ylabel='-Log P-value', xlabel='Coordinate', title='Manhattan Plot',
             only_altticks=False, rasterized=False, linewidth=0.5):
    
    if ax is None:
        (fig,ax) = plt.subplots(1,1)
        
    if (minusLogP is None) and (P is None):
        sys.exit('Requires minusLogP or P')
    elif minusLogP is None:
        Y = -np.log10(P)
    elif P is None:
        Y = minusLogP
        
    if chrlist is None:
        chrlist = sorted(CHR.unique())

    if chrlabels is None:
        chrlabels = [str(x) for x in chrlist]     

    if plotheight is None:
        plotheight = np.max(Y)*1.1
    
    numchroms = len(chrlist)
    numcolors = len(colors)

    chrstart = [0]
    Ycut = Y>cut
    tickpositions = []
        
    for (i,c) in enumerate(chrlist):
        ind = (CHR==c) & Ycut
        color = colors[i%numcolors]
        X = np.linspace(0, 1, ind.sum())+chrstart[-1]
        ax.scatter(X, Y[ind], color=color, edgecolor='None', s=dotsize, alpha=.7, rasterized=rasterized)
        chrstart.append(X.max()+spacer)
        if only_altticks:
            if (i % 2) == 0:
                tickpositions.append(X[int(len(X)/2)])
        else:
            tickpositions.append(X[int(len(X)/2)])            
        
        if i != (numchroms-1):
            ax.plot([chrstart[-1] - spacer/2, chrstart[-1] - spacer/2], [cut, plotheight*1.5], '--',
                    lw=linewidth, color='lightgray')

    if only_altticks:
        chrlabels = [x for x in chrlabels if x != '']

            
    ax.set_xticks(tickpositions)
    ax.set_xticklabels(chrlabels)
    ax.set_xlim([0-spacer, X[-1]+spacer])
    ax.set_ylim([cut, plotheight])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)