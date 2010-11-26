from matplotlib import rc
import matplotlib.pyplot as plt


# Tip: Use \showthe\columnwidth to get column width
def prepareCanvas(fig_width_pt = 350.0, s1 = 10, s2 = 8, lw = 0.5, plotmargins = [0.125, 0.15, 0.05, 0.05]): 
    #
    # Configure matplotlib environment
    # http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    #
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]

    rc('figure', figsize=fig_size)
    rc('lines', linewidth=lw)
    rc('font', family='serif', serif=['Latin Modern Roman','Palatino'], size=s1)
    rc('text', usetex=True)
    rc('legend', fontsize=s1)
    rc('axes', labelsize=s1)
    rc('xtick', labelsize=s2)
    rc('ytick', labelsize=s2)
    p = plotmargins
    plt.axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
