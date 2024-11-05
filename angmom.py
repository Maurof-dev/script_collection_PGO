import numpy as np
import matplotlib
from matplotlib import rc
import pylab as plt
from math import pi, sqrt










def set_size(width_pt, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)





    
y = set_size(510.0)
M = 1200
plt.rcParams['figure.figsize'] = set_size(510.0)
plt.rcParams['xtick.labelsize'] = set_size(M)[0]
plt.rcParams['ytick.labelsize'] = set_size(M)[0]
plt.rcParams['axes.labelsize'] = set_size(M)[0]





        





def phang(namefile):
    data = np.loadtxt(namefile)
    data = np.array([data[i][1:] for i in range(0,len(data))])
    qptarr = np.linspace(0,10,len(data[0]))
    Qptarr = np.linspace(0,10,len(data))
    band = []
    for q in range(0,len(qptarr)):
        band_q = np.array([data[i][q] for i in range(0,len(data))])
        band = band + [band_q]
    band = np.array(band)    
    return [Qptarr,band]





def histo_phang(namefile):
    data = phang(namefile)
    qpt = data[0]
    bands = data[1]
    L = len(bands)
    lvals = np.array(bands).reshape(L*len(qpt))
    return lvals


# Material: Pg5Ge3O11
# ABINIT calculations
# PE: space group P-6
# FE: space group P3
    

# GM to K direction
PE_nSOC_K = 'PHANGMOM_GM_to_K_PE_no_SOC'
FE_nSOC_K = 'PHANGMOM_GM_to_K_FE_no_SOC'
PE_wSOC_K = 'PHANGMOM_GM_to_K_PE_w_SOC'

# GM to K direction
PE_nSOC_A = 'PHANGMOM_GM_to_A_PE_no_SOC'
FE_nSOC_A = 'PHANGMOM_GM_to_A_FE_no_SOC'
PE_wSOC_A = 'PHANGMOM_GM_to_A_PE_w_SOC'


data = np.loadtxt(PE_nSOC_K)



check1 = histo_phang(PE_nSOC_K)
check2 = histo_phang(FE_nSOC_K)
check3 = histo_phang(PE_nSOC_A)
check4 = histo_phang(FE_nSOC_A)



num_bins = 60



fig, ax = plt.subplots()

# the histogram of the data
n, bins, patches = ax.hist(check1, num_bins, density=True,histtype='step',
                               color='blue',label=r'PE, $\Gamma$ $\rightarrow$ K')


n, bins, patches = ax.hist(check2, num_bins, density=True,histtype='step',
                               color='red',label=r'FE, $\Gamma$ $\rightarrow$ K')


n, bins, patches = ax.hist(check3, num_bins, density=True,histtype='step',
                               color='cyan',label=r'PE, $\Gamma$ $\rightarrow$ A')


n, bins, patches = ax.hist(check4, num_bins, density=True,histtype='step',
                               color='magenta',label=r'FE, $\Gamma$ $\rightarrow$ A')



# add a 'best fit' line
#y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
#ax.plot(bins, y, '--')
ax.set_xlabel('L$_{z}$ (Ha)')
ax.set_ylabel('counter')
ax.legend(fontsize='15')
#ax.set_title('Histogram of normal distribution sample')
ax.set_yscale('log')

# Tweak spacing to prevent clipping of ylabel
fig.tight_layout()
plt.show()









check = phang(PE_nSOC_K)
check1 = phang(FE_nSOC_K)




zeroline = np.zeros(len(data))
kptboundfile = np.array([0.0,10.0])


fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])


for q in range(0,len(check[1])):
    bndq = check[1][q]# * 0.124  # cm^-1 to meV
    qpts = check[0]
    bndq1 = check1[1][q]
    qpts1 = check1[0]
    ax.plot(qpts,bndq,'b-')
    ax.plot(qpts1,bndq1,'r-')
plt.plot(qpts,zeroline,'b-',label='PE')
plt.plot(qpts,zeroline,'r-',label='FE')



for i in range(0,len(kptboundfile)):
    x = kptboundfile[i]
    plt.axvline(x=kptboundfile[i])
ax.set_xticks([0.0,10.0])
ax.set_xticklabels(['$\Gamma$','K=(1/3,1/3,0)'])


plt.legend(fontsize="15")
plt.xlabel('wave vector')
plt.ylabel('L$_{z}$ (Ha)')
plt.show()


