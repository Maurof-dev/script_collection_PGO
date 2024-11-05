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








# 62 qpoints

# natom = 57 ===> 57x3 = 171 bands x qpoint

# qpt[0] = 0 0 0, qpt[61] = 0 0 1/2

# frequencies in meV





def phbands(namefile,nacfile,flag):
    data = np.loadtxt(namefile)
    data = np.array([data[i][1:] for i in range(0,len(data))])
    qptarr = np.linspace(0,10,len(data[0]))
    Qptarr = np.linspace(0,10,len(data))
    band = []
    for q in range(0,len(qptarr)):
        #band_q = np.array([data[i][q] for i in range(0,len(data))])
        band_q = np.array([data[i][q] for i in range(0,len(data))])
        band = band + [band_q]
    band = np.array(band)*8.06554# meV to cm-1 conversion
    if flag == 0.0:
        band = band
    else:
        nacdata = open(nacfile)
        nacdata = nacdata.readlines()
        S = len(nacdata)-2
        naclist = []
        for i in range(0,S):
            nac = nacdata[i].split()
            nac = np.array([float(nac[j]) for j in range(0,len(nac))])
            naclist = naclist + [nac]
        naclist = np.array(naclist).reshape(S*5)
        naclist = [naclist[i] for i in range(0,len(naclist))] + [3.446761E-03]
        #naclist = np.array(naclist)*27211.3825435 # Ha to meV conversion
        naclist = np.array(naclist)*220000.00000000003  # Ha to cm-1 conversion
        bandnew = []
        for n in range(0,len(band)):
            bandnew = bandnew + [nac_correction(band[n],naclist,n)] # along z, in Hartree ---> convert in meV
        band = np.array(bandnew)
    return [Qptarr,band]



def nac_correction(bandq,naclist,n):
    bandqnew = []
    for j in range(0,len(bandq)):
        if j == 0.0:
            bandqnew = bandqnew + [naclist[n]]
        else:
            bandqnew = bandqnew + [bandq[j]]
    return np.array(bandqnew)
            




freqs = 'PE_freqs'# GM to A
nacfile = 'PE_freq_Hartree_NAC_Gamma.txt'

data = np.loadtxt(freqs)


check = phbands(freqs,nacfile,3.0)

zeroline = np.zeros(len(data))
kptboundfile = np.array([0.0,10.0])


fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])

for q in range(0,len(data[0])-1):
    bndq = check[1][q]# * 0.124  # cm^-1 to meV
    qpts = check[0]
    ax.plot(qpts,bndq,'g-')
plt.plot(qpts,zeroline,'k-')#,label='PGO, LR interp.')



for i in range(0,len(kptboundfile)):
    x = kptboundfile[i]
    plt.axvline(x=kptboundfile[i])
ax.set_xticks([0.0,10.0])
ax.set_xticklabels(['$\Gamma$','A=(0,0,1/2)'])


#plt.legend(fontsize="15")
plt.xlabel('qpt index')
plt.ylabel('$\omega$ (cm$^{-1}$)')
#plt.ylabel('$\omega$ (meV)')
plt.show()

