########################################################
# Plotting: (for the moment)
########################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
from math import log, log10
from scipy.integrate import quad
from pprint import pprint
from scipy import interpolate
import sys

def plotObjects():
    """Plot all the objects coming from ..?"""
    # Open pickle files:
    mchis = pickle.load(open('saveE.p'))
    ObjList = pickle.load(open('saveObj.p'))

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel(r'<$\sigma$v>$_{\textrm{UL}}$ (cm$^{3}$s$^{-1}$)')
    ax1.set_xlabel(r'$m_{\chi}$ (GeV)')
    ax1.set_ylim(1e-26,1e-20)
    ax1.set_xlim(1e2,1e5)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')
    
    #Plot stuff:
    for o in ObjList:
        ax1.plot(mchis,o.UL,label=o.Name)
        ## if o.Name=="SculptorIso":
        ##     ax1.plot(mchis,o.ul_WW,label=o.Name)
        ##     #print o.Aeffdata
        ## elif o.Name=="Segue1V_tautau":
        ##     ax1.plot(mchis,o.ul_tautau,label=r"Segue1Vtautau")
        ## else:
        ##     ax1.plot(mchis,o.ul_bbbar,label=o.Name)
        ##     #if o.Name=="Willman1V":
        ##         #print o.ul_bbbar
        ##         #print o.Aeffdata

    ## for o in PubList:
    ##     ax1.plot(o[:,0],o[:,1])

    #Plot Pubs?
    Pubs = True
    Pubs = False

    # Plot combination?
    Comb=True

    if Pubs:
        plotPubs()
    else:
        if Comb:
            plotComb()
        else:
            plt.show()
            plt.legend(loc=4,ncol=2,prop=matplotlib.font_manager.FontProperties(size='small'))
            plt.show()
    
    return 0

def plotPubs():
    """Plot the published limits for comparison"""

    # Open pickle files:
    PubList = pickle.load(open('savePub.p'))

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel(r'<$\sigma$v>$_{\textrm{UL}}$ (cm$^{3}$s$^{-1}$)')
    ax1.set_xlabel(r'$m_{\chi}$ (GeV)')
    ax1.set_ylim(1e-26,1e-20)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')

    # Reset color map?
    #matplotlib.cm.get_cmap()
    #matplotlib.colors.Normalize()
    ax1.set_color_cycle(('b','g','r','c','m','y'))
    
    # Plot stuff:
    for p in PubList:
        ax1.plot(p.mchis, p.ul, label=p.legend, linestyle='dashed')

    plt.legend(loc=4,ncol=2,prop=matplotlib.font_manager.FontProperties(size='small'))

    Comb=True
    if Comb:
        plotComb()
    else:
        plt.show()

def plotComb():
    """Plot the combined limit(s) for comparison"""

    # Open pickle files:
    mchis = pickle.load(open('saveE.p'))
    CombList = pickle.load(open('saveComb.p'))
    ## print
    ## print mchis
    ## print CombList
    ## print

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel(r'<$\sigma$v>$_{\textrm{UL}}$ (cm$^{3}$s$^{-1}$)')
    ax1.set_xlabel(r'$m_\chi$ (GeV)')
    ax1.set_ylim(1e-27,1e-20)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')

    # Reset color map?
    #matplotlib.cm.get_cmap()
    #matplotlib.colors.Normalize()
    ax1.set_color_cycle(('b','g','r','c','m','y'))

    ax1.plot(mchis, CombList[0], color='black', lw=3, 
             linestyle='dotted', label='Simply combined limit')
    ax1.plot(mchis, CombList[1], color='magenta', lw=4, 
             linestyle='dashed', label='Combined l\'hood limit, with Jbar PDFs')
    
    plt.legend(loc=4,ncol=2,prop=matplotlib.font_manager.FontProperties(size='small'))
    plt.show()

    
def plotAeffs():
    # Open pickle file:
    ObjList = pickle.load(open('saveObj.p'))
    mchis = pickle.load(open('saveE.p'))
    
    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim(1e1,1e5)
    ax1.set_ylim(1e7,1e10)
    ax1.set_xlabel('E (GeV)')
    ax1.set_ylabel(r'A$_{\rm eff}$ (cm$^2$)')
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')
    
    #Plot stuff:
    for o in ObjList:
        print o.Name
        ## # Calculate mean effective area for each instrument: 
        ## print 'Aeff = %e' % ((quad(lambda E: o.Aeff(E),
        ##                            o.Aeffdata[0,0], o.Aeffdata[-1,0])[0])
        ##                      /(o.Aeffdata[-1,0]- o.Aeffdata[0,0]))
        # Plot eff. area:
        ax1.plot(o.Aeffdata[:,0],o.Aeffdata[:,1],label=o.Name)
    
    plt.legend(["MAGIC","VERITAS","HESS Sgr","HESS Scu"],loc=4)
    plt.legend(loc=4)
    plt.show()


def plotSpectra():
    # Import spectra:
    from PhotonSpectra import tautau, bbbar, Bergstrom1998
    # Data:
    mchi=1000.
    mchis = np.linspace(30.,1010.,1000)
    tautaulist = []
    bbbarlist = []
    Berglist = []
    for e in mchis:
        tautaulist.append(tautau(e,mchi))
        bbbarlist.append(bbbar(e,mchi))
        Berglist.append(Bergstrom1998(e,mchi))

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlabel('E (GeV)')
    ax1.set_ylabel(r'dN/dE (GeV$^{-1}$)')
    ax1.set_xlim(30.,1010.)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')
    ax1.set_title(r'Spectra for $m_\chi$ = 1000 GeV')

    plot1 = ax1.plot(mchis, tautaulist, label=r'$\tau\tau$')
    plot2 = ax1.plot(mchis, bbbarlist, label=r'$b\bar{b}$')
    plot3 = ax1.plot(mchis, Berglist, label=r'$WW$')
    ## plot1 = ax1.plot(mchis, mchis**2*tautaulist, label=r'$\tau\tau$')
    ## plot2 = ax1.plot(mchis, mchis**2*bbbarlist, label=r'$b\bar{b}$')
    ## plot3 = ax1.plot(mchis, mchis**2*Berglist, label=r'$WW$')

    plt.legend(loc='lower left')
    plt.show()


def plotCLmins():
    """Plot some combined profile likelihoods
    ... and/or test interpolation etc."""
    # Open pickle files:
    mchis = pickle.load(open('saveE.p'))
    sigmavTestRange = np.arange(-26, -18, (26.-18.)/len(mchis)) #cf. Exclusion.py!
    # Array of minimized CLs:  CLmins[mchi][sv]
    CLmins = pickle.load(open('saveCLmins.p'))

    # Finer sv range:
    svtr2 = np.arange(-25, -19, (26.-18.)/(10.*len(mchis)))

    print 'mchis: ', mchis
    # Print CLmins array:
    #pprint(CLmins)

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
#    ax1.set_yscale('log')
#   ax1.set_xscale('log')
    ax1.set_ylabel(r'$\mathcal{CL}$')
    ax1.set_xlabel(r'$\log$<$\sigma$v>')
    ax1.set_ylim(260.,300.)
    ax1.set_xlim(-25.,-18)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')


    # Calculate some values:
    #ArrayMean = np.mean(CLmins)
    #ArrayStd = np.std(CLmins) # Complete array: Std very large!!!
    ArrayMean = np.mean(CLmins[:][0])
    ArrayStd = np.std(CLmins[:][0]) # first column: Std very small!!
    print 'ArrayMean, ArrayStd: ', ArrayMean, ArrayStd
    ArrayMin = np.min(CLmins)
    print 'ArrayMin =', ArrayMin

    #sys.exit()

    # Plot data from CLmins array:
    for i in range(len(CLmins[:][0])): # length = number of rows
        #ax1.plot(sigmavTestRange, CLmins[i][:])#, label=r'$\tau\tau$')

        # Refer(!) to a specific row:
        CLmZeile = CLmins[i][:] 

        # Create a clipped 1D array:
        # CLmSmooth = np.clip(CLmZeile, 0., ArrayMean+1.*ArrayStd) # not a good idea
        #CLmSmooth = np.clip(CLmZeile, 0., 300.)
        CLmSmooth = np.clip(CLmZeile, 0., ArrayMin+3.)

        #print CLmSmooth


        # Interpolation/Smoothing for each row's elements:
        for j in range(len(CLmSmooth)):
                    ## Mean = np.mean(CLmZeile[j-2:j+2])
            ## Std = np.std(CLmZeile[j-2:j+2])
            ## Mean = np.mean(CLmZeile[j-5:j+5])
            ## Std = np.std(CLmZeile[j-5:j+5])
            Mean = np.mean(CLmSmooth[j-5:j+5]) # Hic sunt leones!!!
            Std = np.std(CLmSmooth[j-5:j+5])   # Indexing difficult ....
            print 'i, j =', i, j
            print 'Mean, Std =', Mean, Std
            if np.isnan(Mean):
                print 'CLmSmooth[j-5:j+5] =', CLmSmooth[j-5:j+5]
            ## print 'CLmSmooth[j] =', CLmSmooth[j]

            # Actual "smoothing":
            if abs(CLmSmooth[j]-Mean) > 3.*Std: # Or any other value?
                CLmSmooth[j] = Mean
                print 'i, j =', i, j
                print 'Mean, Std =', Mean, Std
                print 'CLmZeile[j] =', CLmZeile[j]
                print 'CLmSmooth, corr =', CLmSmooth[j]


        # Interpolate and test smoothing results: 
        #f = interpolate.interp1d(sigmavTestRange, CLmins[i][:], kind=3)
        #f = interpolate.UnivariateSpline(sigmavTestRange, CLmins[i][:], s=1e13, k=4)
        f = interpolate.interp1d(sigmavTestRange, CLmSmooth, kind=3) 
        y = f(svtr2)

        # Unsmoothed CLmins:
        ax1.plot(sigmavTestRange, CLmins[i][:])

        # Smoothed values:
        #ax1.plot(sigmavTestRange, CLmSmooth)

        # Interpolated values:
        #ax1.plot(svtr2, y)#, label=r'$\tau\tau$')

        
    plt.show()

##########################################
# THE END
##########################################

