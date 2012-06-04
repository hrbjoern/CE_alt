########################################################
# Plotting: (for the moment)
########################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
from math import log, log10


def plotObjects():
    """Plot all the objects coming from ..?"""
    # Open pickle files:
    energies = pickle.load(open('saveE.p'))
    ObjList = pickle.load(open('saveObj.p'))

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel(r'<$\sigma$v>$_{\textrm{UL}}$ (cm$^{3}$s$^{-1}$)')
    ax1.set_xlabel('E (GeV)')
    ax1.set_ylim(1e-24,1e-20)
    ax1.set_xlim(1e2,1e5)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')
    
    #Plot stuff:
    for o in ObjList:
        if o.Name=="SculptorIso":
            ax1.plot(energies,o.ul_WW,label=o.Name)
            #print o.Aeffdata

        else:
            ax1.plot(energies,o.ul_bbbar,label=o.Name)
            #if o.Name=="Willman1V":
                #print o.ul_bbbar
                #print o.Aeffdata

    ## for o in PubList:
    ##     ax1.plot(o[:,0],o[:,1])

    #Plot Pubs?
    Pubs = True
    #Pubs = False

    if Pubs:
        plotPubs(Pubs=True)
    else:
        plt.legend(loc=3,ncol=2,prop=matplotlib.font_manager.FontProperties(size='small'))
        plt.show()
    
    return 0

def plotPubs(Pubs=False):
    """Plot the published limits for comparison"""

    # Open pickle files:
    PubList = pickle.load(open('savePub.p'))

    # Prepare plot:
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel(r'<$\sigma$v>$_{\textrm{UL}}$ (cm$^{3}$s$^{-1}$)')
    ax1.set_xlabel('E (GeV)')
    ax1.set_ylim(1e-25,1e-20)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')

    # Reset color map?
    #matplotlib.cm.get_cmap()
    #matplotlib.colors.Normalize()
    ax1.set_color_cycle(('b','g','r','c','m','y'))
    
    # Plot stuff:
    for p in PubList:
        ax1.plot(p.energies, p.ul, label=p.legend, linestyle='dashed')

    plt.legend(loc=3,ncol=2,prop=matplotlib.font_manager.FontProperties(size='small'))
    plt.show()


def plotAeffs():
    # Open pickle file:
    ObjList = pickle.load(open('saveObj.p'))
    print ObjList
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
        ax1.plot(o.Aeffdata[:,0],o.Aeffdata[:,1],label=o.Name)
    plt.legend(["MAGIC","VERITAS","HESS Sgr","HESS Scu"],loc=4)
    plt.show()


def plotSpectra():
    # Import spectra:
    from PhotonSpectra import tautau, bbbar, Bergstrom1998
    # Data:
    mchi=1000.
    energies = np.linspace(30.,1010.,1000)
    tautaulist = []
    bbbarlist = []
    Berglist = []
    for e in energies:
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
    ax1.set_title(r'Spectra for $m_\chi$ = 1000 GeV')

    plot1 = ax1.plot(energies, tautaulist, label=r'$\tau\tau$')
    plot2 = ax1.plot(energies, bbbarlist, label=r'$b\bar{b}$')
    plot3 = ax1.plot(energies, Berglist, label=r'$WW$')
    ## plot1 = ax1.plot(energies, energies**2*tautaulist, label=r'$\tau\tau$')
    ## plot2 = ax1.plot(energies, energies**2*bbbarlist, label=r'$b\bar{b}$')
    ## plot3 = ax1.plot(energies, energies**2*Berglist, label=r'$WW$')

    plt.legend(loc='lower left')
    plt.show()
