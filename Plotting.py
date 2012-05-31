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
        ax1.plot(energies,o.ul_bbbar,label=o.Name)

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
    ax1.set_ylim(1e-24,1e-20)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')
    
    # Plot stuff:
    for o in PubList:
        ax1.plot(o.energies, o.ul, label=o.legend)

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
    ax1.set_xlabel('E (GeV)')
    ax1.set_ylabel(r'A$_{\rm eff}$ (cm$^2$)')

    #Plot stuff:
    for o in ObjList:
        print o.Name
        ax1.plot(o.Aeffdata[:,0],o.Aeffdata[:,1],label=o.Name)
    plt.legend(loc=4)
    plt.show()


def plotSpectra():
    # Import spectra:
    from PhotonSpectra import tautau, bbbar
    # Data:
    mchi=1000.
    energies = np.linspace(30.,1010.,1000)
    tautaulist = []
    bbbarlist = []
    for e in energies:
        tautaulist.append(tautau(e,mchi))
        bbbarlist.append(bbbar(e,mchi))

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

    plt.legend(loc='lower left')
    plt.show()
