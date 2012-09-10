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
    #Pubs = False

    if Pubs:
        plotPubs()
    else:
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
    """Plot some combined profile likelihoods"""
    # Open pickle files:
    mchis = pickle.load(open('saveE.p'))
    sigmavTestRange = np.arange(-26, -18, (26.-18.)/len(mchis)) #cf. Exclusion.py!
    CLmins = pickle.load(open('saveCLmins.p'))

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
    ax1.set_ylim(200.,400.)
    ax1.set_xlim(-23.,-19)
    ax1.yaxis.grid(color='gray', linestyle='dashed')
    ax1.xaxis.grid(color='gray', linestyle='dashed')

    # Plot data from CLmins array:
    for i in range(len(CLmins[:][0])):
    #for i in (0,1):
        ## print i
        ## print 'sigmavTestRange: Length = ',  len(sigmavTestRange)
        ## print sigmavTestRange
        ## print 'CLmins[i][:] - Length = ', len(CLmins[i][:])
        ## print CLmins[i][:]
        #CLvector = CLmins[i][:]
        ax1.plot(sigmavTestRange, CLmins[i][:])#, label=r'$\tau\tau$')


#    ax1.plot(sigmavTestRange, CLmins[-1][:])

    plt.show()
