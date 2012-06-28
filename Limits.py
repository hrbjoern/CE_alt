# Function for calculating limits like TRolke 4.1.5
# HB from 2012-06-26

######################################################
# Importing:

import sys
from math import sqrt
try:
    import ROOT
    print "pyRoot loaded" 
except:
    print "Could not import pyRoot !"
    print "   Make sure you added root to your python path "
    print '   export PYTHONPATH="$ROOTSYS"/lib:"$PYTHONPATH"      '
    print "   Also, make sure your root installation contains pyRoot"
    print "   e.g. you should have the file $ROOTSYS/lib/libPyROOT.so"  
    sys.exit()

try:
    ROOT.gSystem.Load("libTRolke.1.so")
    print "TRolke library loaded" 
except:
    print "TRolke library not found"
    print "  ... this may not work, going on."


######################################################
# Define function: (cf. "LimitsTRolke4.1.5.py")

def UpperLimit(Non=160, Noff=122, alpha=1.0):
    """ Upper limit a la TRolke 4.1.5"""
    lower=ROOT.Double() 
    upper=ROOT.Double()  
    ncrit=ROOT.Long() 
    LDS=ROOT.Double()  

    # make trolke object named tr
    tr = ROOT.TROLKE2.TRolke()

    # Numbers:
    x = Non       # events in the signal region
    y = Noff      # events in background region
    tau = 1/alpha   # ratio of background/size region
    #print 'tau = ', tau
    alpha = 0.95  # Confidence Level
    #print 'CL = ', alpha
    e = 0.75       # efficiency
    #print 'eff = ', e
    # (Notation is (c) J. Lundberg..!)

    ## print 'Non = ', x
    ## print 'Noff = ', y
    ## print 'C.L. = ', alpha
    ## print 'alpha = ', tau
    ## print 'eff = ', e
    ## print 

    # Select model no. 5:
    # ------------
    tr.SetPoissonBkgKnownEff(x,y,tau,e)
    tr.SetCL(alpha)
    tr.SetBounding(True) # consider only physical values for Nsignal

    ## # Try model no. 4.1.2:
    ## bm = y/tau
    ## sdb = sqrt(y/tau)
    ## print 'bm = ', bm
    ## print 'sdb = ', sdb
    
    ## tr.SetGaussBkgKnownEff(x,bm,sdb,e)
    ## tr.SetCL(alpha)
    ## tr.SetBounding(True)
    
    # Calculate limits
    # ----------------
    status = 0
    status=tr.GetLimits(lower,upper)

    #print lower, upper
    
    # return upper (!) limit:
    # ------------
    if(status):
        return upper
        
    else:
        print "Calculation failed. (This should not be)"
        return 0


