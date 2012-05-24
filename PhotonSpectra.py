from math import pi, exp


############################################################
# Spectra:
##############################################

def Bergstrom1998(E,mchi):
    if E>mchi:
        return 0.
    else:
        dNdE = (1./mchi) * 0.73 * (E/mchi)**-1.5 * exp(-7.8*E/mchi) 
        #    dNdE = 0.73 * (E/mchi)**-1.5 * exp(-7.8*E/mchi) 
    return dNdE


def Bringmann2008(E,mchi):
    if E>mchi:
        return 0.
    else:
        x = E/mchi
        eps = 80./mchi
        dndE = ( (1./mchi) * (alpha/pi) * 4.*(1-x+x**2)**2 /
                 ((1-x+eps/2.)*x)*(log(2.*(1-x+eps/2.)/eps) -0.5 +x -x**3 ) )
        if dndE<0:
            dndE=0.

        dndE = dndE+( (1./mchi)*0.73*(E/mchi)**-1.5 *exp(-7.8*E/mchi) )
    return dndE

def DGLAP2009(E,mchi):
    if E>mchi:
        return 0.
    else:
        x = E/mchi
## Cirelli / Bergstrom: me vs. mmu?
##         dndE = (1./mchi) * ( (alpha/(2.*pi)) *
##                              (-1 + np.log((4.*mchi**2/me**2)*(1.-x)))
##                              *(1.+(1.-x)**2)/x )
##         dndE = (1./mchi) * ( (alpha/(2.*pi)) *
##                              (-1 + np.log((4.*mchi**2/mmu**2)*(1.-x)))
##                              *(1.+(1.-x)**2)/x )

## Bovy:
        dndE = ( (alpha/pi) * ((1.+(1.-x)**2)/E) *
                 (-1. + np.log(4.*(1.-x)) - 2.*np.log(mmu/mchi) ) )
                 

        if dndE<0:
            dndE=0.
    return dndE



def Delta(E, mchi):
    if ( fabs(E-mchi)<1.1):
        return 1.
    else:
        return 0.001


    

##############################################
# Integrated annihilation spectra:
##############################################


def BergInt(mchi):
    return quad(lambda E: Bergstrom1998(E,mchi), Eth, 1.1*mchi)[0]
#    result = quad(lambda E: BergstromSmoothed(E,mchi), Eth, 2.*mchi)[0]
#    print ('\ntype(lambda E: BergstromSmoothed(E,mchi)) = ',
#           type(lambda E: BergstromSmoothed(E,mchi)), '\n')
#    result = quadrature(lambda E: BergstromSmoothed(E,mchi), Eth, 2.*mchi,
#                        vec_func=False)[0]

def BringInt(mchi):
    timeBring_i = time()
    result = quad(lambda E: Bringmann2008(E,mchi), Eth, 1.1*mchi,
                  limit=150,full_output=1)[0]
#        result = quad(lambda E: BringmannSmoothed(E,mchi), Eth, 2.*mchi,
    timeBring_f = time()
#    print 'BringmannIntegralTime= %.1f sec' % (timeBring_f-timeBring_i)
    return result

## def BringIntTest(mchi):
##     timeBring_i = time()
##     result = fixed_quad(lambda E: BringmannSmoothedTest(E,mchi), Eth, 2.*mchi)
##     #,                        vec_func=False)[0]
##     timeBring_f = time()
## #    print 'BringmannIntegralTime= %.1f sec' % (timeBring_f-timeBring_i)
##     return result

# CHB: Testing ...
#print '\n CHB: Testing ... Integration\n'
#print 'BringInt(2000.) = ', BringInt(2000.)
#print 'BergInt(2000.) = ', BergInt(2000.)
#BringIntTest(1000.)
#sys.exit()


# 'old' limit calculation: Now replaced by integration including Aeffs!
def BergSigmavUL(mchi):
    result = 8.*pi*FluxUL*mchi**2 / (Jdeltaomega*BergInt(mchi))
    return result

def BringSigmavUL(mchi):
    result = 8.*pi*FluxUL*mchi**2 / (Jdeltaomega*BringInt(mchi))
    return result


