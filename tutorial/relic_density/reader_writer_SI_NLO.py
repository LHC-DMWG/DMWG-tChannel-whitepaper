import cmath
import math
from pylab import *
from numpy import *
from matplotlib import *
import os
import os.path


## here the functions for SI@NLO

def ft():
    ftpu = 0.0153
    ftpd = 0.0191
    ftps = 0.0447
    ftnu = 0.0110
    ftnd = 0.0273
    ftns = 0.0447
    ftq = [ftpu, ftpd, ftps, ftnu, ftnd, ftns]
    return ftq



def momenta():
    g2 = 0.464
    u2p = 0.223
    ub2p = 0.036
    d2p = 0.118
    db2p = 0.037
    d2n = 0.223
    db2n = 0.036
    u2n = 0.118
    ub2n = 0.037
    s2 = 0.0258
    sb2 = 0.0258
    c2 = 0.0187
    cb2 = 0.0187
    b2 = 0.0117
    bb2 = 0.0117
    mmuz = [g2,u2p,ub2p,d2p,db2p,u2n,ub2n,d2n,db2n,s2,sb2,c2,cb2,b2,bb2]
    return mmuz

def massn():
    mp = 0.938
    mn = 0.939
    m_nucl = [mp,mn]
    return m_nucl

def summp():
    value = momenta()[1]+momenta()[2]+momenta()[3]+momenta()[4]+momenta()[9]+momenta()[10]+momenta()[11]+momenta()[12]+momenta()[13]+momenta()[14] 
    return value

def summn():
    value = momenta()[5]+momenta()[6]+momenta()[7]+momenta()[8]+momenta()[9]+momenta()[10]+momenta()[11]+momenta()[12]+momenta()[13]+momenta()[14] 
    return value

def ftgp():
    return  1.-(ft()[0]+ft()[1]+ft()[2])

def ftgn():
    return 1.-(ft()[3]+ft()[4]+ft()[5])

def sumlqp():
    value = ft()[0]+ft()[1]+ft()[2]
    return value

def sumlqn():
    value = ft()[3]+ft()[4]+ft()[5]
    return value

def red_mass(mdm):
    rm_p = 1./(mdm + massn()[0])
    rm_n = 1./(mdm + massn()[1])
    rm_nucl = [rm_p, rm_n]
    return rm_nucl


# majorana Dm and S3 mediator
def wcoeff_majDM(ac, bc, mdm, mmed):
    C_majDM_q =  (ac**2 + bc**2)/8. * mdm /(mmed**2 -mdm**2)**2
    C_majDM_qT1 = 0.5 * (ac**2 + bc**2) * mdm /(mmed**2 -mdm**2)**2
    C_majDM_qT2 = 0. 
    C_majDM_qAV = (ac**2 + bc**2)/4. * 1. /(mmed**2 -mdm**2)
    C_majDM_g = - mdm*(ac**2 + bc**2)/(96.*mmed**2*(mmed**2 -mdm**2)) 
    C_majDM_gT1 = 0.
    C_majDM_gT2 = 0.
    WilsCoeff_majDM = [C_majDM_q, C_majDM_qT1, C_majDM_qT2, C_majDM_qAV, C_majDM_g, C_majDM_gT1, C_majDM_gT2]
    return WilsCoeff_majDM


def ftot_nucl_majDM(ac,bc,mdm,mmed, NLO):
    ftotqp = wcoeff_majDM(ac, bc, mdm, mmed)[0] * ft()[0]
    ftotqn = wcoeff_majDM(ac, bc, mdm, mmed)[0] * ft()[3]
    if NLO:
        ftotgp = -8./9. * wcoeff_majDM(ac, bc, mdm, mmed)[4] * (1.-ft()[0])
        ftotgn = -8./9. * wcoeff_majDM(ac, bc, mdm, mmed)[4] * (1.-ft()[3])
        ftottqp = 3./4. * ( wcoeff_majDM(ac, bc, mdm, mmed)[1]  + wcoeff_majDM(ac, bc, mdm, mmed)[2]) *  (momenta()[1]+momenta()[2])
        ftottqn = 3./4. * ( wcoeff_majDM(ac, bc, mdm, mmed)[1]  + wcoeff_majDM(ac, bc, mdm, mmed)[2] ) * (momenta()[5]+momenta()[6])
        ftottg = - 3./4. * ( wcoeff_majDM(ac, bc, mdm, mmed)[5] + wcoeff_majDM(ac, bc, mdm, mmed)[6] ) * momenta()[0]
    else:
        ftotgp = 0.
        ftotgn = 0.
        ftottqp = 0.
        ftottqn = 0.
        ftottg = 0.
    
    ftot_majDM = [massn()[0] * (ftotqp + ftotgp + ftottqp + ftottg), massn()[1] * (ftotqn + ftotgn + ftottqn + ftottg)]
    return ftot_majDM

def sigma_SI_maj(ac,bc,mdm,mmed,NLO=False):
    value_p = 4. * ftot_nucl_majDM(ac,bc,mdm,mmed,NLO)[0]**2 * massn()[0]**2 * mdm**2 * red_mass(mdm)[0]**2 / pi
    value_n = 4. * ftot_nucl_majDM(ac,bc,mdm,mmed,NLO)[1]**2 * massn()[0]**2 * mdm**2 * red_mass(mdm)[1]**2 / pi
    value = [value_p, value_n]
    return value


## here read out the scan output of maddm and write on the shell the output below for the plots

mx,my,omegah2,si_lo_p,xe_si,si_lo_n,sd_p,pico_p = np.genfromtxt('scan_run_01.txt',unpack=True,usecols=[1,2,3,7,8,9,11,12])

gevm2tocm2 = 0.3894e-27

for i in range(len(mx)): 
    mdm = my[i]   ## dm mass
    mmed = mx[i]   ## med mass
    y = 1.      ## coupling from uberUFO
    a=0.5 * y      ## a and b lagrangian coeff following hisano
    b=0.5 * y
    mu_ss_SI = 0
    mu_ss_SD = 0
    if (mdm == mmed):
        mmed = 1e-10*mdm
        si_nlo_p = sigma_SI_maj(a,b,mdm,mmed, NLO=True)[0] * gevm2tocm2
        si_nlo_n = sigma_SI_maj(a,b,mdm,mmed, NLO=True)[1] * gevm2tocm2
    else:
        si_nlo_p = sigma_SI_maj(a,b,mdm,mmed, NLO=True)[0] * gevm2tocm2
        si_nlo_n = sigma_SI_maj(a,b,mdm,mmed, NLO=True)[1] * gevm2tocm2
    mu_ss_SI =  si_nlo_p / xe_si[i]   ## definition of signal strenght for SI
    mu_ss_SD =  sd_p[i] / pico_p[i]             ## definition of signal strenght for SD
    print mx[i], my[i], omegah2[i],si_nlo_p,xe_si[i],si_nlo_n,sd_p[i],pico_p[i],mu_ss_SI,mu_ss_SD
