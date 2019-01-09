# Wrapper script to make waveform plots

import waveform as wv
from numpy import zeros,linspace,pi
import matplotlib.pyplot as plot
from astropy import constants
from astropy.table import Table
plot.ion()

m1=[5,10,15,20,30,40,50]
m2=[5,10,15,20,30,40,50]

nCycLim=200.

n2m=len(m1)*len(m2)
m1arr=zeros(n2m)
m2arr=zeros(n2m)
mcharr=zeros(n2m)
n=0
for i in range(len(m1)):
    for j in range(len(m2)):
        m1arr[n]=m1[i]
        m2arr[n]=m2[j]
        mcharr[n]=wv.calc_Mch(m1[i],m2[j])
        n+=1

K0=(5**(3./8.)/(8.*pi))*(constants.c**3 / constants.G)**(5./8.)

t0arr=zeros(n2m)
t1arr=zeros(n2m)
f0arr=zeros(n2m)
f1arr=zeros(n2m)
ncarr=zeros(n2m)
for n in range(n2m):
    ncarr[n]=0
    t0arr[n]=0.
    print('Masses: (%g,%g) -> %g Msun'%(m1arr[n],m2arr[n],mcharr[n]))
    while ncarr[n]<nCycLim:
        t0arr[n]+=1.
        f0arr[n]=K0.value * (mcharr[n]*constants.M_sun.value)**(-5./8.) * t0arr[n]**(-3./8.)
        ncarr[n]=t0arr[n]*f0arr[n]
    print('  t=%g s: f0= %g Hz -> nc=%g'%(t0arr[n],f0arr[n],ncarr[n]))
    t1arr[n]=(nCycLim * (mcharr[n]*constants.M_sun.value)**(5./8.) / K0.value)**(8./5.)
    f1arr[n]=nCycLim**(-3./5.) * (mcharr[n]*constants.M_sun.value)**(-1) * (K0.value)**(8./5.)
    print('  t=%g s: f0= %g Hz'%(t1arr[n],f1arr[n]))

plot.figure(1)
plot.clf()
plot.plot(mcharr,f0arr,ls='None',marker='x',c='b')
plot.plot(mcharr,f1arr,ls='None',marker='o',c='g')
plot.xlabel('Chirp mass (Msun)')
plot.ylabel('Cutoff frequency (Hz)')

plot.figure(2)
plot.clf()
plot.plot(mcharr,t0arr,ls='None',marker='x',c='b')
plot.plot(mcharr,t1arr,ls='None',marker='o',c='g')
plot.xlabel('Chirp mass (Msun)')
plot.ylabel('Waveform length (t)')
plot.show()

tabOut=Table([m1arr,m2arr,mcharr,f1arr,t1arr],names=['m1','m2','M_ch','f_low','duration'])
tabOut.write('low-f_calcs.csv',format='ascii.csv')

