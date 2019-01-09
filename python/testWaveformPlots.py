# Wrapper script to make waveform plots

import waveform as wv
from numpy import array,zeros_like,linspace,arange,log10,pi,sqrt,zeros,ceil
from scipy import where
import matplotlib.pyplot as plot
from astropy import constants
plot.ion()

# set masses
m1=[5,10,15,20,30,40,50]
m2=[5,10,15,20,30,40,50]

K0=(5**(3./8.)/(8.*pi))*(constants.c**3 / constants.G)**(5./8.)
K1=2.27e6

# set distance (in Pc)
dist=1.e6
dist_m = dist*constants.pc

n2m=len(m1)*len(m2)
m1arr=zeros(n2m)
m2arr=zeros(n2m)
mcharr=zeros(n2m)

n=0
for i in range(len(m1)):
    for j in range(len(m2)):
        m1arr[n]=m1[i]
        m2arr[n]=m2[j]
        mcharr[n]=wv.calc_Mch(m1[i],m2[j])*constants.M_sun.value
        fileInFits='../Waveform/data/m1-%d-m2-%d.fits'%(m1[i],m2[j])
        fileIn='../Waveform/data/m1-%d-m2-%d.txt'%(m1[i],m2[j])
        # try:
        #     wv.read_fits(fileInFits)
        #     print('read from %s'%fileInFits)
        # except:
        #     wv.txt2fits(fileIn,fileInFits)
        #     print('converted %s to %s'%(fileIn,fileInFits))
        print (m1[i],m2[j],mcharr[n]/constants.M_sun.value)
        n+=1


m1x=m1[0]
m2x=m2[0]
fileInFits='../Waveform/data/m1-%d-m2-%d.fits'%(m1x,m2x)

print('reading data')
dataInAll=wv.read_fits(fileInFits)
dataIn=dataInAll[0:-1:10]
print('data read')
mch=wv.calc_Mch(m1[i],m2[j])*constants.M_sun.value
t0=-1
print('setting up arrays')
tarr=arange(ceil(min(dataIn['t'])),0,1.)
farr=zeros_like(tarr)
harr=zeros_like(tarr)
himarr=zeros_like(tarr)
hmodarr=zeros_like(tarr)
K0arr=zeros_like(tarr)
K1arr=zeros_like(tarr)
distarr_2=zeros_like(tarr)
mcharr_2=zeros_like(tarr)
print('computing K-values')
for i in range(len(tarr)):
    farr[i]=wv.getFreq(dataIn,tarr[i])
    # print 'freq at t=%g s: %g Hz'%(tarr[i],farr[i])

    harr[i]=wv.getHre(dataIn,tarr[i])
    himarr[i]=wv.getHim(dataIn,tarr[i])
    hmodarr[i]=wv.getModH(dataIn,tarr[i])

hmodarr2 = sqrt(harr**2 + himarr**2)
K0arr=farr * mch**(5./8.) * abs(tarr)**(3./8.) / K0
K1arr=farr**2 * abs(tarr) * hmodarr * dist_m / K1

distarr_2 = K1 / (farr**2 * abs(tarr) * hmodarr) / dist_m
mcharr_2 = K0**(8./5) * farr**(-8./5.) * abs(tarr)**(-3./5.) / mch

first_idx=where(distarr_2.value < 1.01)[0][0]
first_t = tarr[first_idx] * 0.9

    # print 'K0 at t=%g s: %g'%(tarr[i],K0arr[i])

print('plotting graphs')
plot.figure(1)
plot.clf()
plot.plot(tarr,farr)
plot.xlabel('Time (s)')
plot.ylabel('Frequency (Hz)')
plot.axvline(first_t)

plot.figure(2)
plot.clf()
plot.plot(tarr,hmodarr,'b-')
plot.plot(tarr,hmodarr2,'g--')
# plot.plot(tarr,harr,'r-')
plot.xlabel('Time (s)')
plot.ylabel('Hmod (strain)')
plot.axvline(first_t)

plot.figure(3)
plot.clf()
plot.plot(tarr,K0arr)
plot.xlabel('Time (s)')
plot.ylabel('K0 relative error')
plot.axvline(first_t)

plot.figure(4)
plot.clf()
plot.plot(tarr,K1arr)
plot.xlabel('Time (s)')
plot.ylabel('K1 relative error')
plot.axvline(first_t)

plot.figure(5)
plot.clf()
plot.plot(tarr,mcharr_2)
plot.xlabel('Time (s)')
plot.ylabel('Chirp mass relative error')
plot.axvline(first_t)

plot.figure(6)
plot.clf()
plot.plot(tarr,distarr_2)
plot.xlabel('Time (s)')
plot.ylabel('Distance relative error')
plot.ylim(0,5)
plot.axvline(first_t)

# plot.figure(4)
# plot.clf()
# plot.plot(log10(abs(tarr)),log10(farr)/((-3./8.)*log10(abs(tarr)))+(5./8.)*log10(mch.value))
# plot.xlabel('log(Time/s)')
# plot.ylabel('log(Frequency/Hz)')
# plot.ylim(0,10)
# plot.xscale('log')
# plot.yscale('log')