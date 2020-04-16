# analyse GW waveforms

import pandas
import matplotlib.pyplot as plot
import numpy as np
from scipy.optimize import curve_fit
import astropy.constants as const
import astropy.units as u
import os

# from astropy.table import Table

# define sine wave to fit
def fn_sin(t,h0,f,ph):
    return h0 * np.sin(2*np.pi*f*t + ph)

# define functions
def calcdist(h,f,t,tc=0*u.s):
    K1=2.27e6*u.m/u.s
    dist=K1 / (h * f**2 * (tc-t))
    return(dist)

def calcmch(f,t,tc=0*u.s):
    K0=(5**(3/8)/(8*np.pi))*(const.c**3/const.G)**(5/8)
    mch=K0**1.6 * f**-1.6 * (tc-t)**-0.6
    return(mch)

def calcm1m2(mch,q):
    m2=mch * (q**3/(1+q))**0.2
    m1=q*m2
    return(m1,m2)

def getparams(dataFile,dataDir='./',fig=None,q=1.,fguess=10):
    # read in data
    data=pandas.read_csv(dataFile)
    try:
        datat=data['time (s)']
    except:
        datat=data['t']
    try:
        datah=data['strain']
    except:
        datah=data['hp']

    fitres={'file':dataFile}


    # set time for fitting
    dt_fit=0.5
    tmin=min(datat)
    infit=np.where(datat<tmin+dt_fit)[0]

    t_fit=datat[infit]
    h_fit=datah[infit]

    pguess=np.array([1e-22,fguess,np.pi])
    bounds=np.array([[1e-24,5,0],[1e-22,15,2*np.pi]])
    # res,cov=curve_fit(fn_sin,t_fit,h_fit,pguess,bounds=bounds)
    res,cov=curve_fit(fn_sin,t_fit,h_fit,pguess)
    t=np.mean(t_fit)*u.s
    h=abs(res[0])*u.dimensionless_unscaled
    f=res[1]*u.Hz

    fitres['t']=t.value
    fitres['h']=h.value
    fitres['f']=f.value

    # dist_m=K1 / (h * f**2 * (tc-t))
    dist_m=calcdist(h,f,t)
    mch_kg=calcmch(f,t)
    m1m2=calcm1m2(mch_kg,q)

    fitres['dist_m']=dist_m.value
    fitres['dist_ly']=dist_m.to(u.lightyear).value
    fitres['mch_kg']=mch_kg.value
    fitres['mch_Msun']=mch_kg.to(u.Msun).value
    fitres['m1_kg']=m1m2[0].value
    fitres['m2_kg']=m1m2[1].value
    fitres['m1_Msun']=m1m2[0].to(u.Msun).value
    fitres['m2_Msun']=m1m2[1].to(u.Msun).value
    fitres['phase']=res[2]


    plot.figure(fig)
    plot.clf()
    plot.plot(t_fit,h_fit)
    plot.plot(t_fit,fn_sin(t_fit,res[0],res[1],res[2]),ls='--')

    return(fitres)

# fit1=getparams('Dataset1_noise22.csv',fig=1,q=1.5)
# fit2=getparams('Dataset2_noise22.csv',fig=2,q=2)
# fit3=getparams('Dataset3_noise22.csv',fig=3,q=1.4)
# fit4=getparams('Dataset4_noise22.csv',fig=4,q=1)
# fit5=getparams('Dataset5_noise22.csv',fig=5,q=1,fguess=30)

fit1=getparams('Dataset1.csv',fig=1,fguess=10)
fit2=getparams('Dataset2.csv',fig=2)
fit3=getparams('Dataset3.csv',fig=3)
fit4=getparams('Dataset4.csv',fig=4)
fit5=getparams('Dataset5.csv',fig=5,fguess=30)

fits=pandas.DataFrame([fit1,fit2,fit3,fit4,fit5])