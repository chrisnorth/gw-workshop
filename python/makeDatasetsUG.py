#make datasets
from astropy.table import Table
from numpy import random
import os
import matplotlib.pyplot as plot

datasets={'1':'waveform_GW170104_10.0.fits',
    '2':'waveform_GW151226_10.0.fits',
    '3':'waveform_GW170729_10.0.fits',
    '4':'waveform_GW150914_10.0.fits',
    '5':'waveform_GW170817.fits'}

dirIn='full-data'
dirOut='dataUG'
for d in datasets:
    file=datasets[d]
    print('creating dataset {} from {}'.format(d,datasets[d]))
    tab=Table.read(os.path.join(dirIn,file))
    tab.remove_column('hc')
    tab.rename_column('t','time (s)')
    tab.rename_column('hp','strain')
    tab.write(os.path.join(dirOut,'Dataset{}_nonoise.fits'.format(d)),overwrite=True)
    tab.write(os.path.join(dirOut,'Dataset{}_nonoise.csv'.format(d)),format='ascii.csv',overwrite=True)
    noise24=(random.random(len(tab))-0.5)*1e-24
    tab24=Table([tab['time (s)'],tab['strain']])
    tab24['strain']=tab24['strain']+noise24
    tab24.write(os.path.join(dirOut,'Dataset{}_noise24.fits'.format(d)),overwrite=True)
    tab24.write(os.path.join(dirOut,'Dataset{}_noise24.csv'.format(d)),format='ascii.csv',overwrite=True)
    tab23=Table([tab['time (s)'],tab['strain']])

    noise23=(random.random(len(tab))-0.5)*1e-23
    tab23['strain']=tab23['strain']+noise23
    tab23.write(os.path.join(dirOut,'Dataset{}_noise23.fits'.format(d)),overwrite=True)
    tab23.write(os.path.join(dirOut,'Dataset{}_noise23.csv'.format(d)),format='ascii.csv',overwrite=True)

    noise22=(random.random(len(tab))-0.5)*1e-22
    tab22=Table([tab['time (s)'],tab['strain']])
    tab22['strain']=tab22['strain']+noise22
    tab22.write(os.path.join(dirOut,'Dataset{}_noise22.fits'.format(d)),overwrite=True)
    tab22.write(os.path.join(dirOut,'Dataset{}_noise22.csv'.format(d)),format='ascii.csv',overwrite=True)

    plot.figure(int(d)+5)
    plot.clf()
    plot.plot(tab['time (s)'][0:5000],tab['strain'][0:5000])
    plot.plot(tab23['time (s)'][0:5000],tab23['strain'][0:5000])
    plot.plot(tab24['time (s)'][0:5000],tab24['strain'][0:5000])