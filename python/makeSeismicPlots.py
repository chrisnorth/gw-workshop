# wrapper to write seismic plots

import seismic as sc
import matplotlib.pyplot as plot
plot.ion()

file='../SeismicIsolation/LLOSeismicData.mat'
data=sc.read_data(file)

col='LLOSeismicData'
data=data[col][0]

i=24
data=data[i]
dataname=data[0][0]

freq=data[5]
ps=data[6]

# print dataname
# print freq
# print ps

plot.figure(1)
plot.clf()
plot.plot(freq,ps[:,2],lw=2)
plot.xlim(2e-2,100)
plot.xscale('log')
plot.yscale('log')
plot.xlabel('Frequency (Hz)')
plot.xlabel('Frequency (Hz)')
plot.ylabel('Amplitude Spectral Density (m Hz$^{-1/2}$)')
plot.title('LIGO Seismic Power Spectrum')
plot.grid()
plot.show()
plot.savefig('../SeismicIsolation/LIGOSeismicData_%s.png'%(str(dataname)))