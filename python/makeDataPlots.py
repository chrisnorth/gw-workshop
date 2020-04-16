# Wrapper script to make waveform plots

import waveform as wv
import localisation as loc
from numpy import array,zeros_like,linspace,arange,log10,pi,sqrt,zeros,ceil,round,mean,random,floor,ceil
from scipy import where
import matplotlib.pyplot as plot
import matplotlib.patches as patches
from matplotlib.ticker import AutoMinorLocator
from astropy import constants
from astropy.table import Table,Column
import pandas as pd
import argparse
import os
import json

plot.ion()

parser=argparse.ArgumentParser(description='Plot waveform plots')
parser.add_argument('-r','--recrop',action='store_true',help='set to re-crop data')
parser.add_argument('--real',action='store_true',help='set to re-crop data')
parser.add_argument('-s','--save',action='store_true',help='set to save data')
parser.add_argument('-l','--latex',action='store_true',help='set to write latex')
parser.add_argument('-n','--noplot',action='store_true',help='set to not make plots')
parser.add_argument('--nonoise',action='store_true',help='set to not make plots')
parser.add_argument('-g','--gcse',action='store_true',help='set to make GCSE versions')
parser.add_argument('-x','--excel',action='store_true',help='set to output Excel files')

args=parser.parse_args()

recrop=args.recrop
save=args.save
realdata=args.real
latex=args.latex
noplot=args.noplot
nonoise=args.nonoise
gcse=args.gcse
excel=args.excel
filesuff=''
if realdata:
    filesuff+='_realdata'
if nonoise:
    filesuff+='_nonoise'
if gcse:
    filesuff+='_gcse'

# read in data
if realdata:
    jsonin=json.load(open('galaxies_real.json'))
else:
    jsonin=json.load(open('galaxies.json'))

H0=70. #km/s/MPc
c=3.e8 #m/s
datasets={'m1':[],
        'm2':[],
        'mchirp':[],
        'mratio':[],
        'dist_pc':[],
        'dist_Mpc':[],
        'labels':[],
        'angs':[],
        'galname':[],
        'redshift':[],
        'vel_kms':[],
        'H0_calc':[]}
for j in jsonin:
    m1=float(j['M1'])
    m2=float(j['M2'])
    datasets['m1'].append(m1)
    datasets['m2'].append(m2)
    if 'Mratio' in j:
        datasets['mratio']=float(j['Mratio'])
    else:
        datasets['mratio']=m1/m2
    datasets['mchirp'].append((m1*m2)**0.6/(m1+m2)**0.2)
    datasets['dist_Mpc'].append(float(j['Distance_Mpc']))
    datasets['dist_pc'].append(float(j['Distance_Mpc'])*1.e6)
    datasets['angs'].append(float(j['Angle']))
    datasets['labels'].append(j['Label'])
    datasets['galname'].append(j['Name'])
    datasets['redshift'].append(float(j['Redshift']))
    datasets['vel_kms'].append(float(j['Redshift'])*c/1.e3)
    datasets['H0_calc'].append((float(j['Redshift'])*c/1.e3)/float(j['Distance_Mpc']))

datasets['reltA']=zeros_like(datasets['angs'])
datasets['reltB']=zeros_like(datasets['angs'])
datasets['reltC']=zeros_like(datasets['angs'])

# initialise detector positions
dets={
    'A':{'xy':[0,0],'order':1,'col':'r'},
    'B':{'xy':[-0.2,0],'order':2,'col':'g'},
    'C':{'xy':[-0.1,0.2*sqrt(3)/2],'order':3,'col':'b'}
}
detpairs={'AB':{'d1':'A','d2':'B','manginv':True,'moff':0.2,'dtoff':0.3,'col':'r','darc':0.07},
        'AC':{'d1':'A','d2':'C','manginv':False,'moff':0.2,'dtoff':-0.3,'col':'g','darc':0.12},
        'BC':{'d1':'B','d2':'C','manginv':False,'moff':0.2,'dtoff':-0.3,'col':'b','darc':0.12}
        }

# do timing calculations
dt3={}
for p in detpairs:
    pair=detpairs[p]
    dt3[p]={'name_dt':r'$t_{%s}-t_{%s}$ [ms]'%(pair['d2'],pair['d1']),'dt':zeros_like(datasets['angs']),
        'ang1':zeros_like(datasets['angs']),'ang2':zeros_like(datasets['angs']),
        'name_ang1':r'$\theta_\mathrm{%s-%s}$ (1) [deg]'%(pair['d2'],pair['d1']),
        'name_ang2':r'$\theta_\mathrm{%s-%s}$ (2) [deg]'%(pair['d2'],pair['d1'])}

for i in range(len(datasets['angs'])):
    random.seed(1)
    ang=datasets['angs'][i]
    for p in detpairs:
        loc.init_dets(dets)
        loc.init_pairs(detpairs,dets)
        pair=detpairs[p]
        dt=loc.do_calc(pair,ang*pi/180.,1.5e4,3.e5)
        dt3[p]['dt'][i]=float('%.2f'%(pair['dt']*1e3))
        dt3[p]['ang1'][i]=float('%.1f'%pair['ang1'])
        dt3[p]['ang2'][i]=float('%.1f'%pair['ang2'])
        d1=dets[pair['d1']]
        d2=dets[pair['d2']]
        if datasets['relt'+d1['name']][i]==0:
            d1['relt']=random.normal(0,5e-3)
            datasets['relt'+d1['name']][i]=d1['relt']
        if datasets['relt'+d2['name']][i]==0:
            d2['relt']=d1['relt'] + pair['dt']
            datasets['relt'+d2['name']][i]=d2['relt']
    #     print(d1['name'],d1)
    #     print(d2['name'],d2)
        print(i,d1['name'],d1['relt']*1e3,d2['name'],d2['relt']*1e3,(d1['relt']-d2['relt'])*1e3)

assert len(datasets['m1'])==len(datasets['m2']) and \
    len(datasets['dist_pc'])==len(datasets['m1']),\
    'Incorrect sizes of m1 / m2 / dist_pc: %d / %d / %d'%\
        (len(datasets['m1']),len(datasets['m2']),len(datasets['dist_pc']))

# nd=2
nd=len(datasets['m1'])

K0=(5**(3./8.)/(8.*pi))*(constants.c**3 / constants.G)**(5./8.)
K1=2.27e6
# light year in metres
ly=9.46e15

datasets['mch']=zeros_like(datasets['m1'])
datasets['dist_m']=zeros_like(datasets['m1'])
datasets['mratio']=zeros_like(datasets['m1'])
datasets['fileIn']=['']*nd
datasets['filePng']=['']*nd
datasets['filePdf']=['']*nd
datasets['filePs']=['']*nd
datasets['fileEps']=['']*nd
datasets['data']=['']*nd

for i in range(nd):
    m1=datasets['m1'][i]
    m2=datasets['m2'][i]

    mch=wv.calc_Mch(m1,m2) * constants.M_sun.value
    datasets['mch'][i]=mch
    dist_m = 1.e6 * constants.pc
    datasets['dist_m'][i]=dist_m.value
    datasets['mratio'][i]=m1/m2
    datasets['filePng'][i]='Temp/Dataset_%s_waveforms%s.png'%(datasets['labels'][i],filesuff)
    datasets['filePdf'][i]='Temp/Dataset_%s_waveforms%s.pdf'%(datasets['labels'][i],filesuff)
    datasets['filePs'][i]='Temp/Dataset_%s_waveforms%s.ps'%(datasets['labels'][i],filesuff)
    datasets['fileEps'][i]='Temp/Dataset_%s_waveforms%s.eps'%(datasets['labels'][i],filesuff)

    if recrop:
        if (realdata):
            fileInFits='full-data/waveform_{}_10.0.fits'.format(datasets['galname'][i])
            fileInFitsCropped='full-data/waveform_{}_cropped_10.0.fits'.format(datasets['galname'][i])
            print('')
            print(('reading data for {} from {}'.format(datasets['galname'][i],fileInFits)))
            dataInAll=wv.read_fits(fileInFits)
            dataInAll.add_column(Column(sqrt(dataInAll['hp']**2 + dataInAll['hc']**2),name='modh'))
            t0=-10.1
            if min(dataInAll['t'])<t0:
                excidx=where(dataInAll['t']<t0)[0]
                dataInAll.remove_rows(excidx)
            # if min(dataInAll['t'])>50:
            #     tarr=arange(ceil(min(dataInAll['t'])),0,5)
            # else:
            #     tarr=arange(ceil(min(dataInAll['t'])),0,1)
            # farr=zeros_like(tarr)
            # hmodarr=zeros_like(tarr)
            # for t in range(len(tarr)):
            #     farr[t]=wv.getFreq(dataInAll,tarr[t])
            #     hmodarr[t]=wv.getModH(dataInAll,tarr[t])
            # distarr_t = K1 / (farr**2 * abs(tarr) * hmodarr) / dist_m.value
            # mcharr_t = K0**(8./5) * farr**(-8./5.) * abs(tarr)**(-3./5.) / mch
        else:
            fileInFits='../Datasets/data/m1-%d-m2-%d.fits'%(m1,m2)
            fileInFitsCropped='../Datasets/data/m1-%d-m2-%d_cropped.fits'%(m1,m2)
            print('')
            print(('reading data for [%d,%d] from %s'%(m1,m2,fileInFits)))
            datasets['fileIn'][i]==fileInFitsCropped
            dataInAll=wv.read_fits(fileInFits)
            print(('data read (%d samples / %.1fs)'%(len(dataInAll),abs(min(dataInAll['t'])))))
            if min(dataInAll['t'])>50:
                tarr=arange(ceil(min(dataInAll['t'])),0,5)
            else:
                tarr=arange(ceil(min(dataInAll['t'])),0,1)
            print(('calculating %d derived quantities'%(len(tarr))))
            farr=zeros_like(tarr)
            hmodarr=zeros_like(tarr)
            for t in range(len(tarr)):
                farr[t]=wv.getFreq(dataInAll,tarr[t])
                hmodarr[t]=wv.getModH(dataInAll,tarr[t])
            distarr_t = K1 / (farr**2 * abs(tarr) * hmodarr) / dist_m.value
            mcharr_t = K0**(8./5) * farr**(-8./5.) * abs(tarr)**(-3./5.) / mch
            firstidx=where(abs(1.-distarr_t) < 1)[0]
            if len(firstidx)>0:
                t0=tarr[firstidx[0]] * 0.8
                excidx=where(dataInAll['t']<t0)[0]
                if len(excidx>0):
                    print(('removing %d/%d rows'%(len(excidx),len(dataInAll))))
                    dataInAll.remove_rows(excidx)
        datasets['data'][i]=dataInAll
        dataInAll.write(fileInFitsCropped,format='fits',overwrite=True)
        # asdfasf

        # if 'modh' in dataInAll.colnames:
        #     hcol='modh'
        # else:
        #     hcol='hp'
        # plot.figure(i+1)
        # plot.clf()
        # plot.subplot(1,2,1)
        # plot.title('%d : %d'%(m1,m2))
        # plot.plot(dataInAll['t'],dataInAll[hcol])
        # plot.plot(dataInAll['t'],dataInAll['hc'])
        # plot.plot(dataInAll['t'],dataInAll['modh'])
        #
        # plot.subplot(1,2,2)
        # plot.plot(tarr,1.-distarr_t)
        # # plot.axvline(dataInAll['t'][0])
        # # plot.ylim(-1e3,1e3)
    else:
        if (realdata):
            fileInFits='full-data/waveform_{}_cropped_10.0.fits'.format(datasets['galname'][i])
            print('')
            print(('reading data for {} from {}'.format(datasets['galname'][i],fileInFits)))
        else:
            fileInFits='../Datasets/data/m1-%d-m2-%d_cropped.fits'%(m1,m2)
            print('')
            print(('reading data for [%d,%d] from %s'%(m1,m2,fileInFits)))
        datasets['fileIn'][i]==fileInFits
        datasets['data'][i]=wv.read_fits(fileInFits)

# reltplot=array([0.9,0.6,0.3,0.1])
if realdata:
    reltplot=0.1
else:
    reltplot=0.6
# np=len(reltplot)
tplotarr=zeros(nd)
fplotarr=zeros_like(tplotarr)
hplotarr=zeros_like(tplotarr)
# mchplotarr=zeros_like(tplotarr)
# mch5plotarr=zeros_like(tplotarr)
# m1plotarr=zeros_like(tplotarr)
# m2plotarr=zeros_like(tplotarr)
for i in range(nd):
    print(('Plotting %s [%d,%d]'%\
        (datasets['labels'][i],datasets['m1'][i],datasets['m2'][i])))
    data=datasets['data'][i]
    if realdata:
        hcol='hp'
    else:
        data=wv.set_dist(data,datasets['dist_pc'][i])
        hcol='h_re'
    # data=wv.set_dist(data,1.e6)
    data_noise=data.copy()
    if not nonoise:
        print('Adding noise to {}'.format(datasets['galname'][i]))
        if realdata:
            if gcse:
                sigma=1.e-23
            else:
                sigma=3.e-23
        else:
            sigma=1.e-21
        data_noise=wv.add_noise(data_noise,sigma=sigma,seed=i)
        if not 'data_noise' in datasets:
            datasets['data_noise']=[]
        datasets['data_noise'].append(data_noise)
    print(('Time range: %.1f:%.1f'%(data['t'][0],data['t'][-1])))
    tctrs=[0]*nd
    trngs=[[]]*nd
    trngsplot=[[]]*nd
    irngs=[[]]*nd
    irngsplot=[[]]*nd
    hctrs=[[]]*nd

    if realdata:
        idxctr=int(reltplot*len(data))
        if idxctr>=len(data):
            idxctr=len(data)-1
        tctr=data['t'][idxctr]
        fctr=wv.getFreq(data,tctr)
        hctr=wv.getModH(data,tctr)
        trng=[tctr , tctr + 10./fctr]
        trngplot=[tctr , tctr + 10./fctr]
    else:
        idxctr=int(reltplot*len(data))
        if idxctr>=len(data):
            idxctr=len(data)-1
        tctr=data['t'][idxctr]
        fctr=wv.getFreq(data,tctr)
        hctr=wv.getModH(data,tctr)
        trng=[tctr - 5./fctr , tctr + 5./fctr]
        trngplot=[tctr - 10./fctr , tctr + 10./fctr]
        if trng[1]>0.:
            fctr=wv.getFreq(data,0.01*data['t'][0])
            hctr=max(data[hcol])
            trng=[-8./fctr,2./fctr]
            trngplot=[-16./fctr,4./fctr]
        if trng[1]>=max(data['t']):
            trng[1]=data['t'][-2]
            trngplot[1]=data['t'][-2]
    tctr=mean([trng[0],trng[1]])
    fctr=wv.getFreq(data,tctr)
    hctr=wv.getModH(data,tctr)
    print(('Time centre [range]: %.2f  [%.2f:%.2f] (%.4g)'%(tctr,trng[0],trng[1],hctr)))
    irng=[wv.getIdx(data,trng[0]),wv.getIdx(data,trng[1])]
    irngplot=[wv.getIdx(data,trngplot[0]),wv.getIdx(data,trngplot[1])]
    trngs[i]=trng
    trngsplot[i]=trngplot
    tctrs[i]=tctr
    irngs[i]=irng
    irngsplot[i]=irngplot
    tmrg=0.03 * sqrt(datasets['mchirp'][i]/5)
    imrg=[wv.getIdx(data,-tmrg),wv.getIdx(data,tmrg)]
    hctrs[i]=hctr
    tplotarr[i]=-tctr
    fplotarr[i]=fctr
    hplotarr[i]=hctr

    if (excel):
        # output Excel files
        if gcse:
            yscale=1e22
        else:
            scales=-round(log10(max(data[hcol][irngs[i][0]:irngs[i][1]])))
            yscale_exp=scales
            yscale=10**yscale_exp
        if nonoise:
            df=pd.DataFrame({
                'time (s)':data['t'],
                'strain':data[hcol]*yscale
            })
        else:
            df=pd.DataFrame({
                'time (s)':data_noise['t'],
                'strain':data_noise[hcol]*yscale
            })
            # excelOut='../Datasets/Dataset_{}_noise.xlsx'.format(i+1)
        excelOut='../Datasets/Dataset_{}{}.xlsx'.format(i+1,filesuff)
        writer = pd.ExcelWriter(excelOut, engine='xlsxwriter')
        
        workbook=writer.book
        chartsheet=workbook.add_worksheet('Plots')
        worksheet=workbook.add_worksheet('Data')
        worksheet.write('A1','time (s)')
        worksheet.write_column('A2',df['time (s)'])
        worksheet.write('B1','strain')
        worksheet.write_column('B2',df['strain'])
        
        # make detector timestreams
        for d in ['A','B','C']:
            didx=wv.getIdx(data,datasets['relt'+d][i])-wv.getIdx(data,0)
            if didx==0:
                dh=df['strain']
            elif didx < 0:
                dh=zeros_like(df['strain'])
                dh[0:didx]=df['strain'][-didx-1:-1]
                # dh[didx:-1]=0;
            elif didx > 0:
                dh=zeros_like(df['strain'])
                # dh[0:didx]=0;
                dh[didx:-1]=df['strain'][0:-didx-1]
            df['det {} plot'.format(d)]=dh
        
        worksheet.write('C1','det A plot')
        worksheet.write_column('C2',df['det A plot'])
        worksheet.write('D1','det B plot')
        worksheet.write_column('D2',df['det B plot'])
        if not gcse:
            worksheet.write('E1','det C plot')
            worksheet.write_column('E2',df['det C plot'])
        
        #######################
        ymaxall=max(abs(df['strain']))
        ymaxall=10**(floor(log10(ymaxall)))*ceil(ymaxall/10**(floor(log10(ymaxall))))
        ymaxsl=max(abs(df['strain'][irngsplot[i][0]:irngsplot[i][1]]))
        ymaxsl=10**(floor(log10(ymaxsl)))*ceil(ymaxsl/10**(floor(log10(ymaxsl))))    
        ymaxmrg=max(abs(df['strain'][imrg[0]:imrg[1]]))
        ymaxmrg=10**(floor(log10(ymaxmrg)))*ceil(ymaxmrg/10**(floor(log10(ymaxmrg))))
        if ymaxmrg > 5:
            ytmrg = 2
        else:
            ytmrg=1
        
        
        worksheet.write_row('F1',['Slice','Min','Max','Boxx1','Boxx2','Boxx3','Boxx4','Boxx5','Boxy1','Boxy2','Boxy3','Boxy4','Boxy5'])
        worksheet.write_row('F2',[1,trngplot[0],trngplot[1],
            trngplot[0],trngplot[0],trngplot[1],trngplot[1],trngplot[0],
            -ymaxsl,ymaxsl,ymaxsl,-ymaxsl,-ymaxsl])
        worksheet.write_row('F3',[2,-tmrg,tmrg,
            -tmrg,-tmrg,tmrg,tmrg,-tmrg,
            -ymaxmrg,ymaxmrg,ymaxmrg,-ymaxmrg,-ymaxmrg])
        
        #######################
        
        if gcse:
            ylabel='Strain'
            lwsl=1
        else:
            ylabel='Strain (x 10^{:.0f})'.format(-yscale_exp)
            lwsl=0.25
        
        #######################
        
        
        chartall = workbook.add_chart({'type': 'scatter','subtype':'straight'})
        chartall.add_series({'name':'Waveform',
            'categories':['Data',1,0,len(df)+1,0],
            'values':['Data',1,1,len(df)+1,1],
            'line':   {'width': 0.25},
        })
        chartall.add_series({'name':'Slice 1 [τ={:.1f}]'.format(tctrs[i]),
            'categories':['Data',1,8,1,13],
            'values':['Data',1,13,1,17],
            'line':{'width':2,'color':'green','transparency':50}
        })
        chartall.add_series({'name':'Slice 2 [τ=0]',
            'categories':['Data',2,8,2,13],
            'values':['Data',2,13,2,17],
            'line':{'width':2,'color':'red','transparency':50}
        })
        chartall.set_x_axis({'name': 'Time (s)',
            'min':ceil(min(df['time (s)'])),'max':0.1,
            'crossing':ceil(min(df['time (s)'])),
            'num_format': '0.0',
            'major_unit': 1.0,'minor_unit': 0.1,
            'major_tick_mark':'outside','minor_tick_mark':'outside',
            'major_gridlines':{'visible':True,'line':{'width':1}},
            'minor_gridlines':{'visible':True,'line':{'width':1,'dash_type':'dash'}},
        })
        chartall.set_y_axis({'name': ylabel,'crossing':-ymaxall})
        chartall.set_title({'name':'Dataset #{:.0f}'.format(i+1)})
        chartall.set_legend({'position':'bottom','delete_series':[0]})
        
        chartsheet.insert_chart('A2', chartall)
        
        
        #######################
        
        chartslice = workbook.add_chart({'type': 'scatter','subtype':'straight'})
        chartslice.add_series({
            'categories':['Data',irngsplot[i][0]+1,0,irngsplot[i][1]+1,0],
            'values':['Data',irngsplot[i][0]+1,1,irngsplot[i][1]+1,1],
            'line':   {'width': lwsl},
        })
        chartslice.set_x_axis({'name': 'Time (s)',
            'min':trngplot[0],'max':trngplot[1],
            'crossing':trngplot[0],
            'num_format': '0.0',
            'major_unit': 0.1,'minor_unit': 0.01,
            'major_tick_mark':'outside','minor_tick_mark':'outside',
            'major_gridlines':{'visible':True,'line':{'width':1}},
            'minor_gridlines':{'visible':True,'line':{'width':1,'dash_type':'dash'}},
        })
        chartslice.set_y_axis({'name': ylabel,
            'min':-ymaxsl,'max':ymaxsl,
            'crossing':-ymaxsl,
        })
        chartslice.set_title({'name':'Slice 1 [τ={:.1f}]'.format(tctrs[i])})
        chartslice.set_legend({'none':True})
        chartsheet.insert_chart('A18', chartslice)
        
        #######################

        #######################

        chartmrgA = workbook.add_chart({'type': 'scatter','subtype':'straight'})
        chartmrgA.add_series({
            'categories':['Data',imrg[0]+1,0,imrg[1]+1,0],
            'values':['Data',imrg[0]+1+1,2,imrg[1]+1,2],
            'line':   {'width': 1,'color':'red'},
        })
        chartmrgA.set_x_axis({'name': 'Time (s)',
            'num_format': '0.00',
            'min':-tmrg,'max':tmrg,
            'crossing':-tmrg,
            'major_unit': 0.01,'minor_unit': 0.001,
            'major_tick_mark':'outside','minor_tick_mark':'outside',
            'major_gridlines':{'visible':True,'line':{'width':1}},
            'minor_gridlines':{'visible':True,'line':{'width':1,'dash_type':'dash'}},
        })
        chartmrgA.set_y_axis({'name': ylabel,
            'min':-ymaxmrg,'max':ymaxmrg,
            'major_unit':ytmrg,
            'crossing':-ymaxmrg
        })
        chartmrgA.set_size({'y_scale':0.7})
        chartmrgA.set_title({'name':'Slice 2 (Detector A) [τ=0]'})
        chartmrgA.set_legend({'none':True})
        chartsheet.insert_chart('A36', chartmrgA)
        
        #######################
        
        chartmrgB = workbook.add_chart({'type': 'scatter','subtype':'straight'})
        chartmrgB.add_series({
            'categories':['Data',imrg[0]+1,0,imrg[1]+1,0],
            'values':['Data',imrg[0]+1+1,3,imrg[1]+1,3],
            'line':   {'width': 1,'color':'green'},
        })
        chartmrgB.set_x_axis({'name': 'Time (s)',
            'num_format': '0.00',
            'min':-tmrg,'max':tmrg,
            'crossing':-tmrg,
            'major_unit': 0.01,'minor_unit': 0.001,
            'major_tick_mark':'outside','minor_tick_mark':'outside',
            'major_gridlines':{'visible':True,'line':{'width':1}},
            'minor_gridlines':{'visible':True,'line':{'width':1,'dash_type':'dash'}},
        })
        chartmrgB.set_y_axis({'name': ylabel,
            'min':-ymaxmrg,'max':ymaxmrg,
            'major_unit':ytmrg,
            'crossing':-ymaxmrg
        })
        chartmrgB.set_size({'y_scale':0.7})
        chartmrgB.set_title({'name':'Slice 2 (Detector A) [τ=0]'})
        chartmrgB.set_legend({'none':True})
        chartsheet.insert_chart('A48', chartmrgB)
        
        #######################
        
        if not gcse:
            chartmrgC = workbook.add_chart({'type': 'scatter','subtype':'straight'})
            chartmrgC.add_series({
                'categories':['Data',imrg[0]+1,0,imrg[1]+1,0],
                'values':['Data',imrg[0]+1+1,4,imrg[1]+1,4],
                'line':   {'width': 1,'color':'blue'},
            })
            chartmrgC.set_x_axis({'name': 'Time (s)',
                'num_format': '0.00',
                'min':-tmrg,'max':tmrg,
                'crossing':-tmrg,
                'major_unit': 0.01,'minor_unit': 0.001,
                'major_tick_mark':'outside','minor_tick_mark':'none',
                'major_gridlines':{'visible':True,'line':{'width':1}},
                'minor_gridlines':{'visible':True,'line':{'width':1,'dash_type':'dash'}},
            })
            chartmrgC.set_y_axis({'name': ylabel,
                'min':-ymaxmrg,'max':ymaxmrg,
                'major_unit':ytmrg,
                'crossing':-ymaxmrg
            })
            chartmrgC.set_size({'y_scale':0.7})
            chartmrgC.set_title({'name':'Slice 2 (Detector C) [τ=0]'})
            chartmrgC.set_legend({'none':True})
            chartsheet.insert_chart('A60', chartmrgC)
        
        #######################
        
        writer.save()
    
    if not noplot:
        fig=plot.figure(i+1,figsize=(8,11))
        fig.clf()
        if nonoise:
            lw=2
        else:
            lw=0.75
        if gcse:
            yscale=1e22
        else:
            scales=-round(log10(max(data[hcol][irngs[i][0]:irngs[i][1]])))
            yscale_exp=scales
            yscale=10**yscale_exp

        # plot full data frame
        axf=fig.add_subplot(3,1,1)
        # yscale_exp=-round(log10(max(data[hcol])))
        # yscale=10**yscale_exp
        axf.plot(data['t'],data[hcol]*yscale,lw=0.5,c='k')
        axf.set_xlim([0.99*min(data['t']),0.1])
        if gcse:
            axf.set_ylabel('Strain')
        else:
            axf.set_ylabel('Strain ($\mathrm{x} 10^{%g}$)'%(-yscale_exp))
        ylimf=max(abs(array(axf.get_ylim())))
        ylimf=1.1*ylimf/2.
        axf.set_ylim(-ylimf,ylimf)
        # axf.plot(-data_noise['t'][irngsplot[r][0]:irngsplot[r][1]] + datasets['reltB'][i],
        #     data_noise[hcol][irngsplot[r][0]:irngsplot[r][1]]*yscale + ylim*1.5,lw=lw)
        # axf.plot(-data_noise['t'][irngsplot[r][0]:irngsplot[r][1]] + datasets['reltC'][i],
        #     data_noise[hcol][irngsplot[r][0]:irngsplot[r][1]]*yscale - ylim*1.5,lw=lw)
        # axf.grid('on',which='major',c='0.7',lw=1,ls='-')
        axf.minorticks_on()
        # axf.grid('on',axis='x',which='minor',c='0.7',lw=1,ls=':')
        axf.xaxis.set_minor_locator(AutoMinorLocator(10))
        # axf.set_xlabel('Time before merger (s)')

        # plot timeslice
        axt=fig.add_subplot(3,1,2)
        # yscale_exp=-20
        # yscale_exp=-round(log10(max(data[hcol][irngs[r][0]:irngs[r][1]])))
        axt.plot(data_noise['t'][irngsplot[i][0]:irngsplot[i][1]],
            data_noise[hcol][irngsplot[i][0]:irngsplot[i][1]]*yscale,lw=lw)

        if gcse:
            axt.set_ylabel('Strain')
        else:
            axt.set_ylabel('Strain ($\mathrm{x} 10^{%g}$)'%(-yscale_exp))

        axt.set_xlim(trngs[i][0],trngs[i][1])
        ylimt=max(abs(array(axt.get_ylim())))*1.1
        # axt.set_ylim(-ylimt,ylimt)
        # axt.plot(-data_noise['t'][irngsplot[r][0]:irngsplot[r][1]] + datasets['reltB'][i],
        #     data_noise[hcol][irngsplot[r][0]:irngsplot[r][1]]*yscale + ylim*yfac1,lw=lw,c=dets['B']['col'])
        # axt.plot(-data_noise['t'][irngsplot[r][0]:irngsplot[r][1]] + datasets['reltC'][i],
        #     data_noise[hcol][irngsplot[r][0]:irngsplot[r][1]]*yscale - ylim*yfac1,lw=lw,c=dets['C']['col'])

        axt.grid('on',which='major',c='0.7',lw=1,ls='-')
        axt.minorticks_on()
        axt.grid('on',axis='x',which='minor',c='k',lw=1,ls=':')
        axt.xaxis.set_minor_locator(AutoMinorLocator(10))
        plot.title('Slice 1 [$\\tau$=%.1f s]'%(abs(tctrs[i])),color='k',fontweight='bold')
        # if r==np-1:
        #     axt.set_xlabel('Time before merger (s)')

        axt.set_ylim(-ylimt,ylimt)
        # axt.annotate(' Plot 1 [$\\tau$=%.1f s]'%\
        #     (abs(tctrs[i])),(trngs[i][0],-0.9*ylimt),color='g',fontweight='bold')

        # plot merger frame
        axm=fig.add_subplot(3,1,3)
        # yscale_exp=-round(log10(max(data[hcol])))
        # yscale=10**yscale_exp
        axm.plot((data['t'] + datasets['reltA'][i]),data_noise[hcol]*yscale,lw=lw,c=dets['A']['col'],label='Detector %s'%(dets['A']['name']))
        xrng=0.03 * sqrt(datasets['mchirp'][i]/5)
        axm.set_xlim([-xrng, xrng])
        axm.set_ylabel('Strain')
        ylimm=max(abs(array(axm.get_ylim())))
        # ylim=ylim/2.
        if gcse:
            axm.set_ylim(-ylimm*3,ylimm)
            axm.set_yticks([0,-ylimm*2])
            axm.set_yticklabels(['A','B'])
        else:
            axm.set_ylim(-ylimm*6,ylimm)
            axm.set_yticks([0,-ylimm*2,-ylimm*4])
            axm.set_yticklabels(['A','B','C'])
        axm.plot((data['t'] + datasets['reltB'][i]),data_noise[hcol]*yscale - ylimm*2,\
            lw=lw,c=dets['B']['col'],label='Detector %s'%(dets['B']['name']))
        axm.axhline(0,c=dets['A']['col'],lw=lw,alpha=0.5)
        axm.axhline(-ylimm*2,c=dets['B']['col'],lw=lw,alpha=0.5)
        if not gcse:
            axm.plot((data['t'] + datasets['reltC'][i]),data_noise[hcol]*yscale - ylimm*4,\
                lw=lw,c=dets['C']['col'],label='Detector %s'%(dets['C']['name']))
            axm.axhline(-ylimm*4,c=dets['C']['col'],lw=lw,alpha=0.5)

        # axm.grid('on',which='major',c='0.7',lw=1,ls='-')
        # axm.minorticks_on()
        # axm.annotate(' Merger [$\\tau$=0 s]',(-xrng,-5.9*ylimm),color='g',fontweight='bold')
        axm.grid('on',axis='x',which='minor',c='k',lw=1,ls=':')
        axm.xaxis.set_minor_locator(AutoMinorLocator(10))
        axm.set_xlabel('Time (s)')
        axm.set_title(' Slice 2 [$\\tau$=0 s]',color='k',fontweight='bold')
        plot.legend(loc='lower right',frameon=True,ncol=3)

        # add lines and hatching to full dataset
        axf.axvline(trngs[i][0],color='g')
        axf.axvline(trngs[i][1],color='g')
        axf.annotate(' Slice 1',(trngs[i][1],-0.9*ylimf),color='g',fontweight='bold')
        axf.add_patch(patches.Rectangle((trngs[i][0],-ylimf),trngs[i][1]-trngs[i][0],\
            2*ylimf,hatch='///',fill=False,color='g'))

        axf.axvline(xrng,color='g')
        axf.annotate('Slice 2  ',(-xrng,-0.9*ylimf),color='g',fontweight='bold',ha='right')
        axf.add_patch(patches.Rectangle((-xrng,-ylimf),2*xrng,\
            2*ylimf,hatch='///',fill=False,color='g'))

        # plot axis background
        axbg=fig.add_axes([0,0,1,1],frameon=False,facecolor='None')
        axbg.annotate('Dataset #%s'%(datasets['labels'][i]),(0.5,0.98),\
            ha='center',va='top',fontsize='x-large')
        axbg.annotate(r'Mass ratio: %.2f'%(datasets['mratio'][i]),(0.98,0.98),\
            ha='right',va='top',fontsize='large')
        if gcse:
            axbg.annotate('GCSE',(0.02,0.98),\
                ha='left',va='top',fontsize='x-large')

        axbg.set_xticks([])
        axbg.set_yticks([])
        # if r==0:
        #     plot.title()

        if save:
            plot.savefig(datasets['filePng'][i],dpi=300)
            plot.savefig(datasets['filePdf'][i])
            plot.savefig(datasets['filePs'][i])

        # plot.ylim(-hctrs[-1],hctrs[-1])

        plot.show()
    
mchplotarr = K0.value**(8./5) * fplotarr**(-8./5.) * tplotarr**(-3./5.) / constants.M_sun.value
mch5plotarr = mchplotarr**5.
mch5plotarr_round = zeros_like(mch5plotarr)
for m in range(len(mch5plotarr)):
    mch5plotarr_round[m]=round(mch5plotarr[m], -int(floor(log10(abs(mch5plotarr[m]))))+1)
rplotarr_SI = (K1 / (hplotarr * tplotarr * fplotarr**2.))
rplotarr_SI_disp = zeros_like(rplotarr_SI)
rplotarr = rplotarr_SI / (ly * 1e6)
for i in range(nd):
    r=rplotarr_SI[i]
    # r2=round(r / 10**int(log10(r)),2) * 10**int(log10(r))
    r2=float('%.3g'%(r))
    rplotarr_SI_disp[i]=r2
m1plotarr=zeros_like(tplotarr)
m2plotarr=zeros_like(tplotarr)
datasets['calctab']=[[]]*nd
for i in range(nd):
    calcTab=Table([arange(0,1)+1],names=['Panel'])
    mr=datasets['mratio'][i]
    m2plotarr[i] = (mch5plotarr[i] * (1. + mr) / mr**3)**(1./5.)
    m1plotarr[i] = m2plotarr[i] * mr
    # calcTab.add_column(Column(round(tplotarr[i,:],2),name=r'Time $(\tau)$ [s]'))
    # calcTab.add_column(Column(round(fplotarr[i,:],2),name=r'Frequency ($f$) [Hz]'))
    # calcTab.add_column(Column(round(1.e20*hplotarr[i,:],2),name=r'Amplitude ($h$) $\times 10^{20}$'))
    # # calcTab.add_column(Column(round(mch5plotarr[i,:],2),name=r'$M_{ch}^5$ [$M_\odot^5$]'))
    # calcTab.add_column(Column(round(mchplotarr[i,:],2),name=r'$M_{ch}$ [$M_\odot$]'))
    # calcTab.add_column(Column(round(m1plotarr[i,:],2),name=r'$M_{1}$ [$M_\odot$]'))
    # calcTab.add_column(Column(round(m2plotarr[i,:],2),name=r'$M_{2}$ [$M_\odot$]'))
    # calcTab.add_column(Column(rplotarr_SI_disp[i,:],name=r'$D$ [m]'))
    # calcTab.add_column(Column(round(rplotarr[i,:],2),name=r'$D$ [Mly]'))
    # datasets['calctab'][i]=calcTab
    # calcTab.write('New_Waveform-table-%s.tex'%datasets['labels'][i],format='latex',overwrite=True,\
    #     latexdict={'tablealign':'h','preamble': r'\begin{center}','tablefoot': r'\end{center}'})

if latex:
    latexFile='tex/New_Waveform_datasets%s.tex'%(filesuff)
    lf=open(latexFile,'w')
    hdr=\
'''\documentclass[14pt,a4paper]{extarticle}
\\usepackage{fullpage}
\\usepackage[margin=0in]{geometry}
\\usepackage{pdfpages}
\\usepackage{graphicx}
\\usepackage{epsfig}
\\usepackage{nopageno}
\\title{Datasets}
\\author{}
\date{}
\\begin{document}
'''
    lf.write(hdr)
    for i in range(nd):
        txt='\\begin{figure}\n\centering\n\includegraphics[]{%s}\n\end{figure}\n\n'%\
            (datasets['filePdf'][i])
        # txt='\\begin{figure}\n\centering\n\epsfig{width=\\textwidth,file=%s}\n\end{figure}\n\n'%\
        #     (datasets['fileEps'][i].replace('.eps',''))
        # txt='\includepdf[width=\\textwidth]{%s}\n\n'%(datasets['filePdf'][i])
        lf.write(txt)
    lf.write('\end{document}')
    lf.close()
    os.system('pdflatex -output-directory=Temp tex/New_Waveform_datasets%s.tex'%(filesuff))
    os.system('cp Temp/New_Waveform_datasets%s.pdf ../Datasets/'%(filesuff))

    # outpur Galaxies dataset
    tabGals=Table([datasets['galname'],datasets['redshift'],datasets['angs']],\
        names=['Galaxy Name','Redshift','Angle'])
    # write to CSV
    if realdata:
        tabGals.write('../Datasets/Galaxies-table_real.csv',format='ascii.csv',overwrite=True)

        tabGals.write('Temp/Galaxies-table-only_real.tex',format='latex',overwrite=True)
        os.system('pdflatex -output-directory=Temp tex/Galaxies-table_real.tex ')
        os.system('cp -f Temp/Galaxies-table_real.pdf ../Datasets/')


    tabAnsLoc1=Table([datasets['labels'],datasets['angs']],names=['Dataset',r'True $\theta_{\mathrm{GW}}$'])
    for p in ['AB']:
        tabAnsLoc1.add_column(Column(dt3[p]['dt'],name=dt3[p]['name_dt']))
        tabAnsLoc1.add_column(Column(dt3[p]['ang1'],name=dt3[p]['name_ang1']))
        tabAnsLoc1.add_column(Column(dt3[p]['ang2'],name=dt3[p]['name_ang2']))

    tabAnsLoc2=Table([datasets['labels'],datasets['angs']],names=['Dataset',r'True $\theta_{\mathrm{GW}}$'])
    for p in ['BC','AC']:
        tabAnsLoc2.add_column(Column(dt3[p]['dt'],name=dt3[p]['name_dt']))
        tabAnsLoc2.add_column(Column(dt3[p]['ang1'],name=dt3[p]['name_ang1']))
        tabAnsLoc2.add_column(Column(dt3[p]['ang2'],name=dt3[p]['name_ang2']))

    tabAnsDist=Table([datasets['labels']],names=['Dataset'])
    tabAnsDist.add_column(Column(round(tplotarr,2),name=r'Time $(\tau)$ [s]'))
    tabAnsDist.add_column(Column(round(fplotarr,2),name=r'Frequency ($f$) [Hz]'))
    tabAnsDist.add_column(Column(round(1.e22*hplotarr,2),name=r'Amplitude ($h$) $\times 10^{22}$'))
    tabAnsDist.add_column(Column(rplotarr_SI_disp,name=r'$D$ [m]'))
    tabAnsDist.add_column(Column(round(rplotarr,2),name=r'$D$ [Mly]'))

    tabAnsMass=Table([datasets['labels']],names=['Dataset'])
    tabAnsMass.add_column(Column(round(mch5plotarr_round,0),name=r'$M_{ch}^5$ [$M_\odot^5$]'))
    tabAnsMass.add_column(Column(round(mchplotarr,2),name=r'$M_{ch}$ [$M_\odot$]'))
    tabAnsMass.add_column(Column(round(datasets['mratio'],2),name='Mass ratio'))
    tabAnsMass.add_column(Column(round(m1plotarr,2),name=r'$M_{1}$ [$M_\odot$]'))
    tabAnsMass.add_column(Column(round(m2plotarr,2),name=r'$M_{2}$ [$M_\odot$]'))

    tabAnsGals=Table([datasets['labels'],datasets['galname'],datasets['angs'],\
            datasets['redshift'],round(datasets['vel_kms'],0),datasets['dist_Mpc'],round(datasets['H0_calc'],2)],\
        names=['Dataset','Galaxy Name','Angle',\
            'Redshift','Velocity (km/s)','Distance (Mpc)','H0 (km/s/MPc)'])
    if realdata:
        fReal='_real'
    else:
        fReal=''
    tabAnsLoc1.write('Temp/Datasets_Loc1-answer-table-only{}.tex'.format(fReal),format='latex',\
        latexdict={'col_align':"cc|ccc|ccc|ccc|",'tablealign':'!h'})
    tabAnsLoc2.write('Temp/Datasets_Loc2-answer-table-only{}.tex'.format(fReal),format='latex',\
        latexdict={'col_align':"cc|ccc|ccc|ccc|",'tablealign':'!h'})
    tabAnsDist.write('Temp/Datasets_Dist-answer-table-only{}.tex'.format(fReal),format='latex',\
        latexdict={'tablealign':'!h'})
    tabAnsMass.write('Temp/Datasets_Mass-answer-table-only{}.tex'.format(fReal),format='latex',\
        latexdict={'tablealign':'!h'})
    tabAnsGals.write('Temp/Datasets_Gals-answer-table-only{}.tex'.format(fReal),format='latex',\
        latexdict={'tablealign':'!h'})

    os.system('pdflatex -output-directory=Temp tex/Datasets_answer-table{}.tex'.format(fReal))
    os.system('cp -f Temp/Datasets_answer-table{}.pdf ../Datasets/'.format(fReal))


#     latexFileAnswers='New_Waveform_answers.tex'
#     lfa=open(latexFileAnswers,'w')
#     hdra=\
# '''\documentclass[14pt,a4paper]{extarticle}
# \\usepackage{fullpage}
# \\usepackage[margin=0.5in,landscape]{geometry}
# \
# \\usepackage{nopageno}
# \\title{Datasets}
# \\author{}
# \date{}
# \\begin{document}
# '''
#     lfa.write(hdra)
#     for i in range(nd):
#         txta='\section*{Dataset %s}\n\subsection*{$M_1=%d M_\odot$, $M_2 = %d M_\odot$, $M_1/M_2=%.2f$, $D=%.2f Mpc$}'%\
#             (datasets['labels'][i],datasets['m1'][i],datasets['m2'][i],datasets['mratio'][i],datasets['dist_pc'][i]/1.e6)
#         lfa.write(txta)
#         txta='\input{./Waveform-table-%s.tex}\n'%(datasets['labels'][i])
#         lfa.write(txta)
#         if (i % 2)!=0:
#             lfa.write('\\clearpage\n')
#
#     lfa.write('\end{document}')
#     lfa.close()
#     os.system('pdflatex %s'%(latexFileAnswers))
#     os.system('cp %s ../Datasets/'%(latexFileAnswers.replace('.tex','.pdf')))

