# code to plot graphs relevant to the Michelson interferometer
from numpy import linspace,pi,sin,mean
import matplotlib.pyplot as plot
from matplotlib.ticker import FormatStrFormatter
plot.ion()

def make_wave(xarr=None,ncyc=None,nsamp=None,relphase=0.,amp=1.):
    # plot ncyc sinewave cycles, with nsamp samples, with a phase offset of relphase*period
    try:
        xarr
        noX=True
    except:
        noX=False
    if noX:
        if ncyc==None: ncyc=5.
        if nsamp==None: nsamp=1000
        xarr=linspace(0.,ncyc,nsamp)
    yarr=amp*sin(2*pi*(xarr + relphase))
    return(xarr,yarr)


def sum_waves(wave1,wave2):
    assert sum(wave1[0]-wave2[0])==0,'xarr not the same'
    xarr=wave1[0]
    yarr=wave1[1]+wave2[1]
    return(xarr,yarr)

def add_wave_axes(fig,ncyc=None,wide=True):
    # define tick formats
    lformat=FormatStrFormatter(r'$%g\lambda$')
    nformat=FormatStrFormatter(r'$%g$')

    # add background axis (for labels)
    axbg=fig.add_axes([0,0,1,1],frameon=True,axisbg='#ffffff')
    axbg.set_xticks([])
    axbg.set_yticks([])

    # ad first input wave axis
    if wide:
        ax1=fig.add_axes([0.05,0.25,0.25,0.65],frameon=False)
        ax2=fig.add_axes([0.35,0.25,0.25,0.65],frameon=False)
        axsum=fig.add_axes([0.7,0.25,0.25,0.65],frameon=False)
    else:
        ax1=fig.add_axes([0.05,0.55,0.4,0.3],frameon=False)
        ax2=fig.add_axes([0.05,0.15,0.4,0.3],frameon=False)
        axsum=fig.add_axes([0.55,0.2,0.4,0.6],frameon=False)
    # ax1.axhline(0.,c='0.7')
    # ax1.set_xticks([])
    ax1.xaxis.set_major_formatter(lformat)
    ax1.yaxis.set_major_formatter(nformat)
    ax1.set_yticks([-1,0,1])
    ax1.set_ylim(-1.1,1.1)
    ax1.grid(axis='both',ls='-',color='0.7')
    ax1.tick_params('both',length=0)


    # ax2.axhline(0.,c='0.7')
    ax2.xaxis.set_major_formatter(lformat)
    ax2.yaxis.set_major_formatter(nformat)
    # ax2.set_xticks([])
    ax2.set_yticks([-1,0,1])
    ax2.set_ylim(-1.1,1.1)
    ax2.grid(axis='both',ls='-',color='0.7')
    ax2.tick_params('both',length=0)

    # add output wave axis
    # axsum.axhline(0,c='0.7')
    # axsum.set_xticks([])
    axsum.xaxis.set_major_formatter(lformat)
    axsum.yaxis.set_major_formatter(nformat)
    axsum.set_yticks([-2,-1,0,1,2])
    axsum.set_ylim(-2.1,2.1)
    axsum.grid(axis='both',ls='-',color='0.7')
    axsum.tick_params('both',length=0)

    # add + and = signs to background axis
    if wide:
        locs={'+':(0.325,0.5),'=':(0.65,0.5)}
    else:
        locs={'+':(0.25,0.5),'=':(0.5,0.5)}
    axbg.annotate('+',locs['+'],ha='center',va='center',fontsize='xx-large')
    axbg.annotate('=',locs['='],ha='center',va='center',fontsize='xx-large')
    return(axbg,ax1,ax2,axsum)

def makesumplot(offset=0.,amp=[1,1],fignum=None,showsum=True,title='',save=False,fileOut=None,wide=False):

    if fileOut==None:
        tstr=''
        wstr=''
        if title!='':
            tstr='_labelled'
        if wide:
            wstr='_wide'
        fileOut='../Interferometer/wave_sum_offset{0:.2f}{1}{2}.png'.format(offset,tstr,wstr)
        # else:
        #     fileOut='../Interferometer/wave_sum_offset%g_labelled.png'%offset
    wave1 = make_wave(xarr=None,ncyc=None,nsamp=None,amp=amp[0])
    wave2 = make_wave(xarr=wave1[0],relphase=offset,amp=amp[1])

    wave_sum=sum_waves(wave1,wave2)

    if wide:
        fig1=plot.figure(fignum,figsize=(10,2))
    else:
        fig1=plot.figure(fignum,figsize=(10,5))
    fig1.clf()
    (axbg,ax1,ax2,axsum)=add_wave_axes(fig1,wide=wide)

    ax1.plot(wave1[0],wave1[1],lw=2,c='r')
    ax2.plot(wave2[0],wave2[1],lw=2,c='b')
    axsum.plot(wave_sum[0],wave_sum[1],lw=2,c='g')
    if showsum and offset!=-1.:
        axsum.plot(wave1[0],wave1[1],lw=2,c='r',alpha=0.5,ls='-')
        axsum.plot(wave2[0],wave2[1],lw=2,c='b',alpha=0.5,ls='-')

    if title!=None:
        axbg.annotate(title,(0.5,0.99),fontsize='xx-large',ha='center',va='top')
    axbg.annotate('A',(0.175,0.01),fontsize='xx-large',ha='center',va='bottom',color='r')
    axbg.annotate('B',(0.475,0.01),fontsize='xx-large',ha='center',va='bottom',color='b')
    axbg.annotate('A+B',(0.825,0.01),fontsize='xx-large',ha='center',va='bottom',color='g')

    plot.show()

    if save:
        print(('Saving to %s'%(fileOut)))
        fig1.savefig(fileOut,dpi=300.)
