# Module containing scripts to make localisation plots

from numpy import pi,linspace,cos,sin,arctan2,sqrt,round,min,mean,abs,array,arange,log10,mod
import matplotlib.pyplot as plot
plot.ion()
import matplotlib
import waveform as wv
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

# plot.close('all')

# define conversion functions
def rth2xy(r,th):
    return(array([r*cos(th),r*sin(th)]))
def xy2rth(x,y):
    return(array([sqrt(x**2+y**2),arctan2(y,x)]))

def init_dets(dets):
    for d in dets:
        dets[d]['name']=d
        if not ('rth' in dets[d]):
            dets[d]['xy']=array(dets[d]['xy'])
            dets[d]['rth']=xy2rth(dets[d]['xy'][0],dets[d]['xy'][1])
        elif not ('xy' in dets[d]):
            dets[d]['rth']=array(dets[d]['rth'])
            dets[d]['xy']=rth2xy(dets[d]['rth'][0],dets[d]['rth'][1])
    return(dets)

def init_pairs(pairs,dets):
    for p in pairs:
        dd=pairs[p]
        dd['dets']=[dets[dd['d1']],dets[dd['d2']]]
        dd['ang']=arctan2(dd['dets'][0]['xy'][1]-dd['dets'][1]['xy'][1],
            dd['dets'][0]['xy'][0]-dd['dets'][1]['xy'][0])
        dd['sep']=sqrt((dd['dets'][1]['xy'][1]-dd['dets'][0]['xy'][1])**2 +
            (dd['dets'][1]['xy'][0]-dd['dets'][0]['xy'][0])**2)
    return(pairs)

def add_ax_pol(fig,fillfact=0.8,axisbg='#ffffff',frameon=True,ha='centre',va='centre'):
    figw=fig.get_figwidth()
    figh=fig.get_figheight()

    # axis is square
    axsize=fillfact*min([figw,figh])

    margin=(min([figw,figh])-axsize)/2.

    # set axis width relative to figure size
    axwidth=axsize/figw
    axheight=axsize/figh
    marw=margin/figw
    marh=margin/figh
    assert (ha=='centre' or ha=='left' or ha=='right'),'ha must be one of centre|left|right'
    assert (va=='centre' or va=='top' or va=='bottom'),'va must be one of centre|top|bottom'
    if ha=='centre':
        axleft=(1.-axwidth)/2.
    elif ha=='right':
        axleft=1.-(axwidth+marw)
    elif ha=='left':
        axleft=marw
    if va=='centre':
        axbottom=(1.-axheight)/2.
    elif va=='top':
        axbottom=1-(axheight+margin)
    elif va=='bottom':
        axbottom=margin
    # set up polar axes
    ax = fig.add_axes([axleft,axbottom,axwidth,axheight],
        polar=True, axisbg=axisbg)
    ax.set_ylim(0,1)

    # remove ytick labels
    yt=ax.get_yticks()
    ax.set_yticklabels([' ']*len(yt))
    # add more x ticks
    ax.set_thetagrids(linspace(0,360,24,endpoint=False),frac=1.07)

    # turn radial grid off
    ax.xaxis.grid()
    # add radial lines from r=0.2 outwards
    for xg in linspace(0,2*pi,24,endpoint=False):
        ax.plot([xg,xg],[0.2,1.1],'k:')

    # centre point
    ax.plot(0,0,ls='None',marker='x',ms=5,color='k')
    # add_ctr(ax)

    return(ax)

def add_ctr(fig,fillfact=0.8,ha='centre',va='centre'):
    ax_pol2=add_ax_pol(fig,fillfact=fillfact,axisbg='None',frameon=False,ha=ha,va=va)
    ax_pol2.set_xticks([])
    ax_pol2.set_yticks([])
    # ax_pol2.yaxis.grid()
    ax_pol2.plot(0,0,ls='None',marker='o',ms=10,mfc='w')
    ax_pol2.plot(0,0,ls='None',marker='x',ms=5,color='k')
    return()

def add_ax_xy(fig,fillfact=0.8,axisbg='None',frameon=False,ha='centre',va='centre'):
    figw=fig.get_figwidth()
    figh=fig.get_figheight()
    # set axis range to fill figure
    # axis size in inches
    axsize=fillfact*min([figw,figh])
    # margin size in inches
    margin=(min([figw,figh])-axsize)/2.
    # inch2units conversion
    in2u=2/axsize
    axw=2*figw/axsize
    axh=2*figh/axsize
    if ha=='centre':
        xrange=[-axw/2.,axw/2.]
    elif ha=='left':
        xrange=[-(1+margin*in2u),axw-(1+margin*in2u)]
    elif ha=='right':
        xrange=[-(axw-(1+margin*in2u)),1+margin*in2u]

    if va=='centre':
        yrange=[-figh/axsize,figh/axsize]

    # set up xy figure
    ax=fig.add_axes([0,0,1,1],
        frameon=frameon,axisbg=axisbg)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(xrange[0],xrange[1])
    ax.set_ylim(yrange[0],yrange[1])
    return(ax)

def add_ax_xy_inset(fig,fillfact=0.2,axisbg='None',frameon=False,ha='left',va='bottom',ctr=[0,0]):
    figw=fig.get_figwidth()
    figh=fig.get_figheight()
    # set axis range to fill figure
    # axis size in inches
    axsize=fillfact*min([figw,figh])
    # margin size in inches
    margin=(min([figw,figh])-axsize)/2.
    # inch2units conversion
    in2u=2/axsize
    axw=2*figw/axsize
    axh=2*figh/axsize

    # set up xy figure
    ax=fig.add_axes([0,0,fillfact,fillfact*figw/figh],
        frameon=frameon,axisbg=axisbg)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(ctr[0]-0.5,ctr[0]+0.5)
    ax.set_ylim(ctr[1]-0.5,ctr[1]+0.5)
    return(ax)

def add_det(det,ax,pol=False,col='b'):
    # add detectors to plot
    quad=int((det['rth'][1]+pi/4)%(2*pi)/(pi/2))+1
    if (quad==1):
        ha='left'
        va='center'
    elif quad==2:
        ha='center'
        va='bottom'
    elif quad==3:
        ha='right'
        va='center'
    elif quad==4:
        ha='center'
        va='top'
    if pol:
        ax.plot(det['rth'][1],det['rth'][0],ls='None',marker='x',
            ms=10,markeredgewidth=3,color='b',lw=3)
        ax.annotate(d,(det['rth'][1],det['rth'][0]+0.05),
            color='b',fontweight='bold',fontsize='large',ha=ha,va=va)
    else:
        ax.plot(det['xy'][0],det['xy'][1],ls='None',marker='x',
            ms=10,markeredgewidth=3,color='b',lw=3)
        labxy=rth2xy(det['rth'][0]+0.05,det['rth'][1])
        ax.annotate(det['name'],labxy,
            color='b',fontweight='bold',fontsize='large',ha=ha,va=va)

def add_det_pair(pair,ax,scale,col='0.7',lab=True,extlines=False):
    (mang,lang,junk)=getAngs(pair['ang'],offset=0,inv=pair['manginv'])
    pair['mang']=mang
    # get xy vector between pair
    pair_dxy=rth2xy(pair['sep'],pair['ang'])
    # get xy vector of line
    line_dxy=rth2xy(pair['moff'],mang)
    # get xy vector from centre to each detector
    dets_dxy=[]
    for det in pair['dets']:
        dets_dxy.append(det['xy'])
    mid_dxy=array([mean([pair['dets'][0]['xy'][0],pair['dets'][1]['xy'][0]]),
            mean([pair['dets'][0]['xy'][1],pair['dets'][1]['xy'][1]])])

    # ax.plot([dets_dxy[0][0],dets_dxy[0][0]+line_dxy[0]],
    #     [dets_dxy[0][1],dets_dxy[0][1]+line_dxy[1]],c='m',lw=3)

    # set arrow and text offsets
    moff=pair['moff']
    moffxy=array([moff*cos(mang),-moff*sin(mang)])
    mr=sqrt(moffxy[0]**2+moffxy[1]**2)

    for det in pair['dets']:
        # add arrow
        arr0=(1-0.05/mr)*line_dxy+mid_dxy
        darr=det['xy']-mid_dxy
        ax.arrow(arr0[0],arr0[1],darr[0],darr[1],
            length_includes_head=True,
            head_length=0.02,head_width=0.02,lw=2,color=col)
        # add measure line
        ax.plot([det['xy'][0],det['xy'][0]+line_dxy[0]],
                [det['xy'][1],det['xy'][1]+line_dxy[1]],
                ls='-',lw=1,c=col)

    # add connecting line
    ang=arctan2(pair['dets'][1]['xy'][1]-pair['dets'][0]['xy'][1],\
        pair['dets'][1]['xy'][0]-pair['dets'][0]['xy'][0])
    ax.plot([pair['dets'][0]['xy'][0],pair['dets'][1]['xy'][0]],
        [pair['dets'][0]['xy'][1],pair['dets'][1]['xy'][1]],
        lw=2,c=pair['col'],ls='--')
    if extlines:
        ax.plot([pair['dets'][0]['xy'][0],pair['dets'][0]['xy'][0]-sin(ang+pi/2)],\
            [pair['dets'][0]['xy'][1],pair['dets'][0]['xy'][1]+cos(ang+pi/2)],\
            lw=2,c=pair['col'],ls='--')
        ax.plot([pair['dets'][1]['xy'][0],pair['dets'][1]['xy'][0]-sin(ang-pi/2)],\
            [pair['dets'][1]['xy'][1],pair['dets'][1]['xy'][1]+cos(ang-pi/2)],\
            lw=2,c=pair['col'],ls='--')

    # add label
    toff=(1+0.05/mr)*moff
    tang=pair['ang']+pair['mang']+pi/2
    toffxy=rth2xy(toff,mang)

    txt='%dkm'%(int(pair['sep']*scale))
    ax.annotate(txt,(toffxy[0]+mid_dxy[0],toffxy[1]+mid_dxy[1])
        ,ha='center',va='center',color=col,rotation=lang*180./pi,
        fontweight='bold',fontsize='large')

def add_det_lines(pair,dets,ax):
    dy=dets[pair['d1']]['xy'][1]-dets[pair['d2']]['xy'][1]
    dx=dets[pair['d1']]['xy'][0]-dets[pair['d2']]['xy'][0]
    ang=arctan2(dy,dx)-(pi/2.)
    label='Detectors %s-%s'%(pair['d1'],pair['d2'])
    ax.plot([-sin(ang),sin(ang)],[cos(ang),-cos(ang)],'--',c=pair['col'],lw=3)
    # ax.plot([0,0.5*sin(-ang+pi),0.5*sin(-ang+pi)+0.05*sin(-ang+pi/2)],[0,0.5*cos(-ang+pi),0.5*cos(-ang+pi)+0.05*cos(-ang+pi/2)],'k-')
    textrot=mod(ang*180/pi+90,180)
    if textrot>90:
        textrot+=180
    ax.annotate(label,(0.5*sin(-ang+pi)+0.05*sin(-ang+pi/2),0.5*cos(-ang+pi)+0.05*cos(-ang+pi/2)),\
        rotation=textrot,\
        color=pair['col'],fontweight='bold',fontsize='large',ha='center',va='center')

def getAngs(ang,offset=0,inv=False):
    quad=int((ang+offset)%(2*pi)/(pi/2))+1
    if inv:
        s=-1
    else:
        s=1
    if (quad==1):
        mang=ang+s*pi/2
        swlang=ang
        gwlang=ang+pi/45
    elif (quad==2):
        mang=ang-s*pi/2
        swlang=ang+pi
        gwlang=ang-pi/45
    elif (quad==3):
        mang=ang-s*pi/2
        swlang=ang+pi
        gwlang=ang-pi/45
    else:
        mang=ang+s*pi/2
        swlang=ang
        gwlang=ang+pi/45
    return(mang,swlang,gwlang)

def add_sinwave(ang,ax,col='g'):
    (mang,swlang,gwlang)=getAngs(ang,offset=0)

    # add sin wave
    sw0=0.4
    swx=linspace(sw0,1,500)
    swy=0.02*sin((2*pi*(swx-sw0)/0.07))
    pswx=swx*cos(ang) - swy*sin(ang)
    pswy=swx*sin(ang) + swy*cos(ang)
    # plot wave
    ax.plot(pswx,pswy,lw=2,c=col)
    # plot arrow
    ax.arrow(sw0*cos(ang),sw0*sin(ang),0.5*sw0*cos(ang+pi),0.5*sw0*sin(ang+pi),
        color=col,lw=2)
    # plot dashed line
    ax.plot([0.5*sw0*cos(ang),-cos(ang)],[0.5*sw0*sin(ang),-sin(ang)],
        color=col,ls='--',lw=2)
    # add "Gravitational Wave" label

    gwlxy=rth2xy(0.75,gwlang)
    ax.annotate('Gravitational Wave',gwlxy,ha='center',va='center',color=col,
        rotation=swlang*180./pi,fontweight='bold')

def add_linepair(pair,ax,wang,scale,c,col='r',txt=None,outerArrows=True,alen=0.1,
            fontsize='x-large'):

    moff=pair['dtoff']

    (mang,lang,junk)=getAngs(wang,offset=0)
    # (mang,lang,junk)=getAngs(pair['ang'],offset=0,inv=True)
    # get xy vector between pair
    pair_dxy=rth2xy(pair['sep'],pair['ang'])
    # get xy vector of line
    line_dxy=rth2xy(pair['moff'],mang)
    # get xy vector from centre to each detector
    dets_dxy=[]
    for det in pair['dets']:
        dets_dxy.append(det['xy'])

    for det in pair['dets']:
        # calculate nearest point on wave path
        # mr0=det['rth'][0]*cos(mang)
        det['mxy0']=rth2xy(det['rth'][0]*cos(det['rth'][1]-wang),wang)
    mid_mxy0=array([mean([pair['dets'][0]['mxy0'][0],pair['dets'][1]['mxy0'][0]]),
            mean([pair['dets'][0]['mxy0'][1],pair['dets'][1]['mxy0'][1]])])
    for det in pair['dets']:
        mxy0=det['mxy0']
        d_mxy0=mxy0-mid_mxy0
        mr0=sqrt(d_mxy0[0]**2+d_mxy0[1]**2)
        # calculate end of label marker
        mxy1=rth2xy(moff,mang)
        mr1=sqrt(mxy1[0]**2+mxy1[1]**2)

        # add line from wave
        ax.plot([mxy0[0],mxy0[0]+mxy1[0]],[mxy0[1],mxy0[1]+mxy1[1]],
            c=col,lw=1)
        # add line from wave to det
        ax.plot([mxy0[0],det['xy'][0]],[mxy0[1],det['xy'][1]],
            c=col,lw=1)
        # add arrow
        # print mxy0,mxy1
        if outerArrows:
            arr0=mxy1*(1-0.05/mr1)+mid_mxy0+(1.+(alen/mr0))*d_mxy0
            darr=-(alen/mr0)*d_mxy0
        else:
            arr0=[mxy1[0]*(1-0.1/mr1),mxy1[1]*(1-0.1/mr1)]
            darr[mxy0[0],mxy0[1]]
        ax.arrow(arr0[0],arr0[1],darr[0],darr[1],
                length_includes_head=True,color=col,lw=2)

    # add \Delta t label
    lxy=mid_mxy0 + rth2xy((1+0.05/mr1)*moff,mang)
    dt=do_calc(pair,wang,scale,c)
    if txt==None:
        txt=r'$t_{%s}-t_{%s}=%.2f\,\mathrm{ms}$'%(pair['d1'],pair['d2'],dt*1000)
    ax.annotate(txt,lxy,ha='center',va='center',color=col,
        rotation=lang*180./pi,fontweight='bold',fontsize=fontsize)

def do_calc(pair,swang,scale,c=None):
    if c==None:
        c==3.e5
    dt=pair['sep']*scale*cos(swang-pair['ang'])/c
    pair['dt']=dt
    dang=abs(swang-pair['ang'])
    if dang>pi:
        dang=2*pi-dang
    pair['ang1'] = ((pair['ang'] + dang)*180./pi)%(360.)
    pair['ang2'] = ((pair['ang'] - dang)*180./pi)%(360.)
    return(dt)

def add_pair_angles(pair,ax,swang):
    dang=abs(swang-pair['ang'])
    if dang>pi:
        dang=2*pi-dang
    ax.plot([pair['dets'][0]['xy'][0],pair['dets'][1]['xy'][0]],
        [pair['dets'][0]['xy'][1],pair['dets'][1]['xy'][1]],
        lw=2,c=pair['col'],ls='-',alpha=0.5)
    for s in [+1,-1]:
        angplot=pair['ang'] + s*dang
        linexy = rth2xy(1,angplot)
        ax.plot([0,linexy[0]],[0,linexy[1]],c=pair['col'],ls='--',lw=2,alpha=0.5)
        labang=pair['ang']+s*dang/2
        labtxt=r'$\Delta\theta_{\mathregular{%s-%s}}$'%(pair['d1'],pair['d2'])
        ax.annotate(labtxt,pair['dets'][0]['xy']+rth2xy(pair['darc']+0.02+0.01*s,labang),
            ha='center',va='center',rotation=(labang-pi/2)*180/pi,fontsize='x-small',color=pair['col'])
        if pair['ang']<angplot:
            dth=arange(pair['ang'],angplot,0.5*pi/180)
        else:
            dth=arange(angplot,pair['ang'],0.5*pi/180)
        arcx=(pair['darc']+0.01*s)*cos(dth)
        arcy=(pair['darc']+0.01*s)*sin(dth)
        ax.plot(arcx,arcy,c=pair['col'])
    pair['ang1'] = ((pair['ang'] + dang)*180./pi)%(360.)
    pair['ang2'] = ((pair['ang'] - dang)*180./pi)%(360.)
    return()

def add_dt_text(detpairs,ax,swang,scale,c,loc="upper left"):
    text=''
    for p in detpairs:
        pair=detpairs[p]
        pid=r'\mathregular{%s-%s}'%(pair['d1'],pair['d2'])
        dt=do_calc(pair,swang,scale,c)

        text=text+r'Detector pair $%s$:'%(pid)+'\n'
        text=text+r'  Pair Separation $\Delta R_{%s}=%gkm$'%(pid,pair['sep']*scale)+' \n'
        text=text+r'  Time delay $\Delta t_{%s}$:'%(pid)+' \n'
        text=text+r'    $\Delta t_{%s} =  t_{\mathregular{%s}}-t_{\mathregular{%s}}$:'%\
            (pid,pair['d1'],pair['d2'])+'\n'
        text=text+r'    $\Delta t_{%s}=%.2f\,\mathrm{ms}$'%(pid,(pair['dt']*1000))+' \n\n'
    # print text
    if loc=='upper left':
        va='top'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='upper right':
        va='top'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower left':
        va='bottom'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower right':
        va='bottom'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.annotate(text,[locx,locy],va=va,ha=ha)

def add_calc_text(detpairs,ax,swang,scale,c,loc="upper left"):
    text=r'GW direction: $\theta_{GW}=%.1f^{\circ}$'%(swang*180/pi)+' \n\n'
    for p in detpairs:
        pair=detpairs[p]
        pid=r'\mathregular{%s-%s}'%(pair['d1'],pair['d2'])
        dt=do_calc(pair,swang,scale,c)
        dang=abs(swang-pair['ang'])
        if dang>pi:
            dang=2*pi-dang
        # pair['dt']=dt
        # print pair,pair['dt']*1000
        # text=text+'\n'
        text=text+r'Detector pair $%s$:'%(pid)+'\n'
        text=text+r'  Pair Separation $\Delta R_{%s}=%g\,\mathrm{km}$'%(pid,pair['sep']*scale)+' \n'
        text=text+r'  Pair Orientation $\theta_{%s}=%.1f^{\circ}$'%(pid,pair['ang']*180/pi)+' \n'
        text=text+r'  Relative angle $\Delta\theta_{%s}=%.1f^{\circ}$'%(pid,abs(dang)*180/pi)+' \n'
        text=text+r'  Time delay $\Delta t_{%s}$:'%(pid)+'\n'
        text=text+r'    $\Delta t_{%s}=\Delta R_{%s}\ \cos(\Delta\theta_{%s})/c$'%\
            (pid,pid,pid)+' \n'
        text=text+r'    $\Delta t_{%s}=%.2f\,\mathrm{ms}$'%\
            (pid,(pair['dt']*1000))+' \n\n'
    # print text
    if loc=='upper left':
        va='top'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='upper right':
        va='top'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower left':
        va='bottom'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower right':
        va='bottom'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.annotate(text,[locx,locy],va=va,ha=ha)

def add_invcalc_text(detpairs,ax,swang,scale,c,loc="upper left"):
    text=''
    for p in detpairs:
        pair=detpairs[p]
        pid=r'\mathregular{%s-%s}'%(pair['d1'],pair['d2'])
        dt=do_calc(pair,swang,scale,c)
        dang=abs(swang-pair['ang'])
        if dang>pi:
            dang=2*pi-dang
        ang1=pair['ang']+dang
        ang2=pair['ang']-dang
        # pair['dt']=dt
        # print pair,pair['dt']*1000
        # text=text+'\n'
        text=text+r'Detector pair $%s$:'%(pid)+'\n'
        text=text+r' Pair Separation $\Delta R_{%s}=%g\,\mathrm{km}$'%(pid,pair['sep']*scale)+' \n'
        text=text+r' Pair Orientation $\theta_{%s}=%.1f^{\circ}$'%(pid,pair['ang']*180/pi)+' \n'
        text=text+r' Time delay $\Delta t_{%s}=%.2f\,\mathrm{ms}$'%(pid,(pair['dt']*1000))+' \n'
        text=text+r' Relative angle $\Delta\theta_{%s}$:'%(pid)+'\n'
        text=text+r'   $\Delta\theta_{%s}=\cos^{-1}(c\Delta t_{%s}/\Delta R_{%s})$'%(pid,pid,pid)+' \n'
        text=text+r'   $\Delta\theta_{%s}=\pm%.1f^{\circ}$'%(pid,abs(dang)*180/pi)+' \n'
        text=text+r' GW direction: $\theta_{GW}=%.1f^{\circ}$ or $%.1f^{\circ}$'%\
            ((ang1*180/pi)%(360.),(ang2*180/pi)%(360.))+' \n\n'
    if len(list(detpairs.keys()))>1:
        text=text+r'GW direction: $%.1f^{\circ}$'%(swang*180/pi)

    # print text
    if loc=='upper left':
        va='top'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='upper right':
        va='top'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower left':
        va='bottom'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower right':
        va='bottom'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.annotate(text,(locx,locy),va=va,ha=ha)

def add_title(text,ax,loc='upper left'):
    if loc=='upper left':
        va='top'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='upper right':
        va='top'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.95*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower left':
        va='bottom'
        ha='left'
        locx=ax.get_xlim()[0] + 0.05*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    if loc=='lower right':
        va='bottom'
        ha='right'
        locx=ax.get_xlim()[0] + 0.95*(ax.get_xlim()[1] - ax.get_xlim()[0])
        locy=ax.get_ylim()[0] + 0.05*(ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.annotate(text,(locx,locy),va=va,ha=ha,fontsize='x-large')
    return()


################################################################################
################################################################################
################################################################################

def makeplot(dets,detpairs,waveang,scale=3e4,title=None,fignum=None,fileOut=None,save=False,showcalc=False):
    # draw detectors on

    swang=waveang*pi/180.

    cols={'det':'b','pair':'k','wave':'g','dt':'r','calc':'k'}
    # speed of light in km/s
    c=3.e5

    if fileOut==None:
        if showcalc:
            fileOut='Temp/%ddet_plot_withcalc.pdf'%(len(list(dets.keys())))
        else:
            fileOut='Temp/Outputs/%ddet_plot.pdf'%(len(list(dets.keys())))

    init_dets(dets)

    # define detector-pairs
    init_pairs(detpairs,dets)

    # initialise figure
    fig1=plot.figure(fignum,figsize=(11,8))
    fig1.clf()

    # add axes
    if showcalc:
        ha='right'
    else:
        ha='centre'
    ax_pol=add_ax_pol(fig1,ha=ha)

    xyavg=[0,0]
    for d in dets:
        xyavg[0]+=dets[d]['xy'][0]
        xyavg[1]+=dets[d]['xy'][1]
    xyavg[0] /= len(dets)
    xyavg[1] /= len(dets)

    ax_xy=add_ax_xy(fig1,fillfact=0.8,ha="centre",va="centre")
    ff_ins=0.25
    ax_xy_ins=add_ax_xy_inset(fig1,fillfact=ff_ins,ha="left",va="bottom",frameon=False,ctr=xyavg)

    xl_xy=ax_xy.get_xlim()
    yl_xy=ax_xy.get_ylim()
    xr=xl_xy[1]-xl_xy[0]
    yr=xl_xy[1]-yl_xy[0]
    print(xl_xy,yl_xy)
    thr=arange(0,2*pi,2*pi/360)
    ax_xy_ins.plot(xyavg[0]+0.45*sin(thr),xyavg[1]-0.45*cos(thr),lw=3,c='gray',alpha=0.5)
    ax_xy.plot(0.05*sin(thr),0.05*cos(thr),lw=3,c='gray',alpha=0.5)
    ax_xy.plot([-0.05*sin(pi/4),xl_xy[0]+0.5*ff_ins*(xr)],\
        [0.05*cos(pi/4),yl_xy[0]+0.95*ff_ins*(xr)],lw=3,c='gray',alpha=0.5)
    ax_xy.plot([0,xl_xy[0]+ff_ins*(xr)*(0.5-0.45*sin(5*pi/4))],\
        [-0.05,yl_xy[0]+ff_ins*(xr)*(0.5+0.45*cos(5*pi/4))],lw=3,c='gray',alpha=0.5)

    # add detector marks
    for d in dets:
        add_det(dets[d],ax_xy_ins,pol=False,col=cols['det'])

    # add detpair marker(s)
    for p in detpairs:
        add_det_pair(detpairs[p],ax_xy_ins,scale,col=cols['pair'],extlines=True)
        add_det_lines(detpairs[p],dets,ax_xy)

    # add_ctr(fig1,ha=ha)

    if title!=None:
        add_title(title,ax_xy,loc='upper left')
    # if showcalc:
    #     add_dt_text(detpairs,ax_xy,swang,scale,c,loc='lower left')

    # output unannotated plot
    if save:
        print(('saving to %s'%(fileOut)))
        plot.savefig(fileOut)



def makeplot_wv(dets,detpairs,waveang,scale=3e4,fignum=None,fileOut=None,noise=0,m1=5,m2=5,save=False,title=''):

    from matplotlib.ticker import AutoMinorLocator

    # plot waveforms
    wvFileIn='../Waveform/data/m1-%d-m2-%d.txt'%(m1,m2)
    data=wv.read_txt(wvFileIn)
    data=wv.add_zeros(data,1)

    if fileOut==None:
        fileOut='Temp/%ddet_waveforms.pdf'%(len(list(dets.keys())))

    tcoal=-min(data['t'])
    data['t']=data['t']-min(data['t'])
    swang=waveang*pi/180.
    c=3.e5

    # set relative time indices for waveforms
    for p in detpairs:
        pair=detpairs[p]
        # if not pair.has_key('dt'):
        # init_dets(dets)
        # init_pairs(detpairs,dets)
        # dt=do_calc(pair,swang,scale,c)
        # print 'dt',dt
        d1=dets[pair['d1']]
        d2=dets[pair['d2']]
        # if (not d1.has_key('relt')) and (not d2.has_key('d2')):
        #     # neither has relt key
        #     d1['relt']=0.
        #     d2['relt']=d1['relt'] - pair['dt']
        # elif (not d1.has_key('relt')):
        #     # only d1 is missing relt
        #     d1['relt']=d2['relt'] + pair['dt']
        # elif (not d2.has_key('relt')):
        #     # only d2 is missing relt
        #     d2['relt']=d1['relt'] - pair['dt']
        d1['relt']=0.
        d2['relt']=d1['relt'] - pair['dt']

    ndet=len(dets)

    # plot waveform
    fig1=plot.figure(fignum,figsize=(11,8))
    fig1.clf()

    for d in dets:
        fig1.add_subplot(ndet,1,dets[d]['order'])
        ax=fig1.gca()
        xscale=1.e3
        yscale_exp=-round(log10(max(data['h_re'])))
        yscale=10**yscale_exp
        ax.grid('on',which='major',c='0.7',lw=1,ls='-')
        ax.minorticks_on()
        ax.grid('on',axis='x',which='minor',c='0.7',lw=1,ls=':')

        tsc=(data['t']+dets[d]['relt'])*xscale
        # print ('Det %s: t=%.2g ms'%(d,dets[d]['relt']*1e3))
        # print (d1,d2)
        if noise>0:
            data_noise=wv.add_noise(data,noise,seed=dets[d]['order'])
            hsc=data_noise['h_re']*yscale
        else:
            hsc=data['h_re']*yscale
        ax.plot(tsc,hsc,lw=2,alpha=1,c=dets[d]['col'],
            label='Detector %s'%(d))

        ax.set_xlim((tcoal-0.02)*xscale,(tcoal+0.02)*xscale)
        ax.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax.set_ylabel('Strain ($\mathrm{x} 10^{%g}$)'%(-yscale_exp))
        if dets[d]['order']==ndet:
            ax.set_xlabel('Time (ms)')

        ylim=max(abs(array(ax.get_ylim())))
        ax.set_ylim(-ylim,ylim)

        plot.legend(loc='best',frameon=False)

        if title!='' and dets[d]['order']==1:
            # print title
            plot.title(title,fontsize='x-large')

    # output annotated plot
    if save:
        print((waveang,'saving to %s'%(fileOut)))
        plot.savefig(fileOut)

def makeplot_ex(dets,detpairs,waveang,scale=3e4,title=None,fignum=None,fileOut=None,save=False,showcalc=True):

    # draw detectors, GW and time differences

    if fileOut==None:
        if showcalc:
            fileOut='Temp/%ddet_plot_ex%g.pdf'%(len(list(dets.keys())),waveang)
        else:
            fileOut='Temp/%ddet_plot_ex%g_nocalc.pdf'%(len(list(dets.keys())),waveang)

    swang=waveang*pi/180.

    cols={'det':'b','pair':'0.7','wave':'g','dt':'r','calc':'k'}
    # speed of light in km/s
    c=3.e5

    init_dets(dets)

    # define detector-pairs
    init_pairs(detpairs,dets)

    # initialise figure
    fig1=plot.figure(fignum,figsize=(11,8))
    fig1.clf()

    # add axes
    if showcalc:
        ha='right'
    else:
        ha='centre'
    ax_pol=add_ax_pol(fig1,ha=ha)
    ax_xy=add_ax_xy(fig1,ha=ha)

    # add detector marks
    for d in dets:
        add_det(dets[d],ax_xy,pol=False,col=cols['det'])

    # add detpair marker(s)
    # for p in detpairs:
    #     add_det_pair(detpairs[p],ax_xy,scale,col=cols['pair'])

    # add sinusoidal grAv wave
    swmoff=0.2
    add_sinwave(swang,ax_xy,col=cols['wave'])

    # add deltaT for det pair(s)
    if showcalc:
        for p in detpairs:
            add_linepair(detpairs[p],ax_xy,swang,scale,c,col=detpairs[p]['col'],alen=0.07)
            # do_calc(detpairs[p])
        add_calc_text(detpairs,ax_xy,swang,scale,c,loc='lower left')
    else:
        for p in detpairs:
            add_linepair(detpairs[p],ax_xy,swang,scale,c,col=detpairs[p]['col'],alen=0.07,
                txt=r'$t_{%s}-t_{%s}$'%(detpairs[p]['d1'],detpairs[p]['d2']))

    # add_ctr(fig1,ha=ha)

    if title!=None:
        add_title(title,ax_xy,loc='upper left')

    # output annotated plot
    if save:
        print(('saving to %s'%(fileOut)))
        plot.savefig(fileOut)

    # plot.show()

    return()

def makeplot_invex(dets,detpairs,waveang,scale=3e4,title=None,fignum=None,fileOut=None,save=False):

    # draw detectors and possible angles

    if fileOut==None:
        fileOut='Temp/%ddet_plot_invex%g.pdf'%(len(list(dets.keys())),waveang)

    swang=waveang*pi/180.

    cols={'det':'b','pair':'0.7','wave':'g','dt':'r','calc':'k'}
    # speed of light in km/s
    c=3.e5

    init_dets(dets)

    # define detector-pairs
    init_pairs(detpairs,dets)

    # initialise figure
    fig1=plot.figure(fignum,figsize=(11,8))
    fig1.clf()

    # add axes
    ax_pol=add_ax_pol(fig1,ha='right')
    ax_xy=add_ax_xy(fig1,ha='right')

    # add detector marks
    for d in dets:
        add_det(dets[d],ax_xy,pol=False,col=cols['det'])

    # add detpair marker(s)
    # for p in detpairs:
    #     add_det_pair(detpairs[p],ax_xy,scale,col=cols['pair'])

    # add lines
    for p in detpairs:
        add_pair_angles(detpairs[p],ax_xy,swang)

    if title!=None:
        add_title(title,ax_xy,loc='upper left')
    add_invcalc_text(detpairs,ax_xy,swang,scale,c,loc='lower left')

    # output annotated plot
    if save:
        print(('saving to %s'%(fileOut)))
        plot.savefig(fileOut)

    # plot.show()

    return()
