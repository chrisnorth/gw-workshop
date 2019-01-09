# wrapper script to make interferometer plots

import interferometer as inter

save=True
wide=True

# plot labelled plots
inter.makesumplot(offset=0,amp=[1,0.5],fignum=1,title=r'Constructive Interference',save=save,wide=wide)
inter.makesumplot(offset=0.25,fignum=2,title=r'Offset: $\lambda/4$',save=save,wide=wide)
inter.makesumplot(offset=0.5,fignum=3,title=r'Destructive Interference',save=save,wide=wide)

# plot unlabelled plots
inter.makesumplot(offset=0,amp=[1,0.5],fignum=4,save=save,wide=wide)
inter.makesumplot(offset=0.25,fignum=5,save=save,wide=wide)
inter.makesumplot(offset=0.5,fignum=6,save=save,wide=wide)