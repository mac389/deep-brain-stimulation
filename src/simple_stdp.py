import numpy as np
from sandbox import weight_change,poisson_neurons
from progressbar import ProgressBar

import Graphics as artist
import matplotlib.pyplot as plt

dt = 0.001 
tau_m = 0.010
tau_w = 0.1
duration = 5
rate = 30
threshold ={'u':-2,'v':-2}

#units and normalization are important
#remember thresholds

u = np.zeros((int(duration/dt),))
v = np.zeros_like(u)
w = np.zeros_like(u)
theta = np.zeros_like(u)
I = np.zeros_like(u)

stimulation_rate = 30 #impulses per second
stimulation_delay = 0
stimulation_period = 1./stimulation_rate
I,_=poisson_neurons(rate,duration)
#Initial conditions
u[0] = 10*np.random.random()
v[0] = 10*np.random.random()
w[0] = np.random.random()
theta[0] = np.random.random()

pbar = ProgressBar(maxval=duration/dt)
for t in range(1,int(duration/dt)):
	v[t] = v[t-1]+dt/tau_m*(-v[t-1] + w[t-1]*u[t-1]-threshold['v'])
	u[t] = u[t-1]+dt/tau_m*(-u[t-1] + I[t-1]-threshold['u'])
    
		pbar.update(t)
pbar.finish()  
fig = plt.figure()
ax = fig.add_subplot(111)
x = np.arange(0,duration,dt)
ax.plot(x[:200],v[:200],'k.',linewidth=2,clip_on=False,label='Presynaptic')
ax.plot(x[:200],u[:200],'k--',linewidth=2,clip_on=False,label='Postsynaptic')
ax.plot(x[:200],w[:200],'r.',clip_on=False,label='Synaptic strength')
ax.plot(x[:200],I[:200],'r--',clip_on=False,label='Input')
artist.adjust_spines(ax)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Firing rate (impulses per second)')
ax.set_ylim(ymin=0,ymax=10)
plt.legend(frameon=False)
plt.savefig('../results/simple_stpd.png',dpi=300)