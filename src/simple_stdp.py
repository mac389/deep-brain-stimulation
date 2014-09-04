import numpy as np
from sandbox import weight_change,poisson_neurons, simple_raster
from progressbar import ProgressBar

import Graphics as artist
import matplotlib.pyplot as plt

dt = 0.001 
tau_m = 0.010
tau_w = 0.10
tau_theta = 0.050
duration = 5
rate = 60

#units and normalization are important
#remember thresholds

u = np.zeros((int(duration/dt),))
#v = np.zeros_like(u)
w = np.zeros_like(u)
theta = np.zeros_like(u)
I = np.zeros_like(u)

threshold = {'v':-4,'u':-2}

I,Itimes=poisson_neurons(rate,duration)
#Initial conditions
u[0] = 10*np.random.random()
#v[0] = 10*np.random.random()
w[0] = np.random.random()
theta[0] = np.random.random()

def sigmoid(x):
	return 1 / (1 + np.exp(-x))

pbar = ProgressBar(maxval=duration/dt)
for t in range(1,int(duration/dt)):
	#forward Euler unstable for this system, multi-step method unstable
	
	u_inf = I[t]-threshold['u']
	u[t] = u_inf + (u[0]-u_inf)*np.exp(-t/tau_m)
	v = w[t-1]*u[t]
    #Assume v at steady state, so not modeling the millisecond to millisecond dynamics
	    
	theta_inf = v**2 
	theta[t] = theta_inf + (theta[0]-theta_inf)*np.exp(-t/tau_theta)
    
	w[t] = w[t-1] + dt/tau_w*(u[t]*v*(v-theta[t]))
    
	pbar.update(t)
pbar.finish()  
fig,axs = plt.subplots(sharex=True,nrows=4,ncols=1)
x = np.arange(0,duration,dt)
simple_raster(Itimes,axs[0])
artist.adjust_spines(axs[0],[])
axs[0].set_ylabel('Input')
for ax,(var,label) in zip(axs[1:],[(u,'Presynpatic'),(w,'Weight'),(theta,'threshold')]):
    ax.plot(x,var,'k',clip_on=False)
    artist.adjust_spines(ax,['left'])
    ax.set_ylabel(label)
#ax.plot(x[:200],v[:200],'k.',linewidth=2,clip_on=False,label='Postynaptic')
#ax.plot(x,u,'k--',linewidth=2,clip_on=False,label='Presynaptic')
#ax.plot(x,w,'r.',linewidth=2,clip_on=False,label='Weight')
#ax.plot(x,theta,'r--',linewidth=2,clip_on=False,label='Threshold')
#ax.plot(x[:200],I[:200],'g',linewidth=2,clip_on=False,label='Input')
#artist.adjust_spines(ax)
#ax.set_ylabel('Firing rate (impulses per second)')
#plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('../results/simple_stpd.png',dpi=300)

'''u[t] = u[t-1]+dt/tau_m*(-u[t-1] + I[t-1]-threshold['u'])
	v[t] = v[t-1]+dt/tau_m*(-v[t-1] + w[t-1]*u[t-1]-threshold['v'])
	    
	w[t] = w[t-1]+dt/tau_w*(-w[t]+u[t]*v[t]*(v[t]-theta[t]))
	theta[t] = theta[t-1] + dt/tau_theta*(v[t]**2-theta[t-1])
'''