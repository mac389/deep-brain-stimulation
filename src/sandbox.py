 import numpy as np

import matplotlib.pyplot as plt
import Graphics as artist

from matplotlib import mlab
from numpy import random
dt = 0.001
H = lambda t: np.exp(-t) if t > 0 else -np.exp(-t) 

#restrict kernel to surrounding 50 ms 
def weight_change(pre,post,dt,t,kernel=H,window=50): #Assuming kernel is a lambda function
	return sum([dt*(kernel(i*dt)*post[t]*pre[min(0,t-i)]+kernel(-t*dt)*post[min(0,t-i)]*pre[t]) for i in np.arange(min(0,t-window),min(t+window,len(post)),1.)])

def tonic_stimulus(pulse_duration,interpulse_interval,stimulus_duration):
	stimulus = np.zeros((stimulus_duration,))
	pulse = np.ones((pulse_duration,))
	silence = np.zeros((interpulse_interval,))
	for i in range(int(stimulus_duration/(pulse_duration+interpulse_interval))):
		stimulus[i*(pulse_duration+interpulse_interval):(i+1)*(pulse_duration+interpulse_interval)] = np.r_[pulse,silence]
	return stimulus

def phasic_stimulus(interpulse_interval_mu,interpulse_interval_sigma,interburst_interval_mu,interburst_interval_sigma,stimulus_duration):
    #Build underlying gamma spike train
	interpulse_intervals = random.gamma(interpulse_interval_mu,interpulse_interval_sigma,stimulus_duration-1)
	interpulse_intervals/=1000.#Ignoring refractory period
	interburst_intervals = random.normal(interburst_interval_mu,interburst_interval_sigma,stimulus_duration-1)
	interburst_intervals/=1000.
    #Make sure to mention spike train analysis, ISI, Kreuz, 
    #Assume bins at ms
	
	spiketimes = np.cumsum(np.sort(np.r_[interburst_intervals,[x for x in interpulse_intervals if x > 0.01]]))
	print spiketimes
	return spiketimes

def poisson_neurons(rate,stimulus_duration,dt=0.001):
	stimulus = random.uniform(size=(int(stimulus_duration/dt),))
	stimulus[stimulus<(rate*dt)]=1
	stimulus[stimulus!=1]=0
	return (stimulus,stimulus.nonzero()[0]/1000.)

def simple_raster(spikes,ax):
	ax.plot([spikes], [np.zeros((len(spikes)))],'|k', ms=20)
	ax.set_xlim([np.min(spikes), np.max(spikes)])

'''
rate = 30
duration = 1
#spiketimes = phasic_stimulus(5,5,900,71,rate*duration)
spiketimes = poisson_neurons(30,3)
fig = plt.figure()
ax = fig.add_subplot(111)
simple_raster(spiketimes,ax)
artist.adjust_spines(ax,['bottom'])
plt.savefig('../results/sample_raster.png',dpi=300)
'''