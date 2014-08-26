import numpy as np
import Graphics as artist
import matplotlib.pyplot as plt

from numpy.random import random_sample
from matplotlib import rcParams
from numpy.random import randint
from random import random
from itertools import izip_longest
from pprint import pprint

rcParams['text.usetex'] = True


#Connection matrix

params = {'duration':{'pattern':1000,'interpattern':100},'noise_level':0.1,'pattern':['background','gpi','stn']}
total_duration = len(params['pattern'])*params['duration']['pattern'] + (len(params['pattern'])+1)*params['duration']['interpattern']
format = lambda txt: r'\Large \textbf{%s}'%txt
def sgn(activity):
	activity[activity>=0]=1
	activity[activity!=1]=-1
	return activity 

M = np.loadtxt('connections.csv',delimiter=',')
n_brain_areas = M.shape[1]

#Quick meaningful access to parts of connection matrix  
stn = 3
gpi = 2

area = {0:'Cortex',1:'Striatum',2:'Gpi',3:'STN',4:'Gpe',5:'Thalamus'}
inverse_area = {v: k for k, v in area.items()}

current = {'stn':np.array([0,0,0,-2,0,0]),
		   'gpi':np.array([0,0,-2,0,0,0]),
		   'background':np.zeros((n_brain_areas,))}

current['interpattern_interval'] = current['background']
rectify = lambda current: current*np.greater(current,0)

#--Create patterns of stimulation
pattern_ = filter(None,[val for pair in izip_longest(['interpattern_interval']*(len(params['pattern'])+1),
										params['pattern'],fillvalue=None) for val in pair])

stimulation_pattern = np.concatenate([np.tile(current[condition],(params['duration']['interpattern' if 'interpattern' in condition else 'pattern'],1)).T 
						for condition in pattern_],axis=1)

afferent_activation_pattern = np.concatenate([2*np.tile(rectify(M[eval(condition),:]) if 'gpi' in condition  or 'stn' in condition else np.zeros((n_brain_areas,)),
								(params['duration']['interpattern' if 'interpattern' in condition else 'pattern'],1)).T 
									for condition in pattern_],axis=1)

#--Initial conditions
idx = randint(low=0,high=n_brain_areas,size=(total_duration,))
background_noise = 2*randint(2,size=(total_duration))-1

#initial_activity = 2*randint(2,size=(n_brain_areas,))-1
initial_activity = -np.ones((6,))
#Simulate  Parkinson's patient by making GPi,STN active
initial_activity[inverse_area['Gpi']] = 1
initial_activity[inverse_area['STN']] = 1

v = np.tile(initial_activity,(total_duration,1)).T
#v = np.zeros((n_brain_areas,total_duration))
v[:,0] = initial_activity

threshold = random_sample(size=(total_duration,))
F = lambda current: 1./(1+np.exp(-current))

debugging=False
for t in range(1,total_duration):
		# Determine the state of unit idx[t] at time t
		I = afferent_activation_pattern[:,t]#*v[:,t] # Afferent activation 
		I[idx[t]] += M[:,idx[t]].dot(v[:,t-1])+stimulation_pattern[idx[t],t] #Target inhibition

		#Stochastically update the state of unit idx[t]
		v[idx[t],t] = 1 if F(I[idx[t]]) > threshold[t] else -1
		v[:,t] += stimulation_pattern[:,t] + afferent_activation_pattern[:,t]
		v[:,t] = sgn(v[:,t])	
		#---Debugging
		if debugging:
			print '------------'
			print 'Unit %d selected, corresponding to %s'%(idx[t],area[idx[t]])
			print 'Kernel: ',M[:,idx[t]]
			print 'Network activity: ', v[:,t-1]
			print 'Stimulation pattern: ', stimulation_pattern[idx[t],t]
			print 'Total input %.02f'%I[5]
			print 'Probability: %.02f, threshold:%.02f'%(F(I[idx[t]]),threshold[t])
			print 'Acitvity %d -> %d'%(v[idx[t],t-1],v[idx[t],t])
			print '--------'
		#---
	
artist.raster_plot(v,duration=[params['duration']['interpattern' if 'interpattern' in condition else 'pattern'] for condition in pattern_])
#	parasite_labels=[condition if 'interpattern' not in condition else '' for condition in pattern_])


artist.firing_rate_plot(v,duration=[params['duration']['interpattern' if 'interpattern' in condition else 'pattern'] for condition in pattern_],
	parasite_labels=[condition if 'interpattern' not in condition else '' for condition in pattern_])