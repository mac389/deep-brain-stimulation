import numpy as np
import Graphics as artist
import matplotlib.pyplot as plt

from numpy.random import random_sample
from matplotlib import rcParams
from numpy.random import randint
from random import random

rcParams['text.usetex'] = True


M = np.array([[0, 1, 0, 1, 0, 1], 
			  [0, 0,-1, 0,-1, 0],
			  [0, 0, 0,-1, 0,-1],
			  [0, 0, 1, 0, 0, 0],
			  [0, 0, 0,-1, 0,-1],
			  [1, 0, 0, 0, 0, 0]])


duration = 10000

v = np.zeros((6,3*duration))

v[:,0] = random_sample(size=(6,))
K = M-np.eye(6)
noise_level = 0.1

format = lambda txt: r'\Large \textbf{%s}'%txt
I = 0.01*np.ones_like(v[:,0].shape)
def sgn(activity):
	activity[activity>=0]=1
	activity[activity!=1]=-1
	return activity 

#Need currents to be biphasic
current = {'STN':np.array([0,0,1,-1,0,0]),'GPi':np.array([0,0,-1,0,0,-1])}

idx = randint(low=0,high=6,size=(v.shape[1],))
for t in range(1,duration+1): 
	I = M[:,idx[t]].dot(v[:,t-1])
	v[idx[t],t] = 1 if I>0 or random() < noise_level else 0

idx = randint(low=0,high=6,size=(v.shape[1],))
for t in range(duration,2*duration+1):
	I = M[:,idx[t]].dot(v[:,t-1]+current['STN'])
	v[idx[t],t] = 1 if I>0 or random() < noise_level else 0


idx = randint(low=0,high=6,size=(v.shape[1],))
for t in range(2*duration,3*duration): 
	I = M[:,idx[t]].dot(v[:,t-1]+current['GPi'])
	v[idx[t],t] = 1 if I>0 or random() < noise_level else 0

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(v,interpolation='nearest',aspect='auto',cmap=plt.cm.binary)

ax.set_yticks(range(6))
ax.set_yticklabels(map(format,['Cortex','Striatum',r'$\textrm{GP}_\textrm{i}$/SNr','STN', r'$\textrm{GP}_\textrm{e}$',r'Thalamus']))
ax.set_xlabel(format('Time (arbitrary units)'))
ax.axvline(x=duration,color='r',linestyle='--', linewidth=2)
ax.axvline(x=2*duration,color='r',linestyle='--',linewidth=2)
plt.colorbar(cax)
plt.tight_layout()
plt.show()