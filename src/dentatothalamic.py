import numpy as np

from numpy import random

'''
   Cerebellum --< Thalamus --< Cortex      --< Excitation
      ---						 |        
       |_________________________|         --| Inhibition

'''

inputs = {'cortex':1,'thalamus':1}
T = 1000
tau = {'membrane':0.01,'synapses':0.1}
window = {'tau':0.01,'a':1}
H = lambda t: -window['a']*np.exp(-t/window['tau']) if t >0 else window['a']*np.exp(t/window['tau']) 

I = #input
v = np.zeros((1,T))
w = np.zeros_like(v)
x = np.zeros_like(w)
y = np.zeros_like(w)

#initial conditions

for t in range(1,T):
	x[t] = x[t-1]+ 1/tau['synapses']*(-x[t-1]+0.1*I[:t].sum())

	y[t] = y[t-1]+ 1/tau['synapses']*(-y[t-1]+0.1*I[:t].sum())
	
	w[t] = aplus(w[t-1])*x[t]*