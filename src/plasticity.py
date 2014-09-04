from brian import *
from time import time

from matplotlib import rcParams
import Graphics as artist

rcParams['text.usetex'] = True


'''
   Cerebellum --< Thalamus --< Cortex      --< Excitation
      ---						 |        
       |_________________________|         --| Inhibition

'''

N = 1000
taum = 10 * ms
tau_pre = 20 * ms
tau_post = tau_pre
Ee = 0 * mV
vt = -54 * mV
vr = -60 * mV
El = -74 * mV
taue = 5 * ms
gmax = 0.01
F = 15 * Hz
dA_pre = .01
dA_post = -dA_pre * tau_pre / tau_post * 2.5

H = lambda t: dA_pre*exp(-t/tau_pre) if t >0 else dA_post*exp(t/tau_post) 

eqs_neurons = '''
dv/dt=(ge*(Ee-vr)+El-v)/taum : volt   # the synaptic current is linearized
dge/dt=-ge/taue : 1
'''

input = PoissonGroup(N, rates=F)
neurons = NeuronGroup(1, model=eqs_neurons, threshold=vt, reset=vr)
synapses = Connection(input, neurons, 'ge', weight=rand(len(input), len(neurons)) * gmax,
                    structure='dense')
neurons.v = vr

stdp = ExponentialSTDP(synapses, tau_pre, tau_post, dA_pre, dA_post, wmax=gmax, update='mixed')

rate = PopulationRateMonitor(neurons)
weights = StateMonitor(neurons,'ge',record=True)
start_time = time()
run(10 * second, report='text')
print "Simulation time:", time() - start_time


summary = figure()
rate_plot = summary.add_subplot(211)
rate_plot.plot(rate.times / second, rate.smooth_rate(100 * ms),'k',linewidth=2)
artist.adjust_spines(rate_plot)
rate_plot.set_xlabel(r'\Large Time $\left(s\right)$')
rate_plot.set_ylabel(r'\Large Impulses per second')
weight_dist = summary.add_subplot(212)
weight_dist.hist(synapses.W.todense(), 20,color='k')
artist.adjust_pinse(weight_dist)
tight_layout()

fig2 = figure()
weight_course = fig2.add_subplot(111)
plot(weights,'k')



#visualize STDP curve
fig = figure()
ax = fig.add_subplot(111)
ax.plot(np.arange(-40,40),[H(x*ms)/gmax for x in np.arange(-.40,.40,.01)],'k',linewidth=2)
artist.adjust_spines(ax)
ax.set_ylabel(r'\Large $\frac{\Delta w_{ij}}{w_{ij}}$',rotation='horizontal')
ax.set_xlabel(r'\Large $t_{post}-t_{pre} \; \left(\textrm{ms}\right)$ ')
show()

