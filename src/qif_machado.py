from brian import *
from brian.library.IF import quadratic_IF
import Graphics as artist
n = 2


taum = 20*ms
taue = 5*ms
taui = 10*ms

gmax=1.

Ee = 60.*mV
Ei = -20.*mV

Vt = 10.*mV
Vr = 0.*mV
refractory = 5*ms


tau_pre = 20*ms
tau_post = 20*ms
dA_pre = 0.01
dA_post = -0.01

eqs = Equations('''
	dV/dt = (Vr-V + ge*(Ee-V)+I)*(1./taum):volt
	dge/dt = -ge/taue : 1
	I:mV
	''') 
# NB 1: conductances are in units of the leak conductance
brain = NeuronGroup(n,model=eqs, reset=Vr,threshold=Vt,refractory=refractory,order=1)

synapses = Connection(brain,brain,'ge')
synapses[0,1] = rand()*gmax
stdp = ExponentialSTDP(synapses,tau_pre,tau_post,dA_pre,dA_post,wmax=gmax,update='mixed')

trace = StateMonitor(brain,'V',record=True)
M = StateMonitor(brain,'ge',record=True)
raster = SpikeMonitor(brain)

run(0.25*second)
brain[0].I = 11*mV
run(0.5*second)

figure()
plot(M.times,M[0]/gmax)
print M[0]
fig = figure()
ax = fig.add_subplot(111)
for i in range(n):
	ax.plot(trace.times/ms,trace[i]/mV)
artist.adjust_spines(ax)

'''
fig = figure()
ax = fig.add_subplot(111)
artist.simple_raster(raster.spiketimes,ax)
artist.adjust_spines(ax)
ax.set_xlabel('Time')
ax.set_yticks(np.arange(n)+1)
ax.set_yticklabels(range(n))
'''
show()