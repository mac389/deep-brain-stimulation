import brian_no_units
import cPickle

from brian import *
from time import time
import os

#Regular stimulation at cerebellum doesn't mean that thalamocortical synapses are firing regularly 
#Perhaps irregular stimulation two synapses away leads to regular stimulation

from optparse import OptionParser
from numpy import random

op = OptionParser()
op.add_option('--rate', dest='rate', type='int', 
      help='Average firing rate of input stimuli')
op.add_option('--rhythm',dest='rhythm',type='str',
      help='REGULAR for tonic firing, IRREGULAR for phasic firing')
op.add_option('--duration',dest='duration',type='int',
      help='duration of recording in seconds')
op.print_help()

opts,args = op.parse_args()
if len(args) > 0:
      op.error('This script only takes arguments preceded by command line options.')
      sys.exit(1)
defaultclock.dt = 1*ms
N = 1000
taum = 10 * ms
tau_pre = 20 * ms
tau_post = tau_pre
Ee = 0 * mV
vt = -54 * mV
vr = -60 * mV
El = -74 * mV
taue = 5 * ms
F = 15 * Hz
gmax = .05
dA_pre = .01
dA_post = -dA_pre * tau_pre / tau_post * 1.05
dA_post *= gmax
dA_pre *= gmax

eqs_neurons = '''
dv/dt=(ge*(Ee-vr)+El-v)/taum : volt   # the synaptic current is linearized
dge/dt=-ge/taue : 1
'''
os.system('say "starting simulation with stimulation at %d hurts and %s"'%(opts.rate,opts.rhythm))

'''
Accurate stimulation pattern
Across all frequencies
'''

def lognormal_spiketime(start,stop,avg_rate):
      increments = random.lognormal(0.2,0.2,(stop-start)*opts.rate) #assumes rate is in units of Hz
      return np.arange(start,stop,avg_rate) + increments

if opts.rhythm != 'IRREGULAR': #Default action if rhythm type not recognized is to make regular rhythm
      spiketimes = [(i,t) for i in xrange(N) for t in np.arange(0,opts.duration,1./opts.rate)]
else:
      spiketimes = [(i,t) for i in xrange(N) for t in lognormal_spiketime(0,opts.duration,1./opts.rate)]
input = SpikeGeneratorGroup(N,spiketimes)
neurons = NeuronGroup(1, model=eqs_neurons, threshold=vt, reset=vr)
S = Synapses(input, neurons,
             model='''w:1
             A_pre:1
             A_post:1''',
             pre='''ge+=w
             A_pre=A_pre*exp((lastupdate-t)/tau_pre)+dA_pre
             A_post=A_post*exp((lastupdate-t)/tau_post)
             w=clip(w+A_post,0,gmax)''',
             post='''
             A_pre=A_pre*exp((lastupdate-t)/tau_pre)
             A_post=A_post*exp((lastupdate-t)/tau_post)+dA_post
             w=clip(w+A_pre,0,gmax)''')
neurons.v = vr
S[:,:]=True
S.w='rand()*gmax'

rate = PopulationRateMonitor(neurons)
M = StateMonitor(S,'w',record=True,clock=Clock(dt=100*ms)) # monitors synapses number 0 and 1

start_time = time()
run(opts.duration * second, report='text')
cPickle.dump({'synapses':M.values,'postsynaptic_rates':rate.rate,'postsynaptic_rate_smoothed':rate.smooth_rate(100*ms), 
      'postsynaptic_times':rate.times/second},open('../results/output-%s-%s.pkl'%(str(opts.rate),opts.rhythm),'wb'))
os.system('say "weightmonitor is finished" ')