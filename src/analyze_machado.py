import cPickle,os,json

import numpy as np
import matplotlib.pyplot as plt
import Graphics as artist 

from scipy.optimize import curve_fit
from matplotlib import rcParams, mlab
from optparse import OptionParser

op = OptionParser()
op.add_option('--rate', dest='rate', type='int', 
      help='Average firing rate of input stimuli')
op.add_option('--rhythm',dest='rhythm',type='str',
      help='REGULAR for tonic firing, IRREGULAR for phasic firing')
op.print_help()

opts,args = op.parse_args()
if len(args) > 0:
      op.error('This script only takes arguments preceded by command line options.')
      sys.exit(1)

rcParams['text.usetex'] = True
audio = False

READ = 'rb'
WRITE = 'wb'
if audio:
	os.system('say "Beginning analysis by opening %s"'%'output dot pickle')
data = cPickle.load(open('../results/output.pkl',READ))

format = lambda text: r'\Large \textbf{\textsc{%s}}'%text

'''
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data['postsynaptic_times'],data['postsynaptic_rate_smoothed'],'k',linewidth=2)
artist.adjust_spines(ax)
ax.set_ylabel(format('Firing rate (impulses per second)'))
ax.set_xlabel(format('Time (s)'))
plt.tight_layout()
'''

wfig = plt.figure()
wax = wfig.add_subplot(111)
wax.hist(data['synapses'][0,:],bins=20, color='b',alpha=0.8, histtype='stepfilled',normed=True,label=format('Before'))
plt.hold(True)
wax.hist(data['synapses'][-1,:],bins=20, color='r',alpha=0.5, histtype='stepfilled',normed=True,label=format('After'))

artist.adjust_spines(wax)
wax.set_xlabel(format('Synaptic strength'))
wax.set_ylabel(format('Prevalence'))
plt.legend(frameon=False,loc='upper left')
plt.tight_layout()
plt.savefig('../results/weight_shift-%d-%s.png'%(opts.rate,opts.rhythm),dpi=300)
print [data['synapses'][0,:].mean(),data['synapses'][-1,:].mean()]
print [data['synapses'][0,:].var(),data['synapses'][-1,:].var()]
json.dump({'mu':[data['synapses'][0,:].mean(),data['synapses'][-1,:].mean()],'var':[data['synapses'][0,:].var(),data['synapses'][-1,:].var()]},
	open('../results/output-weights-%d-%s.json'%(opts.rate,opts.rhythm),WRITE))