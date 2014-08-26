import itertools

import numpy as np
import networkx as nx
import numpy.linalg as la
import matplotlib.pyplot as plt 
import Graphics as artist

from pprint import pprint
from numpy import mean,cov,cumsum,dot,linalg,size,flipud,argsort
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['text.usetex'] = True
format = lambda txt: r'\huge \textbf{%s}'%txt
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['ytick.major.pad'] = 12 
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']

def princomp(A,numpc=0):
	# computing eigenvalues and eigenvectors of covariance matrix
	M = (A-mean(A.T,axis=1)).T # subtract the mean (along columns)
	[latent,coeff] = linalg.eig(cov(M))
	p = size(coeff,axis=1)
	idx = argsort(latent) # sorting the eigenvalues
	idx = idx[::-1]       # in ascending order
	# sorting eigenvectors according to the sorted eigenvalues
	coeff = coeff[:,idx]
	latent = latent[idx] # sorting eigenvalues
	if numpc < p and numpc >= 0:
		coeff = coeff[:,range(numpc)] # cutting some PCs if needed
		score = dot(coeff.T,M) # projection of the data in the new space
	return coeff,score,latent

M = np.loadtxt('connections.csv',delimiter=',')
G = nx.from_numpy_matrix(M,create_using=nx.DiGraph())

Lap = nx.directed_laplacian_matrix(G)
test = princomp(Lap,numpc=5)

coeff,score,latent = princomp(M,numpc=4)


nuclei = ['Cortex','Striatum','GPi/SNr','STN', 'GPe','Thalamus']
inputs = dict(zip(nuclei,np.eye(len(nuclei))))

def distance(a,b):
	a = np.array(a)
	b = np.array(b)
	return 1-abs(a.dot(b)/(linalg.norm(a)*linalg.norm(b)))

stn_pattern = [0,1,1,-1,0,0]
gpi_pattern = [0,0,-1,0,0,0]


stn_projection = np.array(score).dot(stn_pattern)
gpi_projection = np.array(score).dot(gpi_pattern)


overlaps = [distance(stn_projection[n-2:n],gpi_projection[n-2:n]) for n in [2,3,4]]

overlaps = [distance(stn_pattern,gpi_pattern)] + overlaps


fig,ax = plt.subplots()
width=0.35
ax.bar(np.arange(len(overlaps))-0.5*width,overlaps,width,color='k')
artist.adjust_spines(ax)
ax.set_xticks(range(len(overlaps)))
ax.set_xticklabels([format('Overall'),r'\huge $\mathbf{e}_1$',r'\huge $\mathbf{e}_2$',r'\huge $\mathbf{e}_3$'])
ax.set_ylabel(format('Similarity'))
plt.tick_params(axis='y',which='major',labelsize=20)
plt.tight_layout()
plt.show()

#fig,(tracts,scree) = plt.subplots(nrows=2,ncols=1,sharex=True)
#--Eigentracts
'''
fig = plt.figure()
tracts = fig.add_subplot(111)
cax = tracts.imshow(score.T[:,:4]*latent[:4]/latent.sum(),interpolation='nearest',aspect='auto', vmin=-2,vmax=2)
artist.adjust_spines(tracts)
tracts.set_xlabel(format('Eigentract'))
tracts.set_yticks(range(6))
tracts.set_xticks(range(4))
tracts.set_yticklabels(map(format,['Cortex','Striatum',r'$\textrm{GP}_\textrm{i}$/SNr','STN', r'$\textrm{GP}_\textrm{e}$',r'Thalamus']))
plt.tick_params(axis='x',which='major',labelsize=20)
cbar = plt.colorbar(cax,use_gridspec=True)
cbar.ax.tick_params(labelsize=20)
tracts.xaxis.labelpad=20
plt.tight_layout()
plt.show()
'''
#--scree plot
'''
scree.stem(range(4),latent[:4]/latent.sum(),linefmt='k-',markerfmt='ko',basefmt='k-', bottom=0,clip_on=False)
artist.adjust_spines(scree)
scree.set_ylabel(r'\Large $\lambda$',rotation='horizontal')
scree.set_xlabel(format('Eigentract'))
scree.set_xticks(range(4))
plt.tight_layout()
plt.show() 
'''

#--create different stimulation parameters

energy = lambda activity: activity.dot(M).dot(activity.T)
pattern= np.array(list(itertools.product([-1,1],repeat=6))).T

'''
fig, (patterns,stability) = plt.subplots(nrows=2,ncols=1,sharex=True)

patterns.imshow(pattern,interpolation='nearest',aspect='auto',cmap=plt.cm.binary)
artist.adjust_spines(patterns)
patterns.set_yticks(np.arange(6))
patterns.set_yticklabels(map(format,['Cortex','Striatum',r'$\textrm{GP}_\textrm{i}$/SNr','STN', r'$\textrm{GP}_\textrm{e}$',r'Thalamus']))
patterns.set_xlabel(format('Active'))

stability.plot(-0.5*np.array([energy(pat) for pat in pattern.T]),'k',linewidth=2)
artist.adjust_spines(stability)
stability.set_xlabel(format('Pattern'))
stability.set_ylabel(format('Energy'))
plt.tight_layout()
plt.show()
'''
#--project different patterns of stimulation onto eigenvalues of connection matrix
#--TODO: Color ones where STN or GPi are active, low alpha for all others

'''
eig1 = coeff[:,0]
eig2 = coeff[:,1]

x = np.array([pat.dot(eig1) for pat in pattern.T])
x /= latent[0]

y = np.array([pat.dot(eig2) for pat in pattern.T])
y /= latent[1]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x,y,c=['g' if pat[-1] == 1 else 'r' for pat in pattern])
artist.adjust_spines(ax)
ax.set_xlabel(format('PC 1'))
ax.set_ylabel(format('PC 2'))
plt.tight_layout()
plt.show()
'''
