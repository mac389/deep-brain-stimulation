import numpy as np
import matplotlib.pyplot as plt
import Graphics as artist

from numpy import mean,cov,cumsum,dot,linalg,size,flipud,argsort
from matplotlib import rcParams

rcParams['text.usetex'] = True
data = np.loadtxt('./pain-v-spare.tsv',delimiter='\t',skiprows=1)
M = np.loadtxt('../data/connections-pain.csv',delimiter=',')

format = lambda text:r'\Large \textbf{\textsc{%s}}'%text
def princomp(A,numpc=3):
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

npc = 3
eigenvectors,projection,eigenvalues = princomp(M,numpc=npc)
area_names = ['Spinal Cord','Medulla','Medullary RF','Pons','Pontine RF','CL Nuc.','VPL','NAc','MD','Cortex']



data_projections = data.T.dot(eigenvectors)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(data_projections[99:1099,0],data_projections[99:1099,1],c='k', s=50, label=format('Background'))
ax.scatter(data_projections[1199:2199,0],data_projections[1199:2199,1],c='r',s=50,label=format('NAc'))
ax.scatter(data_projections[2299:3299,0],data_projections[2299:3299,1],c='g',s=50,label=format('Spine'))
artist.adjust_spines(ax)
ax.set_xlabel(format('Principal component 1'))
ax.set_ylabel(format('Principal component 2'))
plt.legend(frameon=False,numpoints=1)
plt.tight_layout()
plt.show()
#Print contribution of each area to eigenvectors

'''
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(np.absolute(projection.T),interpolation='nearest',aspect='auto', cmap=plt.cm.binary)
artist.adjust_spines(ax)
ax.set_xticks(range(npc))
ax.set_xlabel(format('Principal component'))
ax.set_yticks(range(projection.shape[1]))
ax.set_yticklabels(map(format,area_names))
cbar = plt.colorbar(cax)
plt.tight_layout()
plt.show()
'''