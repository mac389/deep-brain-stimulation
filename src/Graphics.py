import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1 import host_subplot

from pprint import pprint
from matplotlib import rcParams
from matplotlib.mlab import psd
from mpl_toolkits.mplot3d import Axes3D

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['ytick.major.pad'] = 12 
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']

box = dict(facecolor=None, pad=15, alpha=0)

def angle_plot(one,two=None):
	if not two:
		two=one

	#must take the product of the columns
	angles = np.array([np.inner(first,second)/(np.inner(first,first)*np.inner(second,second)) 
			for first,second in zip(one.transpose(),two.transpose())])
	print angles

def adjust_spines(ax,spines=['bottom','left']):
	''' Taken from http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html '''
	for loc, spine in ax.spines.iteritems():
		if loc in spines:
			spine.set_position(('outward',10))
			spine.set_smart_bounds(True) #Doesn't work for log log plots
			spine.set_linewidth(1)
		else:
			spine.set_color('none') 
	if 'left' in spines:
		ax.yaxis.set_ticks_position('left')
	else:
		ax.yaxis.set_ticks([])

	if 'bottom' in spines:
		ax.xaxis.set_ticks_position('bottom')
	else:
		ax.xaxis.set_ticks([])

def power_spectrum(data,Fs=20000, savename=None,show=True, cutoff=50):
	p = Periodogram(data,sampling=Fs)
	p.run()
	p.plot()
	'''
	#stop = np.where(freqs>cutoff)[0][0]
	#print stop
	fig = plt.figure()
	ax = fig.add_subplot(111)
	spec, = ax.plot(freqs,db,'o-')
	adjust_spines(ax,['bottom','left'])
	ax.set_xlabel(r'frequency $\left(Hz\right)$')
	ax.set_ylabel(r'Power $\left(dB\right)$')
	'''
	if show:
		plt.show()
	if savename:
		plt.savefig(savename,dpi=72)

def scree_plot(eigVals,cutoff=0.95,savename=None, show=False,save=True,savebase=None):
	#Assume the list is all of the eigenvalues
	rel = np.cumsum(eigVals)/eigVals.sum()
	x = np.arange(len(rel))+1
	print eigVals.shape
	fig = plt.figure()
	ax = fig.add_subplot(111)
	line, = ax.plot(x,rel)
	line.set_clip_on(False)
	adjust_spines(ax,['bottom','left'])
	ax.set_xlabel(r'$\LARGE \lambda$')
	ax.set_ylabel('Fraction of variance')
	ax.set_xlim(0,len(eigVals))
	
	cutoff_idx = np.where(rel>cutoff)[0][0]
	
	ax.axvline(x=cutoff_idx, color='r',linestyle='--', linewidth=2)
	ax.axhline(y=rel[cutoff_idx],color='r',linestyle='--',linewidth=2)
	ax.tick_params(direction='in')
	ax.annotate(r" {\Large $\mathbf{\lambda=%d}$}" % cutoff_idx,xy=(.25, .9), xycoords='axes fraction', 
											horizontalalignment='center', verticalalignment='center')
	plt.tight_layout()
	if save:
		print savebase
		plt.savefig(savebase+'_scree.png',dpi=100)
			
	if show:
		plt.show()
	plt.close()

def spike_validation(data,clusters,spiketimes=None,eiglist=None,nclus=None,savebase='res',waveforms=None,multi=False, show=False
					,save=True, adj=False):
	best = clusters['models'][np.argmax(clusters['silhouettes'])]
	nclus = best.n_clusters if not nclus else nclus
	fig = plt.figure()
	plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=.97)
	#Clusters of waveforms projected onto the first two principal components
	ax = fig.add_subplot(2,2,1, projection='3d')
	ax.set_axis_bgcolor('white')
	colors = ['#4EACC5', '#FF9C34', '#4E9A06']
	labels_ = best.labels_
	centers = best.cluster_centers_
	unique_labels = np.unique(labels_)
	for n,col in zip(range(nclus),colors):
		my_members = labels_ == n 
		cluster_center = centers[n]
		ax.scatter(data[0,my_members],data[1,my_members],data[2,my_members],c=col, s=6)
		plt.hold(True)
		ax.scatter(cluster_center[0],cluster_center[1],cluster_center[2],c=col,s=8)
	#adjust_spines(ax,['bottom','left'])
	ax.set_ylabel(r'\Large \textbf{PC2}')
	ax.set_xlabel(r'\Large \textbf{PC1}')
	ax.set_zlabel(r'\Large \textbf{PC3}')
	
	ax.set_xticklabels('')
	ax.set_yticklabels('')
	ax.set_zticklabels('')
	
	plt.tight_layout()
	if waveforms is not None:
		print 'drawing wfs'
		wfs = fig.add_subplot(2,2,3)#axes([0.37, 0.65, 0.1, 0.15])
		wfs.set_axis_bgcolor('none')
		artists = []
		for n,col in zip(range(nclus),colors):
			#my_members = labels_[:-300]== n
			my_members = labels_ == n
			print len(my_members)
			print waveforms.shape
			line, = wfs.plot(np.average(waveforms[my_members,:],axis=0),col,linewidth=2)
			line.set_clip_on(False)	
		adjust_spines(wfs,['bottom','left'])
		
		if not adj:
			wfs.set_yticks([0,100])
			wfs.set_yticklabels([r'$0$', r'$100 \; \mu V$'],rotation='vertical')
			wfs.set_xticks([0,100])
			wfs.set_xticklabels([r'$0$',r'$800 \; \mu s$'])
			wfs.spines['bottom'].set_bounds(0,100)
		else:
			wfs.set_yticks([-1000,0,1000])
			wfs.set_yticklabels([r'$-100 \; \mu V$',r'$0$', r'$100 \; \mu V$'],rotation='vertical')
			wfs.set_xticks([0,16,32])
			wfs.set_xticklabels([r'$0$',r'$400 \; \mu s$',r'$800 \; \mu s$'])
			wfs.spines['bottom'].set_bounds(0,32)
		
	sils = fig.add_subplot(2,2,2)
	sils.set_axis_bgcolor('none')
	markerline, stemlines,baseline =sils.stem(np.arange(len(clusters['silhouettes'])),clusters['silhouettes'])
	sils.tick_params(direction='in')
	#sils.axhline(y=0.5,color='r',linestyle='--',linewidth=2)
	adjust_spines(sils,['bottom','left'])
	sils.set_xticks(np.arange(len(clusters['silhouettes']))+1)
	sils.set_yticks([-1,0,1])
	sils.set_ylabel('Silhouette coefficient')
	sils.set_xlabel('Number of clusters')
	sils.set_xlim((0.5,len(clusters['silhouettes'])))
	
	xmx=100
	if spiketimes is not None:
		#break of up the spiketime vector based on clustering
		short_isi = fig.add_axes([0.8, 0.26, 0.15, 0.20])
		isi = fig.add_subplot(2,2,4)
		for n,col in zip(range(nclus),colors):
			#my_members = labels_[:-300]== n #Always add 3000 noise spikes
			my_members = labels_ == n
			these_isis = 0.1*np.diff(spiketimes[my_members])
			these_isis = these_isis[these_isis>1]
			if these_isis.size:
				
				_,_,patches=isi.hist(these_isis, histtype='stepfilled', range=(0,1000),
					alpha=0.5, bins=50)
				adjust_spines(isi,['bottom','left'])
				plt.setp(patches,'facecolor',col)

				_,_,spatches=short_isi.hist(these_isis,range=(0,100), histtype='stepfilled')
				plt.setp(spatches,'facecolor',col)
		isi.tick_params(direction='in')
		isi.set_axis_bgcolor('none')
		isi.set_ylabel(r'\# of spikes')
		isi.set_xlabel(r'Interspike interval $(ms)$')
		isi.set_xlim(xmax=xmx)
		
		
		short_isi.set_axis_bgcolor('none')
		adjust_spines(short_isi,['bottom','left'])
		short_isi.tick_params(direction='in')
		short_isi.set_ylabel(r'\# of Spikes')
		#short_isi.set_yticks(np.arange(8))
		short_isi.axvline(x=2,c='r',linewidth=2)
		short_isi.set_xlabel(r'ISI $(ms)$')
		short_isi.set_xticklabels(np.arange(0,12)[::2])

		
	if eiglist is not None and eiglist.size:
		eigfxns = fig.add_subplot(2,2,3)
		eigfxns.set_axis_bgcolor('none')
		eigfxns.tick_params(direction='in')
		#Assume 6 eigenfunctions
		nfxns =6
		span = len(eiglist[0,:])/2
		print span
		x = arange(2*span) if multi else np.arange(-span,span)
		for i in range(nfxns):
			eigfxns.plot(x,i+eiglist[i,:],'b',linewidth=2)
			plt.hold(True)
		adjust_spines(eigfxns,['bottom','left'])
		if multi:
			eigfxns.set_xlabel(r' $\left(\mu sec\right)$')
		else:
			eigfxns.set_xlabel(r'Time from spike peak $\left(\mu sec\right)$')
			eigfxns.set_xticklabels([r'\textbf{%d}'%(32*(i-5)) for i in range(10)])
			eigfxns.set_yticklabels([' '] + [r' $e_{%d}$' %i for i in range(1,nfxns+1) ])
		eigfxns.set_ylabel(r'Eigenfunctions')
		#draw_sizebar(eigfxns)

	plt.tight_layout()
	plt.savefig(savebase+'_validation.png', bbox_inches='tight')
	if show:
		plt.show()

def voltage_trace(unfiltered=None,filtered=None,threshold = 0, roi=30000,spread=10000,save=None, 
					show=False, fs = 20000, downsampling= 10,savebase=None):
					
	fig = plt.figure()
	trace_panel = fig.add_subplot(211,axisbg='none')
	start = roi-spread
	stop = roi+spread

	traces, = trace_panel.plot(unfiltered[start:stop][::downsampling],'b') #Downsample just for display
		
	spike_panel = fig.add_subplot(212,axisbg='none',sharex=trace_panel)
	spikes, = spike_panel.plot(filtered[start:stop][::downsampling],'b')
		
	panels = [trace_panel,spike_panel]
	
	for panel in panels:
		adjust_spines(panel,['bottom','left'])
	
	trace_panel.set_xlabel(r'time $\left(s\right)$')
	trace_panel.set_ylabel(r'voltage $ \left(\mu V \right)$')
	trace_panel.set_xticklabels(np.arange(start/fs,1.5+stop/fs,0.5).astype(str))

	spike_panel.set_xlabel(r'time $\left(s\right)$')
	spike_panel.set_xticklabels(np.arange(start/fs,1.5+stop/fs,0.5).astype(str))
	spike_panel.set_ylabel(r'voltage $\left(\mu V \right)$')
	#Draw threshold
	spike_panel.axhline(y=threshold,linewidth=1,color='r',linestyle='--')
	spike_panel.axhline(y=-threshold,linewidth=1,color='r',linestyle='--')
	
	plt.tight_layout()
	if save:
		print savebase
		plt.savefig(savebase+'_voltage.png',dpi=100)
			
	if show:
		plt.show()

	def ccf():	
		print 'Calculated'
		rowL=len(filenames)
		colL=rowL
		
		acf_panel,ax=subplots(rowL,colL, sharex=True, sharey=True) 
		#Should use absolute not relative normalization
		#Currently use absolute motivation
		for j in range(rowL):
			for i in range(colL):
				line, = ax[i,j].plot(arange(-w,w),ccfs[i+j], linewidth=2)
				line.set_clip_on(False)
				ax[i,j].axvline(x=0,color='r',linestyle='--', linewidth=2)
				postdoc.adjust_spines(ax[i,j],['bottom','left'])
				ax[i,j].spines['left'].set_smart_bounds(True)
				ax[i,j].spines['bottom'].set_smart_bounds(True)
				ax[i,j].set_ylabel('Covariance')
				ax[i,j].set_xlabel(r'Time $\left(ms\right)$')
				ax[i,j].set_axis_bgcolor('none')
				ax[i,j].tick_params(direction='in')
				ax[i,j].locator_params(nbins=(60/w))
				ax[i,j].annotate(r" {\Large $\mathbf{%s,%s}$}" %(tech.get_channel_id(filenames[i]),tech.get_channel_id(filenames[j])), 
								 xy=(.2, .8), xycoords='axes fraction',horizontalalignment='center', verticalalignment='center')
		tight_layout()
		savefig('test_ccf.png')

format = lambda text: r'\huge \textbf{\textsc{%s}}'%text if '$' not in text else r'\Large %s'%text

def biplot(data):
	import rpy2.robjects as robjects
	


def raster_plot(data, duration=200, parasite_labels=['','STN','GPi'],filename='test.png',
	break_pattern_ratio=0.1,axes_labels=None):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	length = float(len([label for label in parasite_labels if 'interpattern' not in label])+1)
	cax = ax.imshow(data,interpolation='nearest',aspect='auto',cmap=plt.cm.binary)
	for i,condition in enumerate(parasite_labels):
		ax.annotate(format(condition), xy=(1.05/length*i+(1+break_pattern_ratio)/length, 0.97),
	            xycoords='figure fraction',
	            horizontalalignment='left', verticalalignment='top')
	adjust_spines(ax)
	ax.yaxis.grid(True)
	ax.set_yticks(range(len(axes_labels) if axes_labels else 6))
	ax.set_yticklabels(map(format,axes_labels if axes_labels else ['Cortex','Striatum','GPi','STN', 'GPe','Thalamus']))
	ax.set_xlabel(format('Time (arbitrary units)'))

	for demarcation in np.cumsum(duration):
		ax.axvline(x=demarcation,color='r',linestyle='--', linewidth=2)
	ax.set_xlim(xmin=0)
	ax.xaxis.labelpad=20
	plt.tight_layout()
	plt.subplots_adjust(top=0.9)
	plt.savefig(filename if '.' in filename else '%s.png'%filename,dpi=400)

def firing_rate_plot(data, duration=200, parasite_labels=['Background','STN','GPi'],
	filename='test_rate.png',break_pattern_ratio=0.1):
	fig,axs = plt.subplots(ncols=1,nrows=data.shape[0],sharex=True)
	offset = 1.5
	length = float(len([label for label in parasite_labels if 'interpattern' not in label])+1)
	rates =[1/50.*np.sum(np.reshape(0.5*(1+row),(-1,50)),axis=1) for row in data]

	for_bar_graphs = np.array([map(np.average,np.array_split(0.5*(1+row),3)) for row in data])
	pprint(for_bar_graphs)

	for ax,rate,label in zip(axs,rates,map(lambda text: r'\huge \textbf{%s}'%text,
		['Cortex','Striatum','GPi/SNr','STN', 'GPe','Thalamus'][::-1])):
		ax.plot(rate,'k')
		adjust_spines(ax)
		for i,condition in enumerate(parasite_labels):
			ax.annotate(format(condition.capitalize()), xy=(1/length*i+(1+break_pattern_ratio)/length, 0.97),
		            xycoords='figure fraction',
		            horizontalalignment='left', verticalalignment='top')
		ax.set_ylabel(label,rotation='horizontal',bbox=box)
		ax.set_yticks([])
		for demarcation in np.cumsum(duration)/50.:
			ax.axvline(x=demarcation,color='r',linestyle='--', linewidth=2)
		ax.set_xlim(xmin=0)
	axs[-1].set_xlabel(format('Time (arbitrary units)'))
	plt.tight_layout()
	plt.subplots_adjust(top=0.9)
	plt.savefig(filename if '.' in filename else '%s.png'%filename,dpi=400)