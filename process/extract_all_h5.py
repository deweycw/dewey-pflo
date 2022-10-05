import io 
import os 
import sys
import h5py 
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import csv
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
''' 
	

'''

def extractComponentList(h5):
	file = h5py.File(h5,'r')
	groups = list(file.keys())	

	componentList = []
	for gg in groups:
		if 'Time' in gg:
			componentList = list(file[gg].keys())	

	return componentList


def extractTimes(h5):
	file = h5py.File(h5,'r')
	groups = list(file.keys())	

	times = []
	for gg in groups:
		if 'Time' in gg: times.append(float(gg[7:18])) # gets time step 

	times = np.array([times]).T

	return times

def extractComponent_transect(h5,component,dy,return_times,rm_bc_cells):
	file = h5py.File(h5,'r')
	groups = list(file.keys())	
	dx, dz = 1, 1
	nx, nz = int(1/dx), int(1/dz) 

	componentList = []
	times = []
	data = []
	for gg in groups:
		if 'Time' in gg:
			times.append(float(gg[7:18])) # gets time step 
			dset = np.array(file[gg][component])
			if rm_bc_cells == 'y': component_transect = dset[nx,1:-1,nz] # returns an array containing value of component at all cells along model transect at given time
			else: component_transect = dset[nx-1,:,nz-1]  #			else: component_transect = dset[nx,:,nz]

			componentList.append(component_transect)

	componentList = np.asarray(componentList)

	#print(componentList)
	#print(times)
	times = np.array([times]).T
	data = np.append(times,componentList,axis = 1)
	data = data[data[:,0].argsort()]
	times = data[:,0]
	data = data[:,1:]
	
	if return_times == 'y': return times, data
	else: return data 

def extractComponent_atTime(h5,component,dy,distY,time):
	file = h5py.File(h5,'r')
	groups = list(file.keys())	
	dx, dz = 1, 1
	nx, ny, nz = int(0.5/dx), int(distY/dy), int(0.5/dz) 

	componentList = []
	data = []
	
	for gg in groups:
		if str(np.format_float_scientific(time,precision=5,unique=False)).replace('e','E') in gg:
			dset = np.array(file[gg][component])
			componentAtCoords = dset[nx,ny,nz]
	return componentAtCoords

# make list of h5 output files
h5_files = []
for ff in os.listdir(): 
	if 'pflotran-' in ff and '.h5' in ff: 
		h5f = ff
		h5_files.append(h5f)

componentList = extractComponentList(h5_files[0])

#print('\ncondition: '+sys.argv[1])
newDir = 'output' #sys.argv[1]
Path(newDir).mkdir(parents=True, exist_ok=True)


firstComponent = True
transLength = 1
dy = 1
for component in componentList: 
	Path(newDir).mkdir(parents=True, exist_ok=True)
	#os.mkdir(newDir + '/' + component)
	print(component)

	
	dy_vector = np.arange(0,transLength,dy)
	header = list(np.arange(0,transLength,dy))
	header = ['time','C']
	times = extractTimes(h5_files[0]) # only purpose of this is to get the times
	data_holder = times#[:,None]
	for file in h5_files:  # each h5 corresponds to a realization
		#print(component)
		#realization = int(file.split('.')[0].split('-')[1])
		times, component_con = extractComponent_transect(file,component,dy,return_times = 'y',rm_bc_cells = 'n')  # returns an array with number of rows equal to number of output times & number of columns equal to cells in model (ny)
		times = times[:,None]
		allData = np.hstack((times,component_con))
		

		## this chunk of code writes data with rows = to hr and columns = to position along flowpath; uncomment to use 
		fname = newDir + '/' + component+'.csv'
		with open(fname, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow(header)
			writer.writerows(allData)
		'''

		## this code writes data with rows = pos along flow path, cols = hr 
		#fname = newDir + '/transect_'+component+'_'+str(realization)+'.csv'
		fname = component+'.csv'
		pdData = pd.DataFrame(data = allData[:,1:],index = data_holder[:,0],columns = list(np.arange(0.,transLength,dy)))
		pdData = pdData.sort_index()
		pdData = pdData.transpose()

		pdData.to_csv(fname,sep=',',index_label='distance_m')
	'''	



		