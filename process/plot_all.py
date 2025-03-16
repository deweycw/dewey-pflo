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

plottting commands
'''
#print('\ncondition: '+sys.argv[1])
#pltitle = sys.argv[1]
#y_loc = [1,15,31,42,58]
#file_name = pltitle + '_transition_distance'
#outDir = pltitle


#dates = matplotlib.dates.date2num(list_of_datetimes)
#matplotlib.pyplot.plot_date(dates, values)

#s = '2018-04-01 12:00:00'
#start_date = datetime.strptime(s, "%Y-%m-%d %H:%M:%S")

os.chdir('output')
out_files = []
for ff in os.listdir(): 
	if '.csv' in ff : 
		cf = ff
		out_files.append(cf)


for cf in out_files: 
	data = pd.read_csv(cf)
	print(cf)

	fig, host = plt.subplots()
	host.plot(data['time'],data['C'])
	ydata = data['C']
	ymax = max(ydata)*1.1
	#ymax = 1.2e-2
	ymin = 0
	host.legend()
	host.set_ylim(ymin,ymax)

	#host.set_xlim(datemin,datemax)
	#host.format_xdata = mdates.DateFormatter('%Y-%m-%d')
	name = cf.split('.')[0]
	host.set_title(name) 
	host.set_ylabel(name)
	host.ticklabel_format(style='scientific', axis='y',scilimits=(-3,4))
	plt.savefig(name+'.png', dpi=300)
	plt.close(fig)
	#data = 0
	#plt.show()
	
	#pdf = matplotlib.backends.backend_pdf.PdfPages(year+'_'+mt+'/' +cc+'_trans'+ tt+'.pdf')
	#for f in range(1, plt.gcf().number+1):pdf.savefig(f) ## will open an empty extra figure :(
	#plt.close('all')
	#pdf.close()



'''

ax = data.plot(x = 'date', y = '1')
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

ax.set_ylabel("distance (m) to denitrification front")
ax.set_title(pltitle) 
plt.ticklabel_format(style='plain', axis='y', useMathText=True,scilimits=(0,0))

ax.get_legend().remove() 
x_axis = ax.axes.get_xaxis()
x_axis.set_label_text('foo')
x_label = x_axis.get_label()
##print isinstance(x_label, matplotlib.artist.Artist)
x_label.set_visible(False)

# Create a Rectangle patch
production_area = mpl.patches.Rectangle((0,0),60,1,linewidth=1,edgecolor='lightgray',facecolor='lightgray')

# Add the patch to the Axes
ax.add_patch(production_area)

plt.ylim([0, 2.5])
#plt.xlim([-0,60])
#plt.show()
plt.savefig(outDir +'/Da_DO_' + file_name + '.png', dpi=300)
'''

'''



data = pd.read_csv(outDir + '/transect_Total_Tracer [M]_1.csv')
ax = data.plot(x="distance_m", y=y_data)
#ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

ax.set_xlabel("distance (m)")
ax.set_ylabel("Tracer (M)")
ax.set_title(pltitle) 
plt.ticklabel_format(style='sci', axis='y', useMathText=True,scilimits=(0,0))
plt.ylim([0, 0.00001])
plt.xlim([-1,60])
#plt.show()
plt.savefig(outDir +'/Tracer_' + file_name + '.png', dpi=300)



nitrate_data = pd.read_csv(outDir + '/transect_Total_NO3- [M]_1.csv')
ax = nitrate_data.plot(x="distance_m", y=y_data)
#ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

ax.set_xlabel("distance (m)")
ax.set_ylabel("[NO3] (M)")
ax.set_title(pltitle) 
plt.ticklabel_format(style='sci', axis='y', useMathText=True,scilimits=(0,0))
plt.ylim([0, 0.000006])
plt.xlim([-1,30])
#plt.show()
plt.savefig(outDir +'/NO3_' + file_name + '.png', dpi=300)


data = pd.read_csv(outDir + '/transect_Total_O2(aq) [M]_1.csv')

ax = data.plot(x="distance_m", y=y_data)
#ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

ax.set_xlabel("distance (m)")
ax.set_ylabel("DO (M)")
ax.set_title(pltitle) 
plt.ticklabel_format(style='sci', axis='y', useMathText=True,scilimits=(0,0))
plt.ylim([0, 0.000245])
plt.xlim([-1,30])
#plt.show()
plt.savefig(outDir +'/DO_' + file_name + '.png')


data = pd.read_csv(outDir + '/transect_Total_NH3(aq) [M]_1.csv')

ax = data.plot(x="distance_m", y=y_data)
#ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

ax.set_xlabel("distance (m)")
ax.set_ylabel("NH3(aq) (M)")
ax.set_title(pltitle) 
plt.ticklabel_format(style='sci', axis='y', useMathText=True,scilimits=(0,0))
plt.ylim([0, 0.0000075])
plt.xlim([-1,30])
#plt.show()
plt.savefig(outDir +'/NH3_' + file_name + '.png', dpi=300)
'''

