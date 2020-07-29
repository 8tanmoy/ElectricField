#--convergence test for electric field with cutoff radius--
#--generate data for one trajectory frame using field_space_rcut.py--

import pandas as pd
import numpy as np

numRcut 	= 75	#number of r points
fnames 		= [ f'./rcut_{ii}/field_nointra.dat' for ii in range(numRcut) ]
dataframes	= [ pd.read_csv(fname, header=None, sep=' ') for fname in fnames]
means		= [ df.mean(axis=0) for df in dataframes ]
stds		= [ df.std(axis=0) 	for df in dataframes ]
means		= np.array(means)
stds		= np.array(stds)
print(f'shape[mean] : {means.shape}')
print(f'shape[std] : {stds.shape}')
assert len(means) == len(stds), 'means and stds dim error'
#--write plottable datafile--

file 		= open("convergence_rcut.dat", "w")
file.write("#rcut\t#mean_field[x,y,z]\t\t\t#std_field[x,y,z]\n")
for ii in range(len(means)):
	rcut 	= 0.1 + 0.05 * ii
	mean0	= "\t".join(str(x) for x in means[ii])
	std0	= "\t".join(str(y) for y in stds[ii])
	file.write(f'{rcut}\t{mean0} \t{std0}\n')
file.close()

#import matplotlib.pyplot as plt
#
#x = np.array([0.5 + 0.05 * jj for jj in range(numRcut)])
#y = means[:,2]
#e = stds[:,2]
#
#plt.errorbar(x, y, e, linestyle='None', marker='^')
#
#plt.show()
