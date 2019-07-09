import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from os.path import dirname, abspath
from matplotlib.ticker import MultipleLocator, MaxNLocator, AutoMinorLocator
import matplotlib.ticker as ticker

query = 'CjPflA'
parent = dirname(dirname(abspath(__file__)))
table = parent + '/' + query + "-BBHs/" + query + '-BBH_matches.tsv'
head = 'qaccver saccver bitscore evalue qlen slen length qstart qend sstart send qcovs pident'
header = head.split(" ")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET MAX NR OF ROWS IN A FILE:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df = pd.read_csv(table, sep='\t') # open bbh-matches
x = len(df) # number of rows/ entries in bbh-matches
ls = list(df.evalue) #makes a list of the evalues in bbh-matches 
# print('list length: ', len(ls))
# print('this many zeroes: ', ls.count(0)) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVERT DATA TO LOG VALUES:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

zero = [i for i in ls if i == 0] # list of all the zeroes 
data = [i for i in ls if i != 0] # list of all non-zero evalues

logs = np.floor(np.log10(data)) # convert all the numbers to log values
logs = [int(z) for z in logs]

max_point = max(logs)
min_point = min(logs)

bins = list(range(min_point, 1, 1))

def bar_data(data, x_bins=bins): #put in data[], bins, data[1] bins, etc 
	ls = []
	for ele in x_bins:
		nr = data.count(ele)
		ls.append(nr)
	return ls

#~~~~~~~~~~~~~~~~~~~~~~
# MAKE THE BAR CHART:
#~~~~~~~~~~~~~~~~~~~~~~

#set width of bars
width = .9

#set position of bar on x axis
r1 = np.arange(len(bins)) # how many bins there are 
r1_ticks = list(range(len(bins))) # the bins that i want to have as major ticks

major_bins = bins # the labels for the major ticks

#make the plot
fig, ax = plt.subplots()
plt.bar(r1, bar_data(logs), width=width) # input the data. number of bins and the corresponding height of bar for each bin

#add xticks, label every second xtick
plt.xticks(r1_ticks, major_bins)
#plt.minorticks_on()
ax.yaxis.set_minor_locator(AutoMinorLocator(n=2))
plt.yticks(list(range(0, 23, 2)))
plt.tick_params(axis='y', which='minor', length=3)
#plt.tick_params(axis='x', which='minor', length=3)
plt.tick_params(axis='x', which='major', length=6)

# axis labels
plt.xlabel(r'$\log_{10}$(e-value)')
plt.ylabel('number of BBHs')

#turn on y-axis gridlines
ax.yaxis.grid(b=True, which='major')
ax.yaxis.grid(b=True, which='minor', color='0.9')
ax.set_axisbelow(True) # place the gridlines behind the bars 

#figure size
fig.set_size_inches(8, 5)

#position legend, show plot
plt.margins(x=.005)
plt.tight_layout()
plt.savefig('PflA-delta-BBH-logevalue-bar', dpi=300)
plt.show()

