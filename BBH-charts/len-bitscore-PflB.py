import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from os.path import dirname, abspath

#~~~~~~~~~~~~~~~~~~~
# GET DATA FOR PLOT
#~~~~~~~~~~~~~~~~~~~

parentdir = dirname(dirname(abspath(__file__)))
fpath = parentdir+'/CjPflB-BBHs/CjPflB-BBH_matches.tsv'
f = pd.read_csv(fpath, sep='\t')

x = f['slen'].apply(lambda x: x/788)
y = f['bitscore']

#~~~~~~~~~~~~~~~~~~~~~~~
# MAKE THE SCATTER PLOT:
#~~~~~~~~~~~~~~~~~~~~~~~

#make the plot
fig, ax = plt.subplots()

ax.scatter(x, y, color='#cc3333', marker='.', edgecolor='k', linewidth=0.7)
#ax.scatter(x_psi, y_psi, color='C1', marker='^', edgecolor='k', linewidth=0.7, label='psiblast hits')

# axis labels
plt.xlabel('hit length/ query length ratio')
plt.ylabel('bitscore')

#add xticks, label every second xtick
#plt.xticks(np.arange(0.5, 1.1, 0.1))
#for label in ax.xaxis.get_ticklabels()[-1:]:
#	label.set_visible(False)
#plt.yticks(range(50, 110, 10))
ax.xaxis.set_minor_locator(AutoMinorLocator(n=2))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=2)) # turn on minor ticks to be for every 1 bin

#turn on y-axis gridlines
ax.yaxis.grid(b=True, which='major', linestyle='--', color='0.85')
ax.yaxis.grid(b=True, which='minor', linestyle='--', color='0.85')
plt.tick_params(axis='y', which='minor', length=0)
ax.set_axisbelow(True)

#figure size
fig.set_size_inches(6, 4)

#position legend, show plot
plt.legend()
plt.tight_layout()
plt.savefig('PflB-delta-len-bitscore-scatter', dpi=300)
plt.show()
