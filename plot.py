
import numpy as np
import matplotlib
#matplotlib.use('PDF')       
matplotlib.use('agg')       
import matplotlib.pyplot as plt

########################################################################

tau_T = np.loadtxt('tau_T.dat')

matplotlib.rcParams.update({'font.size':18, 'figure.autolayout': True}) 
fig, ax = plt.subplots()
ax.tick_params(direction='in', top=True, right=True)
plt.cla()

ax.plot(tau_T[:,0], tau_T[:,1],   '-',  color='black',  linewidth=1) 

ax.set_xlabel(r"$\tau$ (s)", fontsize=18)
ax.set_ylabel("T (K)",       fontsize=18)
ax.set_xscale('log')
ax.set_xlim([1E-5, 1E10])
ax.set_ylim([0, 2400])
ax.set_title(r'CH$_4$/air, 1 atm, T$_{in}$ 300 K')
ax.grid();

#plt.savefig("tau_T")
plt.savefig("tau_T.png")

