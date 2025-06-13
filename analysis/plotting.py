import numpy as np
import matplotlib.pyplot as plt
import matplotlib

rcParams = matplotlib.rcParams
rcParams["axes.prop_cycle"] = plt.matplotlib.cycler(color=["#FAF0C9","#F4D35E", "#f0c62d", "#EEB902", "#EE964B", "#ea7e1f", "#b0413e", "#72002B", "#4F000B", "#457b9d", "#1D4A5D", "#2D3142", "#061A40", "#040303"])
rcParams["axes.labelsize"]= 14
rcParams["axes.xmargin"]= 0
rcParams["axes.ymargin"]= .1
rcParams["lines.markersize"]= 6
rcParams["lines.linewidth"]= 2
rcParams["figure.autolayout"]= True
rcParams["figure.facecolor"]= "white"
rcParams["font.size"]= 14
rcParams["grid.color"]= "0"
rcParams["grid.linestyle"]= "-"
rcParams["legend.edgecolor"]= "k"
rcParams["legend.fontsize"]= 12
rcParams["xtick.labelsize"]= 12
rcParams["ytick.labelsize"]= 12
rcParams["xtick.direction"]= "in"
rcParams["ytick.direction"]= "in"
# Use LaTeX for text
rcParams["text.usetex"]= True

dir = "../output/example/"

density_file = dir + "density_potential/density_potential_0.dat"
iter_file = dir + "iterations/iterations_0.dat"

######### Density Potential Plotting #########
data = np.loadtxt(density_file)
x = data[:, 0]
rho_plus = data[:, 1]
rho_minus = data[:, 2]
psi = data[:, 3]
Y = data[:, 4]

plt.subplots(2,1, figsize=(5,7), dpi=500, sharex=True)
plt.subplot(2, 1, 1)
plt.plot(x, rho_plus, label=r'$\rho_+$', color='blue')
plt.plot(x, rho_minus, label=r'$\rho_-$', color='red')
plt.ylabel(r'$\rho_i\sigma^3$')
plt.xlim(0,10)
plt.ylim(0,1)
plt.legend(loc='upper center',fancybox=False)

plt.subplot(2, 1, 2)
plt.plot(x, psi, label=r'$\psi$', color='black')
plt.plot(x, Y, label=r'$Y$', color='green')
plt.xlabel(r'$x$')
plt.ylabel(r'$\psi$, $Y$')
plt.xlim(0,10)
plt.legend(loc='upper center',fancybox=False)
plt.savefig('profiles.png', bbox_inches='tight', dpi=500)    