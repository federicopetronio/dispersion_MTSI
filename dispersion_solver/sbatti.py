import matplotlib.pyplot as plt

# x axis
ks = range(5, 14)

# y axis
results = [0.77627718442, 0.779758651376, 0.779831336605, 0.780332508531, 0.780986909681, 0.803343901758, 0.804763267426, 0.799750946218, 0.801493061947]

# use LaTeX fonts in the plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# plot
plt.plot(ks, results)
# set labels (LaTeX can be used)
plt.title(r'\textbf{Mutual Information Feature Selection}', fontsize=11)
plt.xlabel(r'\textbf{Best K features}', fontsize=11)
plt.ylabel(r'\textbf{AUC Score on split11 Dataset}', fontsize=11)
plt.show()
