import numpy as np
import matplotlib.pyplot as plt


v = np.arange(-100,100,1)


n_alpha = ((-0.01*(v+34))/ ( np.exp(-0.1*(v+34)))) * 3
n_beta = 0.125* np.exp(-(v+44.)/80)
n_inf = n_alpha/(n_alpha+n_beta)

plt.figure()
plt.plot(n_alpha)
plt.figure()
plt.plot(n_beta)
plt.figure()
plt.plot(n_inf)
plt.show()
