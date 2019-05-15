import random
import numpy as np
import matplotlib.pyplot as plt
sigma = 1
ray = []
ran = []
x = []
y = []
for i in range(0,10000):
	tmp = np.random.random()
	np.log(tmp)
	tmp2 = np.random.random()
	ray.append(np.sqrt(-2*np.log(tmp)))
	x.append(ray[i]*np.cos(tmp2*2*np.pi))
	y.append(ray[i]*np.sin(tmp2*2*np.pi))

plt.hist(ray, bins=50)
#plt.plot(x,y, '*')

plt.show()