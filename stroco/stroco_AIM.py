# -*- coding: utf-8 -*-
""" 


"""
__author__ = "kha"
__version__ = "0.0.1"


import numpy as np
import matplotlib.pyplot as plt





def n_d(ed, beta=1, U=1):
	""" """

	res = 2*np.exp(-beta*ed) 
	res += 2*np.exp(-beta*(2*ed+U))
	res /= (1 + 2*np.exp(-beta*ed) + np.exp(-beta*(2*ed+U)))
	return res


if __name__ == "__main__":
	ed = np.linspace(-5, 1, 10001)
	nd1 = n_d(ed, beta=10)
	nd2 = n_d(ed, beta=20)
	plt.plot(ed, nd1)
	plt.plot(ed, nd2)