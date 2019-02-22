# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:16:44 2017

@author: kh
"""

import numpy as np
import matplotlib.pyplot as plt




x = np.arange(-10,10,0.01) * 10**-2



def f(x):
    """ """
    return np.log(np.abs(x))/x