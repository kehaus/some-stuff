"""
RHEED pattern simulator


"""
__author__ = "kha"
__version__ = "0.0.1"


import numpy as np
import time

# ==============================================================
# classes
# ==============================================================
class CoordinateGrid(object):
	""" """
	def __init__(self):
		""" """

		self.shape = [10,5]
		self._grid = np.zeros(self.shape)
		self.spacing = 0.1		# [m]; micrometer spacing of detector pixel

	# def calc_coord_grid(self):
	# 	""" """

	# 	shp = self.grid.shape
	# 	for i in range(shp[0]):
	# 		for j in range(shp[1]):
	# 			self.grid[i,j] = np.array([i,j]) * self.spacing
	# 	return


class Screen(CoordinateGrid):
	""" """
	def __init__(self):
		""" """
		self.shape = [20,20]
		self.grid = np.zeros(self.shape)
		self.spacing = 10**-6		# [m]; micrometer spacing of detector pixel
		self.e1 = np.array([0,1,0]) * self.spacing		# grid unit vector
		self.e2 = np.array([0,0,1]) * self.spacing		# grid unit vector


class Plane(CoordinateGrid):
	""" """
	def __init__(self):
		""" """

		self.shape = [100,100]
		self.grid = np.zeros(self.shape)
		self.spacing = 10**-9		# [m]; lattice spacing is in nm range
		self.e1 = np.array([1,0,0]) * self.spacing		# grid unit vector
		self.e2 = np.array([0,1,0]) * self.spacing		# grid unit vector


# ==============================================================
# parameters
# ==============================================================
dist 	= 10
theta 	= 5 * np.pi/180.	# [rad]
phi		= 0 * np.pi/180		# [rad]

z = x_dist * np.tan(theta)


source 	= np.array([
	-1.*dist,
	0,
	z
])
plane 	= np.array([		# center of plane is origin
	0,
	0,
	0
])
detector 	= np.array([
	dist,
	0,
	z
])

# s = source; p = plane; d = detector
r_sp = plane - source
r_pd = detector - plane

wl = 1					# [m]; wavelength in Angstrom range
k_vec = 2*np.pi/wl 				# [rad/m]; wave vector 

# ==============================================================
# calculations
# ==============================================================

def f(r):
	"""norm of a vector"""
	return np.linalg.norm(r)

def f_(r):
	"""first derivatie of the norm of a vector"""
	return -r/f(r)

def Tf(r, r0):
	"""First order Taylor expansion"""
	return f(r0) + np.dot(f_(r0), r-r0)



scrn = Screen()
pln = Plane()


### calc intensity

# one detector pixel

def calc_phi():
	"""phase differences of scatter centers"""
	phi = np.zeros(pln.shape, dtype=np.float64)
	shp = pln.shape
	for i in range(shp[0]):
		for k in range(shp[1]):
			r = i*pln.e1 + k*pln.e2
			phi[i,k] = k_vec * Tf(r, r_sp)
	return phi

def calc_pixel_intensity(phi, pxl_i=0, pxl_j=0):
	"""intensity at one detetor pixel from all scatter centers"""
	intens = np.float64(0)
	shp = pln.shape
	for i in range(shp[0]):
		for k in range(shp[1]):
			r = pxl_i*scrn.e1 + pxl_j*scrn.e2
			r -= i*pln.e1 + k*pln.e2 
			intens += np.exp(1j* (k_vec * Tf(r, r_pd) + phi[i,k]))
	return intens

def calc_detector_intensity(phi):
	"""detector intensity map"""
	intens = np.zeros(scrn.shape, dtype=np.complex64)
	shp = scrn.shape
	for i in range(shp[0]):
		for k in range(shp[1]):
			t1 = time.time()
			intens[i,k] = calc_pixel_intensity(phi, i, k)
			print('    time: ', time.time()-t1)
	return intens





## 

t0 = time.time()

phi = calc_phi()
print('start')
intens = calc_detector_intensity(phi)
print('time: ', time.time()-t0)






# ==============================================================
# main
# ==============================================================
if __name__ == "__main__":

	pass