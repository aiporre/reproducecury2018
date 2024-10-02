import numpy as np
from diffeocentroids import blobimage, match_cpp
import matplotlib.pyplot as plt

# Define parameters
Ns = np.array([30, 30, 3])

# Create grids
h1s = np.linspace(0, 1, Ns[0] + 1)
h2s = np.linspace(0, 1, Ns[1] + 1)
h3s = np.linspace(0, 1, Ns[2] + 1)
hxs, hys, hzs = np.meshgrid(h1s, h2s, h3s, indexing='ij')

# Permute and reshape
hxs = np.transpose(hxs, (2, 1, 0))
hys = np.transpose(hys, (2, 1, 0))
hzs = np.transpose(hzs, (2, 1, 0))
xspec = np.vstack([hxs.ravel(), hys.ravel(), hzs.ravel()])

# Create structure s
s = {}
s['xspec'] = xspec

N = np.array([30, 30, 3])
h1 = np.linspace(0, 1, N[0] + 1)
h2 = np.linspace(0, 1, N[1] + 1)
h3 = np.linspace(0, 1, N[2] + 1)
hx, hy, hz = np.meshgrid(h1, h2, h3, indexing='ij')

# Permute and reshape
hx = np.transpose(hx, (2, 1, 0))
hy = np.transpose(hy, (2, 1, 0))
hz = np.transpose(hz, (2, 1, 0))
x = np.vstack([hx.ravel(), hy.ravel(), hz.ravel()])

s['x'] = x

# Set additional parameters
s['useDef'] = 'LargeDef'
s['typefloat'] = 'double'
s['CppKer'] = {'Type': 'CauchyGpu'}
s['T'] = 10
s['sigmaV'] = 0.2
s['gammaR'] = 0
s['optim_useoptim'] = 'adaptdesc'
s['optim_maxiter'] = 100
s['optim_stepsize'] = 0.1

# Configure target
t = {}
t['method'] = 'l2image'
t['rx'] = np.arange(1, xspec.shape[1] + 1)
t['weight'] = 1
t['imsource'] = blobimage(Ns, 0.2, [0.4, 0.5, 0.5])
t['imtarget'] = blobimage(Ns, 0.2, [0.6, 0.5, 0.5])

t['basetarget'] = np.array([0, 0, 0])
t['voxsizetarget'] = np.array([1, 1, 1]) / Ns

target = [t]
s = matchCpp(s, target)

# Visualization
plt.figure(1)
plt.clf()
plt.hold(True)

a = s['phix'] - s['x']
plt.quiver(s['x'][0, :], s['x'][1, :], s['x'][2, :], a[0, :], a[1, :], a[2, :], color='r')
plt.show()