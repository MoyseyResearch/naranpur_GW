#!/usr/bin/python

import math
import matplotlib.gridspec as gridspec
import itertools
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters
import scipy.sparse.linalg
import scipy.interpolate

nt = 3				# number of timesteps
dt = 0.1			# size of time step [weeks]
nx = 100			# number of grids along EW axis
ny = 100			# number of grids along NS axis
dx = 135.0			# grid size along EW axis [m]
dy = 135.0			# grid size along NS axis [m]
K  = 0.864*7			# hydraulic conductivity [m/wk]
Ss = 0.002			# storativity [dimensionless]
b  = 250			# thickness [m]
CN = 63				# curve number [dimensionless]
infRate = 0.07			# infiltration rate [wk^-1]
soilThickness = 1000		# soil layer thickness [mm]
initialSoilMoisture = 20	# initial soil moisture [mm]
irrigation = 50*7		# irrigation applied [L/wk]

p = (Ss/K)*(dx*dy)/dt		# [m]

def setRainfall(week):
  week = int(np.mod(week,52))
  mean = open('rainfall.txt').read().splitlines()[week].split()[0]
  var  = open('rainfall.txt').read().splitlines()[week].split()[1]
  P = np.random.normal(mean,var) * (12.0/52.0)
  if P>0: return P
  else:   return 0

def setET(week):
  week = int(np.mod(week,52))
  mean = open('ET.txt').read().splitlines()[week].split()[0]
  var  = open('ET.txt').read().splitlines()[week].split()[1]
  ET = np.random.normal(mean,var) * (12.0/52.0)
  if ET>0: return ET
  else:    return 0

def setRu(P,CN):
  S = 1000 / CN - 10
  return (P-0.2*S)**2 / (P+0.8*S)

def setIrrigation(Q):
  Q				# L/acre/day
  Q = Q * 7.0			# L/acre/week
  Q = Q * 1e3			# mL/acre/week, cm^3/acre/week
  Q = Q / 1e6			# m^3/acre/week
  Q = Q / (4046.86)		# m/week
  Q = Q * 1e3			# mm/week
  return Q

def setStorage(S,P,I,ET,Ru,i):
  if i>0:
    storage = ( S[i-1]-(P[i]+I[i]-ET[i]-Ru[i])/infRate ) * math.exp(-infRate) + (P[i]+I[i]-ET[i]-Ru[i])/infRate
    if storage<0: storage=0
    if storage>soilThickness:
      Ru[i] += storage-soilThickness
      storage = soilThickness
  else:  storage = S[i]
  return storage

t     = np.zeros(nt)
P     = np.zeros(nt)
ET    = np.zeros(nt)
Ru    = np.zeros(nt)
I     = np.zeros(nt)
S     = np.zeros(nt)
dS    = np.zeros(nt)
Gr    = np.zeros(nt)
B     = np.zeros(nt)
GW    = np.zeros([nt+1,nx,ny])
raw   = np.ones([100,100])

file = open('elevation.csv')
for i in range(100):
  for j in range(100):
    raw[i,j] = file.readline()
elev = filters.gaussian_filter(raw,5)
elev = np.ones([nx,ny])*1500
GW[0] = elev-400

#pump = np.genfromtxt('pump_data.csv', delimiter=',')
#well_rates = np.array([50,100,200,0]);
#Qwell = np.zeros(pump.shape[0])
#for i in range(pump.shape[0]):
#  Qwell[i]= pump[i,2]*well_rates[pump[i,1]-1]
#well_elev = pump[:,0]
#for i in range(well_elev.size): well_elev[i] = elev.ravel()[pump[i,0]]

nodeN = np.ones([nx,ny])
nodeN[:,0] = 2
nodeN[0,:] = 2
nodeN[:,ny-1] = 2
nodeN[nx-1,:] = 2

A = np.zeros([nx*ny,nx*ny])

for i in range(nx):
  for j in range(ny):
    k = i*nx+j
    neighbors = np.array([k-ny, k+ny, k-1 ,k+1])
    neighbors = np.delete(neighbors,np.where(neighbors<0))
    neighbors = np.delete(neighbors,np.where(neighbors>(nx*ny-1)))
    if nodeN[i,j]==1:
      A[k,k] = -4-p
      A[k,neighbors] = 1
    else: A[k,k] = 1

plt.imshow(A)
plt.colorbar()
plt.savefig('A.eps')

print np.linalg.cond(A)

'''S[0]=initialSoilMoisture
for i in range(nt):
  t[i]  = i*dt
  P[i]  = setRainfall(i*dt)
  ET[i] = setET(i*dt)
  I[i]  = setIrrigation(irrigation)
  Ru[i] = setRu(P[i]+I[i],CN)
  S[i]  = setStorage(S,P,I,ET,Ru,i)
  Gr[i] = infRate*S[i]
  B[i]  = P[i]+I[i]-ET[i]-Ru[i]-Gr[i]
  recharge = (Gr[i]/1000)*Ss*(dx*dy/dt)

  print "Running GW (transient bicg) Model: Iteration %i, week %f..."%(i,i*dt)
  b = GW[i].ravel()*p

  for j in range(nx):
    for k in range(ny):
      if j==0 or j==nx-1 or k==0 or k==ny-1: b[j*nx+k] = GW[0,j,k]

#      else: b[j*nx+k] -= recharge
  GW[i+1] = np.dot( np.linalg.inv(A), b.T ).reshape([nx,ny])
#  GW[i+1] = scipy.sparse.linalg.bicgstab(A,b)[0].reshape([nx,ny])


for i in [0,nt-1]:
  print "Plotting hydrologic results: Iteration %i, week %f..."%(i,i*dt)
  plt.figure(num=None, figsize=(32,24), dpi=2400, facecolor='w', edgecolor='k')
  gs = gridspec.GridSpec(8,7)
  ax1 = plt.subplot(gs[0:4, 0:4],projection='3d')
  X,Y = np.meshgrid(np.linspace(0,nx*dx,nx),np.linspace(0,ny*dy,ny))
  ax1.plot_surface(X,Y,nodeN,linewidth=0.3)
  surf = ax1.plot_surface(X,Y,GW[i+1,:,:], vmin=min(GW[1:-1].ravel()), vmax=max(GW[1:-1].ravel()), rstride=1,cstride=1,cmap=mpl.cm.jet,linewidth=0.01,antialiased=False)
  plt.colorbar(surf,shrink=0.75,aspect=10)
  ax1.plot_surface(X,Y,elev,alpha=0.4,linewidth=0.3)
  ax1.view_init(azim = 180+70,elev = 35)
  plt.title('Week: %f'%((i+1)*dt), fontsize=20)
  ax2 = plt.subplot(gs[0,4:6])
  ax2.plot( t[0:i+1], P[0:i+1],  'b-', label='Precipitation')
  ax2.plot( t[0:i], ET[0:i], 'r-', label='Evapotranspiration')
  ax2.plot( t[0:i], Ru[0:i], 'g-', label='Runoff')
  ax2.plot( t[0:i], I[0:i],  'y-', label='Irrigation')
  ax2.set_ylabel('[mm/wk]',fontsize=20)
  ax2.set_xlim(0,max(t))
  ax2.set_ylim( min( min(P.ravel()), min(ET.ravel()), min(Ru.ravel()), min(I.ravel()) ),  max( max(P.ravel()), max(ET.ravel()), max(Ru.ravel()), max(I.ravel()) ) )
  ax2.legend(loc=1,prop={'size':10})
  ax3 = plt.subplot(gs[1,4:6])
  Br = B.clip(max=0)
  Bb = B.clip(min=0)
  ax3.fill_between( t[0:i], Br[0:i], 0,  color='red')
  ax3.fill_between( t[0:i], Bb[0:i], 0,  color='blue')
  ax3.plot( t[0:i], B[0:i],  'k-', label='Balance')
  ax3.set_ylabel('[mm/wk]',fontsize=20)
  ax3.set_xlim(0,max(t))
  ax3.set_ylim( min(B.ravel()),  max(B.ravel()) )
  ax3.legend(loc=1,prop={'size':10})
  ax4 = plt.subplot(gs[2,4:6])
  ax4.plot( t[0:i], S[0:i],  'k-', label='Storage' )
  ax4.set_ylabel('[mm]',fontsize=20)
  ax4.set_xlim(0,max(t))
  ax4.set_ylim( min(S.ravel()),  max(S.ravel()) )
  ax4.legend(loc=1,prop={'size':10})
  ax5 = plt.subplot(gs[3,4:6])
  ax5.plot( t[0:i], Gr[0:i],  'k-', label='Percolation' )
  ax5.set_xlabel('Time [weeks]',fontsize=20)
  ax5.set_ylabel('[mm/wk]',fontsize=20)
  ax5.set_xlim(0,max(t))
  ax5.set_ylim( min(Gr.ravel()),  max(Gr.ravel()) )
  ax5.legend(loc=1,prop={'size':10})
  plt.savefig('./frames/transient_'+'%05d'%i+'.png', bbox_inches='tight')
  plt.clf()'''
