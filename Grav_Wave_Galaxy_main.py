"""
Name: Gravity_Wave_Galaxy_main.py
Author: Mark Monaghan
Date: 30/11/2017
Description: Main program for project
"""
#Imports / Some unused, represent potential improvements
from __future__ import division #Patched division function
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import cm #colourmap
import moviepy.editor as mpy
from matplotlib.colors import LogNorm

#Constants / Not all used
pc_to_km = 3.08567758129e13
sec_per_year = 365.25 * 86400
deg_to_rad = np.pi / 180.0
rad_to_deg = 180.0 / np.pi
pi = np.pi
#Parameters
rcore = 6e4 #must be less than rgalaxy
rgalaxy = 3e5
rdist = 2 * rgalaxy
e1 = 0.75 #eccentricity bounds => Sb Type Galaxy
e2 = 1.0
N = int(4e5)
velocity = 5e-4
timestep = 1e4
angular_offset = 0.00004 #rad

#Arrays
a = np.random.normal(0.0,rgalaxy/2,N) #semi-major axis array
aT = np.linspace(rcore,rgalaxy,10) #_T implies test element
x = np.zeros(N)
y = np.zeros(N)
b = np.empty(N)
EE = np.empty(N)
orbV = np.empty(N)
theta = np.random.uniform(0.0,2.0*pi,N)
tilt = np.empty(N)
testAng = np.linspace(0.0,2.0*pi,1000)

#Functions
def E(a,rcore,rgalaxy,rdist,e1,e2,N,EE):
    for i in range(N):
        if a[i] < rcore:
            EE[i] = 1.0 + ((a[i]/rcore)*(e1 - 1.0))
        elif ((rcore < a[i]) and (a[i] <= rgalaxy)):
            EE[i] = e1 + (((a[i] - rcore)/(rgalaxy - rcore))*(e2 - e1)) #maybe brackets error
        elif (rgalaxy < a[i]) and (a[i] < rdist):
            EE[i] = e2 +(((a[i] - rgalaxy)/(rdist - rgalaxy))*(1.0 - e2))
        else:
            EE[i] = 1.0
    return EE

def orbVel(a,rcore,velocity,rgalaxy):
    for i in range(len(a)):
        if (a[i] < rcore):
            dv = ((velocity-(2.0*velocity))/rcore)
            orbV[i] = (2.0*velocity)+a[i] * dv
        elif (a[i]>=rcore):
            dv = (0.5*velocity - velocity) / (rgalaxy-rcore)
            orbV[i] = velocity+dv * (a[i] - rcore)
    return orbV

def Tilt(a,angular_offset):
    for k in range(N):
        tiltTmp = a[k] * angular_offset
        if (tiltTmp > pi):
            tiltTmp -= (2.0*pi)
        tilt[k] = tiltTmp
    return tilt

def f(a,b,theta,tilt):
    x = (a*(np.cos(theta)*np.cos(tilt)))-(b*(np.sin(theta)*np.sin(tilt)))
    y = (a*(np.cos(theta)*np.sin(tilt)))+(b*(np.sin(theta)*np.cos(tilt)))
    return x,y

#fn calls
ecc = E(a,rcore,rgalaxy,rdist,e1,e2,N,EE)
b=ecc*a
tilt = Tilt(a,angular_offset)
sx,sy = f(a,b,theta,tilt)
orbV = orbVel(a,rcore,velocity,rgalaxy)

xmin = sx.min()
xmax = sx.max()
ymin = sy.min()
ymax = sy.max()

#-----#Test Code#-----#
#one orbit
# for k in range(len(aT)):
#     xT = None
#     yT = None
#     tiltT = aT[k] * angular_offset
#     if (tiltT > pi):
#         tiltT -= (2.0*pi)
#
#     xT = (aT[k]*(np.cos(testAng)*np.cos(tiltT)))-(b[k]*(np.sin(testAng)*np.sin(tiltT)))
#     yT = (aT[k]*(np.cos(testAng)*np.sin(tiltT)))+(b[k]*(np.sin(testAng)*np.cos(tiltT)))
#
#     plt.plot(xT,yT,'k-')
#
#     plt.show()
#
#-----#Test Code#-----#


# #matplotlib animation ref: http://jakevdp.github.io/blog/2012/08/18/
#                                    matplotlib-animation-tutorial/
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(xmin,xmax), ylim=(ymin,ymax))
line, = ax.plot([], [],'k,')

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    global theta,orbV, tilt
    theta += (orbV*timestep*i)
    tilt += 0.0000001
    sx,sy = f(a,b,theta,tilt)
    line.set_data(sx, sy)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=3000, interval=2, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('DensityWaveGalaxy.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# #hexbins
# fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(50, 50))
# fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
# ax = axs[0]
# hb = ax.hexbin(sx, sy, gridsize=50, cmap='inferno')
# ax.axis([xmin, xmax, ymin, ymax])
# ax.set_title("Hexagon binning")
# cb = fig.colorbar(hb, ax=ax)
# cb.set_label('counts')
#
# ax = axs[1]
# hb = ax.hexbin(sx, sy, gridsize=50, bins='log', cmap='inferno')
# ax.axis([xmin, xmax, ymin, ymax])
# ax.set_title("With a log color scale")
# cb = fig.colorbar(hb, ax=ax)
# cb.set_label('log10(N)')

#simple pyplot
# plt.ion()
# plt.plot(sx,sy,'m*')

#plt.show()
