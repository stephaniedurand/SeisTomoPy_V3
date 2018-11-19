import sys
import os, glob
import os.path
import imp
import math
from obspy.core import AttribDict
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from obspy.taup import TauPyModel
from obspy import taup
import numpy as np
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, DictFormatter
from obspy.geodetics.base import gps2dist_azimuth
from glob import iglob
import shutil
import time
from obspy.geodetics import locations2degrees
from scipy.spatial import distance
import re
# from geographiclib.geodesic import Geodesic
from os.path import expanduser
import mpl_toolkits.axisartist.floating_axes as floating_axes

DIR = os.path.dirname(os.path.abspath(__file__))
home = expanduser("~")
DIR2 = home + '/SeisTomoPy_files'
cwd = os.getcwd()

################################################################
#                     CLEAN DIRECTORIES
################################################################

path = DIR2 + "/output_files_cross"
dirs = os.listdir(path)
os.chdir(DIR2 + '/output_files_cross')
for i in range(len(dirs)):
	os.remove(dirs[i])

path = DIR2 + "/output_files_map"
dirs = os.listdir(path)
os.chdir(DIR2 + '/output_files_map')
for i in range(len(dirs)):
	os.remove(dirs[i])

path = DIR2 + "/output_files_corr"
dirs = os.listdir(path)
os.chdir(DIR2 + '/output_files_corr')
for i in range(len(dirs)):
	os.remove(dirs[i])

path = DIR2 + "/output_files_spectre"
dirs = os.listdir(path)
os.chdir(DIR2 + '/output_files_spectre')
for i in range(len(dirs)):
	os.remove(dirs[i])

path = DIR2 + "/output_files_time"
dirs = os.listdir(path)
os.chdir(DIR2 + '/output_files_time')
for i in range(len(dirs)):
	os.remove(dirs[i])

path = DIR2 + "/input_files"
dirs = os.listdir(path)
os.chdir(DIR2 + '/input_files')
for i in range(len(dirs)):
	os.remove(dirs[i])

os.chdir(cwd)

################################################################
#                      TOOLS FOR CROSS SECTIONS
################################################################

def midpoint(lat1,lon1,lat2,lon2):
  
  # compute the mid point between two coordinates in degrees
  assert -90 <= lat1 <= 90
  assert -90 <= lat2 <= 90
  assert -180 <= lon1 <= 180
  assert -180 <= lon2 <= 180
  
  lat1, lon1, lat2, lon2 = map(math.radians, (lat1, lon1, lat2, lon2))

  dlon = lon2 - lon1
  dx = math.cos(lat2) * math.cos(dlon)
  dy = math.cos(lat2) * math.sin(dlon)
  lat3 = math.atan2(math.sin(lat1) + math.sin(lat2), math.sqrt((math.cos(lat1) + dx) * (math.cos(lat1) + dx) + dy * dy))
  lon3 = lon1 + math.atan2(dy, math.cos(lat1) + dx)
  
  return(math.degrees(lat3), math.degrees(lon3))

def cross_section_GUI(model,para,elat,elon,slat,slon,depth,NSmax):

	lat, lon = midpoint(elat,elon,slat,slon)
	width, azi2, baz2 = gps2dist_azimuth(int(round(elat)),int(round(elon)),int(round(slat)),int(round(slon)))
	width /=  111194.9
	width2, azi, baz = gps2dist_azimuth(int(round(lat)),int(round(lon)),int(round(slat)),int(round(slon)))

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 0
	elif para =='VP':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)
	# os.chdir('../../')

	os.chdir(DIR + '/fortran_files/')
	os.system('csh cross_180.csh ' + str(int(round(lat))) + ' ' + str(int(round(lon))) + ' ' + str(int(round(azi))) +\
        ' ' + str(depth) + ' ' + str(int(round(width))) + ' ' + str(model_choice)  + ' ' + str(NSmax))

	Z = np.loadtxt(DIR2 + '/output_files_cross/profile_' + str(model) + '_' + str(para)  + '.res')
	# if (model_choice == 4) or (model_choice == 3):
	r = np.arange(3481, 6376, 5)
	# else:
	# 	r = np.arange(3481, 6331, 50)
	th = Z[0:-1:len(r),0]

	os.chdir(cwd)

	return Z, th, r

def cross_section(model,para,elat,elon,slat,slon,depth,NSmax):

	if (elon)>180:
		elon = elon-360
	if elon<-180:
		elon=elon+360
	if (slon)>180:
		slon = slon-360
	if slon<-180:
		slon=slon+360

	lat, lon = midpoint(elat,elon,slat,slon)
	width, azi2, baz2 = gps2dist_azimuth(int(round(elat)),int(round(elon)),int(round(slat)),int(round(slon)))
	width /=  111194.9
	width2, azi, baz = gps2dist_azimuth(int(round(lat)),int(round(lon)),int(round(slat)),int(round(slon)))

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 0
	elif para =='VP':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)
	# os.chdir('../../')

	os.chdir(DIR + '/fortran_files/')
	os.system('csh cross_180.csh ' + str(int(round(lat))) + ' ' + str(int(round(lon))) + ' ' + str(int(round(azi))) +\
        ' ' + str(depth) + ' ' + str(int(round(width))) + ' ' + str(model_choice)  + ' ' + str(NSmax))

	Z = np.loadtxt(DIR2 + '/output_files_cross/profile_' + str(model) + '_' + str(para)  + '.res')
	# if (model_choice == 4) or (model_choice == 3):
	r = np.arange(3481, 6376, 5)
	# else:
	# 	r = np.arange(3481, 6331, 50)
	th = Z[0:-1:len(r),0]

	os.chdir(cwd)

	return Z, th, r

def cross_section_midpoint(model,para,lat,lon,azi,depth,width,NSmax):

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 0
	elif para =='VP':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	os.chdir(DIR + '/fortran_files/')
	os.system('csh cross_180.csh ' + str(int(round(lat))) + ' ' + str(int(round(lon))) + ' ' + str(int(round(azi))) +\
        ' ' + str(depth) + ' ' + str(int(round(width))) + ' ' + str(model_choice)  + ' ' + str(NSmax))
	os.chdir(DIR)

	Z = np.loadtxt(DIR2 + '/output_files_cross/profile_' + str(model) + '_' + str(para)  + '.res')
	# th = np.arange(90.-(float(width)/2.), 90.+(float(width)/2.)+1, 1) #*np.pi/180. # deg
	# if (model_choice == 4) or (model_choice == 3):
	
	r = np.arange(3481, 6376, 5)
	# else:
	# 	r = np.arange(3481, 6331, 50)
	th = Z[0:-1:len(r),0]

	return Z, th, r

def cross_section_plot(model,para,elat,elon,slat,slon,depth,NSmax,Vmax,earthquakes=True,hotspots=True):

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)

	Z, th, r = cross_section(model,para,elat,elon,slat,slon,depth,NSmax)
	width = abs(th[-1]-th[0])

	if para == 'VS':
		para1 = 0
	elif para =='VP':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	model2 = Z[:,2]
	model_new = np.reshape(model2, (len(th),len(r)))
	model_new2=np.transpose(model_new)
	model_new2 = np.fliplr(model_new2)
	X, Y = np.meshgrid(th, r)
	fig = plt.figure()

	th0, th1 = (90-(width/2), 90+(width/2))
	r0, r1 = (6371-depth, 6371)
	thstep, rstep = (10, 500)

	# scale degrees to radians:
	tr_scale = Affine2D().scale(np.pi/180., 1.)
	tr = tr_scale + PolarAxes.PolarTransform()
	theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
	theta_tick_formatter = angle_helper.FormatterDMS()
	r_ticks = [(5871.,r"$500$"),(5371.,r"$1000$"),(4871.,r"$1500$"),(4371,r"$2000$"),(3891,r"$2500$")]
	r_grid_locator  = FixedLocator([v for v, s in r_ticks])
	r_tick_formatter = DictFormatter(dict(r_ticks))
	grid_helper = GridHelperCurveLinear(tr,
	                                      extremes=(th0, th1, r0, r1),
	                                      grid_locator1=theta_grid_locator,
	                                      grid_locator2=r_grid_locator,
	                                      tick_formatter1=None,
	                                      tick_formatter2=r_tick_formatter)

	fig.subplots_adjust(wspace=0.3, left=0.05, right=0.95, top=0.84)
	fig.subplots_adjust()

	ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
	fig.add_subplot(ax)

	# adjust x axis (theta):
	ax.axis["bottom"].set_visible(True)
	ax.axis["bottom"].toggle(all=False)
	ax.axis["top"].set_visible(True)
	ax.axis["top"].toggle(all=False)
	# adjust y axis (r):
	ax.axis["left"].set_visible(True)
	ax.axis["left"].set_axis_direction("bottom") # tick direction
	ax.axis["left"].major_ticklabels.set_axis_direction("right")

	# create a parasite axes whose transData is theta, r:
	auxa = ax.get_aux_axes(tr)
	# make aux_ax to have a clip path as in a?:
	auxa.patch = ax.patch 
	# this has a side effect that the patch is drawn twice, and possibly over some other
	# artists. So, we decrease the zorder a bit to prevent this:
	ax.patch.zorder = -2
	__s = auxa.pcolormesh(X, Y, model_new2, cmap='RdYlBu', vmin=-Vmax, vmax=Vmax)
	EQ1=np.loadtxt(DIR2 + '/output_files_cross/eq_profile.xy')
	PC3=np.loadtxt(DIR2 + '/output_files_cross/hotspots_profile.xy')
	PC=np.loadtxt(DIR + '/hotspots.xy')
	PB=np.loadtxt(DIR + '/plate_boundaries.xy')

	if earthquakes:
	  if np.size(EQ1)/2 != 1:
	    auxa.scatter(EQ1[:,0],EQ1[:,1], s=30, zorder=3,\
	                                          color="white", marker="o",\
	                                          edgecolor="k",clip_on=False)

	if hotspots:
	  if np.size(PC3)/2 != 1:
	    auxa.scatter(PC3[:,0],PC3[:,1], s=30, zorder=3,\
	                                          color="magenta", marker="o",\
	                                          edgecolor="k",clip_on=False)

	PL=np.loadtxt(DIR2 + '/output_files_cross/dots_profile.xy')
	auxa.scatter(float(PL[0,0]),6281., s=120,zorder=3,\
	                                          color="red", marker="o",\
	                                          edgecolor="k",clip_on=False)
	auxa.scatter(float(PL[1,0]),6281., s=100, zorder=3,\
	                                  color="yellow", marker="o",\
	                                  edgecolor="k",clip_on=False)
	auxa.scatter(float(PL[2,0]),6281., s=120, zorder=3,\
	                                  color="lime", marker="o",\
	                                  edgecolor="k",clip_on=False)

	fig.subplots_adjust(bottom=0.2, right=0.9, left=0.1, top=0.9)
	cb = fig.add_axes([0.2, 0.07, 0.6, 0.05])
	step2 = float(Vmax)/4
	fig.colorbar(__s,orientation="horizontal", cax = cb, ticks=np.arange(-Vmax,Vmax+step2,step2))
	auxa.add_artist(plt.Circle([0, 0], radius=5961., ls='--', lw=1, color='k', fill=False,
	                           transform=ax.transData._b, zorder=2))
	auxa.add_artist(plt.Circle([0, 0], radius=5701., ls='--', lw=1, color='k', fill=False,
	                           transform=ax.transData._b, zorder=2))

	ax2 = fig.add_axes([0.64, 0.67, .35, .4])
	map = Basemap(projection='moll', lon_0=0, resolution="c",ax=ax2)

	map.drawmapboundary(fill_color='#cccccc')
	map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)

	# fig.patch.set_alpha(0.0)

	xPC, yPC = map(PC[:,0],PC[:,1])
	xPB, yPB = map(PB[:,0],PB[:,1])
	map.scatter(xPC, yPC, s=5, zorder=10,color="magenta", marker="o",edgecolor="magenta")
	map.scatter(xPB, yPB, s=2, zorder=10,color="0.3", marker=".")

	EQ=np.loadtxt(DIR + '/catalogue.xy')
	xEQ, yEQ = map(EQ[:,0],EQ[:,1])
	map.scatter(xEQ, yEQ, s=5, zorder=10,color="white", marker="o",edgecolor="k")

	x2, y2 = map(elon,elat)
	x3, y3 = map(slon,slat)
	map.scatter(x3, y3, s=80, zorder=10,color="red", marker="o",edgecolor="k")
	map.scatter(x2, y2, s=80, zorder=10,color="lime", marker="o",edgecolor="k")

	if (elon)>180:
		elon = elon-360
	if elon<-180:
		elon=elon+360
	if (slon)>180:
		slon = slon-360
	if slon<-180:
		slon=slon+360

	clat, clon = midpoint(elat,elon,slat,slon)
	map.drawgreatcircle(clon,clat,slon,slat,ls='None',marker='o',markersize=2,color='k')
	map.drawgreatcircle(clon,clat,elon,elat,ls='None',marker='o',markersize=2,color='k')
	x1, y1 = map(clon,clat)
	map.scatter(x1, y1, s=80, zorder=10,color="yellow", marker="o", edgecolor="k")

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)

	plt.show()

################################################################
#                      TOOLS FOR MAP
################################################################

def tomomap_plot(model,para,depth,NSmax,Vmax,lon0=140,plate_boundaries=True,hotspots=True):

    path = DIR2 + "/output_files_map"
    dirs = os.listdir(path)
    os.chdir(DIR2 + '/output_files_map')
    for i in range(len(dirs)):
    	os.remove(dirs[i])
    os.chdir(DIR)
#    os.chdir('../../')

    Z2, lat, lon = tomomap(model,para,depth,NSmax)

    fig = plt.figure()
    model2 = Z2[:,2]
    X2, Y2 = np.meshgrid(lon, lat)
    model_new3 = np.reshape(model2, (len(lon),len(lat)))
    model_new4 = np.transpose(model_new3)

    ax  = fig.add_axes([0.1, 0.01, .8, .8])
    map = Basemap(projection='moll', lon_0=lon0, resolution="c",ax=ax)
    __s2 = map.pcolormesh(X2, Y2, model_new4, cmap='RdYlBu',latlon=True, vmin=-Vmax, vmax=Vmax)
    map.drawcoastlines(linewidth=1)
    step = float(Vmax)/4

    PC=np.loadtxt(DIR + '/hotspots.xy')
    PB=np.loadtxt(DIR + '/plate_boundaries.xy')
    xPC, yPC = map(PC[:,0],PC[:,1])
    xPB, yPB = map(PB[:,0],PB[:,1])
    if hotspots:
    	map.scatter(xPC, yPC, s=20, zorder=10,color="magenta", marker="o",edgecolor="k")
    if plate_boundaries:
    	map.scatter(xPB, yPB, s=2, zorder=10,color="gray", marker=".")

    path = DIR2 + "/output_files_map"
    dirs = os.listdir(path)
    os.chdir(DIR2 + '/output_files_map')
    for i in range(len(dirs)):
        os.remove(dirs[i])
    os.chdir(cwd)

    fig.colorbar(__s2,orientation="horizontal", ticks=np.arange(-Vmax,Vmax+step,step))
    plt.show()

def tomomap(model,para,depth,NSmax):

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 0
	elif para =='VP':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth) + ' ' + str(model_choice) +' ' + str(NSmax))
	os.chdir(cwd)

	filename = DIR2 + '/output_files_map/map_NEW_' + str(model) + '_' + str(para) + '.xyz'

	Z = np.loadtxt(filename)
	lat = np.arange(-89,90, 1) #*np.pi/180. # deg
	lon = np.arange(0, 360, 1)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)

	return Z, lat, lon

def tomomap_GUI(model,para,depth,NSmax):

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 0
	elif para =='VP':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth) + ' ' + str(model_choice) +' ' + str(NSmax))
	os.chdir(cwd)

	filename = DIR2 + '/output_files_map/map_NEW_' + str(model) + '_' + str(para) + '.xyz'

	Z = np.loadtxt(filename)
	lat = np.arange(-89,90, 1) #*np.pi/180. # deg
	lon = np.arange(0, 360, 1)

	os.chdir(DIR)

	return Z, lat, lon

################################################################
#                      TOOLS FOR PATH
################################################################

def path_plot(model,para,Vmax,elat,elon,slat,slon,EVT,STA,liste,model1d):

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(DIR)
	# os.chdir('../../')

	Z, th, r = cross_section(model,para,elat,elon,slat,slon,2890,40)
	width = abs(th[-1]-th[0])
	depth = 2890

	if para == 'VP':
		para1 = 0
	elif para =='VS':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	model = Z[:,2]

	model_new = np.reshape(model, (len(th),len(r)))
	model_new2=np.transpose(model_new)
	model_new2 = np.fliplr(model_new2)
	X, Y = np.meshgrid(th, r)
	fig = plt.figure()

	th0, th1 = (90-(width/2), 90+(width/2))
	r0, r1 = (0, 6281)
	thstep, rstep = (10, 500)

	# scale degrees to radians:
	tr_scale = Affine2D().scale(np.pi/180., 1.)
	tr = tr_scale + PolarAxes.PolarTransform()
	theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
	theta_tick_formatter = angle_helper.FormatterDMS()
	r_ticks = [(5871.,r"$500$"),(5371.,r"$1000$"),(4871.,r"$1500$"),(4371,r"$2000$"),(3891,r"$2500$"),(3480,r"$CMB$"),(1221,r"$ICB$"),(6371,r"$0$")]
	r_grid_locator  = FixedLocator([v for v, s in r_ticks])
	r_tick_formatter = DictFormatter(dict(r_ticks))
	grid_helper = GridHelperCurveLinear(tr,
	                                      extremes=(th0, th1, r0, r1),
	                                      grid_locator1=theta_grid_locator,
	                                      grid_locator2=r_grid_locator,
	                                      tick_formatter1=None,
	                                      tick_formatter2=r_tick_formatter)

	fig.subplots_adjust(bottom=0.3, left=0.05, right=0.95, top=0.84)
	fig.subplots_adjust()

	ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
	fig.add_subplot(ax)

	# adjust x axis (theta):
	ax.axis["bottom"].set_visible(True)
	ax.axis["bottom"].toggle(all=False)
	ax.axis["top"].set_visible(True)
	ax.axis["top"].toggle(all=False)
	# adjust y axis (r):
	ax.axis["left"].set_visible(True)
	ax.axis["left"].set_axis_direction("bottom") # tick direction
	ax.axis["left"].major_ticklabels.set_axis_direction("right")

	# create a parasite axes whose transData is theta, r:
	auxa = ax.get_aux_axes(tr)
	# make aux_ax to have a clip path as in a?:
	auxa.patch = ax.patch 
	# this has a side effect that the patch is drawn twice, and possibly over some other
	# artists. So, we decrease the zorder a bit to prevent this:
	ax.patch.zorder = -2
	__s = auxa.pcolormesh(X, Y, model_new2, cmap='RdYlBu', vmin=-Vmax, vmax=Vmax)


	# fig.subplots_adjust(bottom=0.2, right=0.9, left=0.1, top=0.9)
	cb = fig.add_axes([0.2, 0.1, 0.6, 0.05])
	step2 = float(Vmax)/4
	fig.colorbar(__s,orientation="horizontal", cax = cb, ticks=np.arange(-Vmax,Vmax+step2,step2))
	auxa.add_artist(plt.Circle([0, 0], radius=1221, ls='--', lw=1, color='k', fill=False,
                            transform=ax.transData._b, zorder=2))
	auxa.add_artist(plt.Circle([0, 0], radius=3480., ls='--', lw=1, color='k', fill=False,
                            transform=ax.transData._b, zorder=2))

	if STA.ndim == 1:
		auxa.scatter(-(width-180)/2.+width-(STA[0]),6371.-(STA[1]), s=80,zorder=3,\
                                           color="black", marker="^",\
                                           edgecolor="k",clip_on=False)
	else:
		auxa.scatter(-(width-180)/2.+width-(STA[:,0]),6371.-(STA[:,1]), s=80,zorder=3,\
                                           color="black", marker="^",\
                                           edgecolor="k",clip_on=False)
	if EVT.ndim ==1:
		auxa.scatter(-(width-180)/2.+width-(EVT[0]),6371.-(EVT[1]), s=80,zorder=3,\
                                           color="red", marker="*",\
                                           edgecolor="k",clip_on=False)
	else:
		auxa.scatter(-(width-180)/2.+width-(EVT[:,0]),6371.-(EVT[:,1]), s=80,zorder=3,\
                                           color="red", marker="*",\
                                           edgecolor="k",clip_on=False) 	

	if liste:
		liste_ph = liste.split()
		PHASES = eval(str(liste_ph))
	else:
		PHASES = ['ttall']

	if os.path.exists(model1d):
		taup.taup_create.build_taup_model(model1d,output_folder=DIR2 + '/input_files')
		modelpath = glob.glob(DIR2 + '/input_files/*.npz')
		model = TauPyModel(model=modelpath[0]) 	
	else:
		model = TauPyModel(model=model1d)
	for i in range(len(EVT)):
		for j in range(len(STA)):
			if EVT.ndim ==1:
				if STA.ndim ==1:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[1],distance_in_degree=abs(STA[0]-EVT[0]),phase_list=eval(str(PHASES)))
				else:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[1],distance_in_degree=abs(STA[j,0]-EVT[0]),phase_list=eval(str(PHASES)))
			else:
				if STA.ndim ==1:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[i,1],distance_in_degree=abs(STA[0]-EVT[i,0]),phase_list=eval(str(PHASES)))
				else:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[i,1],distance_in_degree=abs(STA[j,0]-EVT[i,0]),phase_list=eval(str(PHASES)))
			nber_phase=np.size(arrivals)
			for k in range(nber_phase):
				arrival = arrivals[k]
				path_depth = arrival.path['depth']
				if EVT.ndim ==1:
					path_deg = arrival.path['dist']*180./np.pi+EVT[0]
				else:
					path_deg = arrival.path['dist']*180./np.pi+EVT[i,0]
				auxa.plot(-(width-180)/2.+width-path_deg,6371.-path_depth,color="black",LineStyle="-",clip_on=False)

	ax2 = fig.add_axes([0.64, 0.67, .35, .4])
	map = Basemap(projection='moll', lon_0=0, resolution="c",ax=ax2)

	map.drawmapboundary(fill_color='#cccccc')
	map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)

	# fig.patch.set_alpha(0.0)

	PC=np.loadtxt(DIR + '/hotspots.xy')
	PB=np.loadtxt(DIR + '/plate_boundaries.xy')
	xPC, yPC = map(PC[:,0],PC[:,1])
	xPB, yPB = map(PB[:,0],PB[:,1])
	map.scatter(xPC, yPC, s=5, zorder=10,color="magenta", marker="o",edgecolor="magenta")
	map.scatter(xPB, yPB, s=2, zorder=10,color="0.3", marker=".")

	EQ=np.loadtxt(DIR + '/catalogue.xy')
	xEQ, yEQ = map(EQ[:,0],EQ[:,1])
	map.scatter(xEQ, yEQ, s=5, zorder=10,color="white", marker="o",edgecolor="k")

	x2, y2 = map(elon,elat)
	x3, y3 = map(slon,slat)
	map.scatter(x3, y3, s=80, zorder=10,color="red", marker="o",edgecolor="k")
	map.scatter(x2, y2, s=80, zorder=10,color="lime", marker="o",edgecolor="k")

	if (elon)>180:
		elon = elon-360
	if elon<-180:
		elon=elon+360
	if (slon)>180:
		slon = slon-360
	if slon<-180:
		slon=slon+360

	clat, clon = midpoint(elat,elon,slat,slon)
	map.drawgreatcircle(clon,clat,slon,slat,ls='None',marker='o',markersize=2,color='k')
	map.drawgreatcircle(clon,clat,elon,elat,ls='None',marker='o',markersize=2,color='k')
	x1, y1 = map(clon,clat)
	map.scatter(x1, y1, s=80, zorder=10,color="yellow", marker="o", edgecolor="k")

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)	

	plt.show()

def path_plot_fromfile(file,th,r,para,elat,elon,slat,slon,Vmax,EVT,STA,liste,model1d):

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(DIR)
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(DIR)
	# os.chdir('../../')

	shutil.copy2(os.path.abspath(file),DIR2 + '/input_files/cross_section_path.xyz')

	Z = np.loadtxt(DIR2 + '/input_files/cross_section_path.xyz')

	if para == 'VP':
		para1 = 0
	elif para =='VS':
		para1 = 1
	elif para == 'RHO':
		para1 = 2

	model = Z[:,para1+2]

	width = abs(th[-1]-th[0])
	depth = 2890

	model_new = np.reshape(model, (len(th),len(r)))
	model_new2=np.transpose(model_new)
	model_new2 = np.fliplr(model_new2)
	X, Y = np.meshgrid(th, r)
	fig = plt.figure()

	th0, th1 = (90-(width/2), 90+(width/2))
	r0, r1 = (0, 6371)
	thstep, rstep = (10, 500)

	# scale degrees to radians:
	tr_scale = Affine2D().scale(np.pi/180., 1.)
	tr = tr_scale + PolarAxes.PolarTransform()
	theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
	theta_tick_formatter = angle_helper.FormatterDMS()
	r_ticks = [(5871.,r"$500$"),(5371.,r"$1000$"),(4871.,r"$1500$"),(4371,r"$2000$"),(3891,r"$2500$"),(3480,r"$CMB$"),(1221,r"$ICB$"),(6371,r"$0$")]
	r_grid_locator  = FixedLocator([v for v, s in r_ticks])
	r_tick_formatter = DictFormatter(dict(r_ticks))
	grid_helper = GridHelperCurveLinear(tr,
	                                      extremes=(th0, th1, r0, r1),
	                                      grid_locator1=theta_grid_locator,
	                                      grid_locator2=r_grid_locator,
	                                      tick_formatter1=None,
	                                      tick_formatter2=r_tick_formatter)

	fig.subplots_adjust(bottom=0.3, left=0.05, right=0.95, top=0.84)
	fig.subplots_adjust()

	ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
	fig.add_subplot(ax)

	# adjust x axis (theta):
	ax.axis["bottom"].set_visible(True)
	ax.axis["bottom"].toggle(all=False)
	ax.axis["top"].set_visible(True)
	ax.axis["top"].toggle(all=False)
	# adjust y axis (r):
	ax.axis["left"].set_visible(True)
	ax.axis["left"].set_axis_direction("bottom") # tick direction
	ax.axis["left"].major_ticklabels.set_axis_direction("right")

	# create a parasite axes whose transData is theta, r:
	auxa = ax.get_aux_axes(tr)
	# make aux_ax to have a clip path as in a?:
	auxa.patch = ax.patch 
	# this has a side effect that the patch is drawn twice, and possibly over some other
	# artists. So, we decrease the zorder a bit to prevent this:
	ax.patch.zorder = -2
	__s = auxa.pcolormesh(X, Y, model_new2, cmap='RdYlBu', vmin=-Vmax, vmax=Vmax)


	# fig.subplots_adjust(bottom=0.2, right=0.9, left=0.1, top=0.9)
	cb = fig.add_axes([0.2, 0.1, 0.6, 0.05])
	step2 = float(Vmax)/4
	fig.colorbar(__s,orientation="horizontal", cax = cb, ticks=np.arange(-Vmax,Vmax+step2,step2))
	auxa.add_artist(plt.Circle([0, 0], radius=1221, ls='--', lw=1, color='k', fill=False,
                            transform=ax.transData._b, zorder=2))
	auxa.add_artist(plt.Circle([0, 0], radius=3480., ls='--', lw=1, color='k', fill=False,
                            transform=ax.transData._b, zorder=2))

	if STA.ndim == 1:
		auxa.scatter(-(width-180)/2.+width-(STA[0]),6371.-(STA[1]), s=80,zorder=3,\
                                           color="black", marker="^",\
                                           edgecolor="k",clip_on=False)
	else:
		auxa.scatter(-(width-180)/2.+width-(STA[:,0]),6371.-(STA[:,1]), s=80,zorder=3,\
                                           color="black", marker="^",\
                                           edgecolor="k",clip_on=False)
	if EVT.ndim ==1:
		auxa.scatter(-(width-180)/2.+width-(EVT[0]),6371.-(EVT[1]), s=80,zorder=3,\
                                           color="red", marker="*",\
                                           edgecolor="k",clip_on=False)
	else:
		auxa.scatter(-(width-180)/2.+width-(EVT[:,0]),6371.-(EVT[:,1]), s=80,zorder=3,\
                                           color="red", marker="*",\
                                           edgecolor="k",clip_on=False) 	

	if liste:
		liste_ph = liste.split()
		PHASES = eval(str(liste_ph))
	else:
		PHASES = ['ttall']

	if os.path.exists(model1d):
		taup.taup_create.build_taup_model(model1d,output_folder=DIR2 + '/input_files')
		modelpath = glob.glob(DIR2 + '/input_files/*.npz')
		model = TauPyModel(model=modelpath[0]) 	
	else:
		model = TauPyModel(model=model1d)
	for i in range(len(EVT)):
		for j in range(len(STA)):
			if EVT.ndim ==1:
				if STA.ndim ==1:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[1],distance_in_degree=abs(STA[0]-EVT[0]),phase_list=eval(str(PHASES)))
				else:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[1],distance_in_degree=abs(STA[j,0]-EVT[0]),phase_list=eval(str(PHASES)))
			else:
				if STA.ndim ==1:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[i,1],distance_in_degree=abs(STA[0]-EVT[i,0]),phase_list=eval(str(PHASES)))
				else:
					arrivals = model.get_ray_paths(source_depth_in_km=EVT[i,1],distance_in_degree=abs(STA[j,0]-EVT[i,0]),phase_list=eval(str(PHASES)))
			nber_phase=np.size(arrivals)
			for k in range(nber_phase):
				arrival = arrivals[k]
				path_depth = arrival.path['depth']
				if EVT.ndim ==1:
					path_deg = arrival.path['dist']*180./np.pi+EVT[0]
				else:
					path_deg = arrival.path['dist']*180./np.pi+EVT[i,0]
				auxa.plot(-(width-180)/2.+width-path_deg,6371.-path_depth,color="black",LineStyle="-",clip_on=False)

	ax2 = fig.add_axes([0.64, 0.67, .35, .4])
	map = Basemap(projection='moll', lon_0=0, resolution="c",ax=ax2)

	map.drawmapboundary(fill_color='#cccccc')
	map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)

	# fig.patch.set_alpha(0.0)

	PC=np.loadtxt('hotspots.xy')
	PB=np.loadtxt('plate_boundaries.xy')
	xPC, yPC = map(PC[:,0],PC[:,1])
	xPB, yPB = map(PB[:,0],PB[:,1])
	map.scatter(xPC, yPC, s=5, zorder=10,color="magenta", marker="o",edgecolor="magenta")
	map.scatter(xPB, yPB, s=2, zorder=10,color="0.3", marker=".")

	EQ=np.loadtxt('catalogue.xy')
	xEQ, yEQ = map(EQ[:,0],EQ[:,1])
	map.scatter(xEQ, yEQ, s=5, zorder=10,color="white", marker="o",edgecolor="k")

	x2, y2 = map(elon,elat)
	x3, y3 = map(slon,slat)
	map.scatter(x3, y3, s=80, zorder=10,color="red", marker="o",edgecolor="k")
	map.scatter(x2, y2, s=80, zorder=10,color="lime", marker="o",edgecolor="k")

	if (elon)>180:
		elon = elon-360
	if elon<-180:
		elon=elon+360
	if (slon)>180:
		slon = slon-360
	if slon<-180:
		slon=slon+360

	clat, clon = midpoint(elat,elon,slat,slon)
	map.drawgreatcircle(clon,clat,slon,slat,ls='None',marker='o',markersize=2,color='k')
	map.drawgreatcircle(clon,clat,elon,elat,ls='None',marker='o',markersize=2,color='k')
	x1, y1 = map(clon,clat)
	map.scatter(x1, y1, s=80, zorder=10,color="yellow", marker="o", edgecolor="k")

	path = DIR2 + "/output_files_cross"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_cross')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)
	# os.chdir('../../')
	
	plt.show()

################################################################
#                      TOOLS FOR SPECTRUM
################################################################

def spectrum(model,para,depth,NSmax):

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 1
	elif para =='VP':
		para1 = 2
	elif para == 'RHO':
		para1 = 3

	path = DIR2 + "/output_files_spectre"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_spectre')
	for i in range(len(dirs)):
		if (os.path.getmtime(dirs[i])-time.time())>10.:
			os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
 # os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth) + ' ' + str(model_choice)+ ' ' + str(NSmax)) 
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model +'_'+ para+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model +'_'+ para+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_spectre.csh ' + str(NSmax) + ' ' + str(model_choice) + ' ' + str(para1) + ' ' + str(depth))
	filename4 = DIR2 + '/output_files_spectre/spectre_' +str(model)+ '_'+str(para)+'_'+str(depth)+'km.xy'
	sp = np.loadtxt(filename4)

	path = DIR2 + "/output_files_spectre"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_spectre')
	for i in range(len(dirs)):
		os.remove(dirs[i])

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	os.chdir(cwd)

	return sp

def spectrum_GUI(model,para,depth,NSmax):

	if model == 'SEISGLOB2':
		model_choice = 1
	elif model == 'S40RTS':
		model_choice = 2
	elif model == 'SEMUCBWM1':
		model_choice = 3
	elif model == 'S362WMANIM':
		model_choice = 4
	elif model == 'SEISGLOB1':
		model_choice = 5
	elif model == 'SP12RTS':
		model_choice = 6
	elif model == 'SGLOBE':
		model_choice = 7
	elif model == '3D2016':
		model_choice = 8
	elif model == 'MYMODEL':
		model_choice = 9

	if para == 'VS':
		para1 = 1
	elif para =='VP':
		para1 = 2
	elif para == 'RHO':
		para1 = 3

	path = DIR2 + "/output_files_spectre"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_spectre')
	for i in range(len(dirs)):
		if (os.path.getmtime(dirs[i])-time.time())>10.:
			os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
 # os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth) + ' ' + str(model_choice)+ ' ' + str(NSmax)) 
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model +'_'+ para+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model +'_'+ para+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_spectre.csh ' + str(NSmax) + ' ' + str(model_choice) + ' ' + str(para1) + ' ' + str(depth))
	filename4 = DIR2 + '/output_files_spectre/spectre_' +str(model)+ '_'+str(para)+'_'+str(depth)+'km.xy'
	sp = np.loadtxt(filename4)

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		if (os.path.getmtime(dirs[i])-time.time())>30.:
			os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
 # os.chdir(DIR)

	os.chdir(DIR)

	return sp

def spectrum_fromfile_GUI(file,NSmax):

	path = DIR2 + "/output_files_spectre"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_spectre')
	for i in range(len(dirs)):
		if (os.path.getmtime(dirs[i])-time.time())>10.:
			os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_spectre_fromfile.csh ' + str(file) + ' ' + str(NSmax))
	filename4 = DIR2 + '/output_files_spectre/spectre_fromfile.xy'
	sp = np.loadtxt(filename4)
	os.chdir(DIR)

	return sp

def spectrum_fromfile(file,NSmax):

	path = DIR2 + "/output_files_spectre"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_spectre')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	shutil.copy2(os.path.abspath(file),DIR2 + '/input_files/map_file.xyz')
	file_new = DIR2 + '/input_files/map_file.xyz'

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_spectre_fromfile.csh ' + str(file_new) + ' ' + str(NSmax))
	filename4 = DIR2 + '/output_files_spectre/spectre_fromfile.xy'
	sp = np.loadtxt(filename4)

	path = DIR2 + "/output_files_spectre"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_spectre')
	for i in range(len(dirs)):
		if (os.path.getmtime(dirs[i])-time.time())>10.:
			os.remove(dirs[i])

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])

	os.chdir(cwd)

	return sp

################################################################
#                      TOOLS FOR CORRELATION
################################################################

def correlation_GUI(model1,model2,depth1,depth2,para1,para2,NSmax):

	if model1 == 'SEISGLOB2':
		model_choice1 = 1
	elif model1 == 'S40RTS':
		model_choice1 = 2
	elif model1 == 'SEMUCBWM1':
		model_choice1 = 3
	elif model1 == 'S362WMANIM':
		model_choice1 = 4
	elif model1 == 'SEISGLOB1':
		model_choice1 = 5
	elif model1 == 'SP12RTS':
		model_choice1 = 6
	elif model1 == 'SGLOBE':
		model_choice1 = 7
	elif model1 == '3D2016':
		model_choice1 = 8
	elif model1 == 'MYMODEL':
		model_choice1 = 9

	if model2 == 'SEISGLOB2':
		model_choice2 = 1
	elif model2 == 'S40RTS':
		model_choice2 = 2
	elif model2 == 'SEMUCBWM1':
		model_choice2 = 3
	elif model2 == 'S362WMANIM':
		model_choice2 = 4
	elif model2 == 'SEISGLOB1':
		model_choice2 = 5
	elif model2 == 'SP12RTS':
		model_choice2 = 6
	elif model2 == 'SGLOBE':
		model_choice2 = 7
	elif model2 == '3D2016':
		model_choice2 = 8
	elif model2 == 'MYMODEL':
		model_choice2 = 9

	if para1 == 'VS':
		para11 = 1
	elif para1 =='VP':
		para11 = 2
	elif para1 == 'RHO':
		para11 = 3

	if para2 == 'VS':
		para22 = 1
	elif para2 =='VP':
		para22 = 2
	elif para2 == 'RHO':
		para22 = 3

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth1) + ' ' + str(model_choice1)+ ' ' + str(NSmax) )
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model1 +'_'+ para1+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model1 +'_'+ para1+'_'+ str(depth1)+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_map.csh ' +  str(depth2) + ' ' + str(model_choice2)+ ' ' + str(NSmax) )
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model2 +'_'+ para2+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model2 +'_'+ para2+'_'+ str(depth2)+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_corr.csh ' + str(NSmax) + ' ' + str(depth1) + ' ' + str(model_choice1) + ' ' + str(para11) + ' ' + str(depth2)+ ' ' + str(model_choice2) + ' ' + str(para22))
	filename4 = DIR2 + '/output_files_corr/corr_'+str(model1)+'_'+str(model2)+'.xy'
	corr = np.loadtxt(filename4)

	os.chdir(cwd)

	return corr

def correlation(model1,model2,depth1,depth2,para1,para2,NSmax):

	if model1 == 'SEISGLOB2':
		model_choice1 = 1
	elif model1 == 'S40RTS':
		model_choice1 = 2
	elif model1 == 'SEMUCBWM1':
		model_choice1 = 3
	elif model1 == 'S362WMANIM':
		model_choice1 = 4
	elif model1 == 'SEISGLOB1':
		model_choice1 = 5
	elif model1 == 'SP12RTS':
		model_choice1 = 6
	elif model1 == 'SGLOBE':
		model_choice1 = 7
	elif model1 == '3D2016':
		model_choice1 = 8
	elif model1 == 'MYMODEL':
		model_choice1 = 9

	if model2 == 'SEISGLOB2':
		model_choice2 = 1
	elif model2 == 'S40RTS':
		model_choice2 = 2
	elif model2 == 'SEMUCBWM1':
		model_choice2 = 3
	elif model2 == 'S362WMANIM':
		model_choice2 = 4
	elif model2 == 'SEISGLOB1':
		model_choice2 = 5
	elif model2 == 'SP12RTS':
		model_choice2 = 6
	elif model2 == 'SGLOBE':
		model_choice2 = 7
	elif model2 == '3D2016':
		model_choice2 = 8
	elif model2 == 'MYMODEL':
		model_choice2 = 9

	if para1 == 'VS':
		para11 = 1
	elif para1 =='VP':
		para11 = 2
	elif para1 == 'RHO':
		para11 = 3

	if para2 == 'VS':
		para22 = 1
	elif para2 =='VP':
		para22 = 2
	elif para2 == 'RHO':
		para22 = 3

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth1) + ' ' + str(model_choice1)+ ' ' + str(NSmax) )
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model1 +'_'+ para1+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model1 +'_'+ para1+'_'+ str(depth1)+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_map.csh ' +  str(depth2) + ' ' + str(model_choice2)+ ' ' + str(NSmax) )
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model2 +'_'+ para2+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model2 +'_'+ para2+'_'+ str(depth2)+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_corr.csh ' + str(NSmax) + ' ' + str(depth1) + ' ' + str(model_choice1) + ' ' + str(para11) + ' ' + str(depth2)+ ' ' + str(model_choice2) + ' ' + str(para22))
	filename4 = DIR2 + '/output_files_corr/corr_'+str(model1)+'_'+str(model2)+'.xy'
	corr = np.loadtxt(filename4)

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)
	os.chdir(cwd)

	return corr

def correlation_fromfile_GUI(model1,depth1,para1,file,NSmax):

	if model1 == 'SEISGLOB2':
		model_choice1 = 1
	elif model1 == 'S40RTS':
		model_choice1 = 2
	elif model1 == 'SEMUCBWM1':
		model_choice1 = 3
	elif model1 == 'S362WMANIM':
		model_choice1 = 4
	elif model1 == 'SEISGLOB1':
		model_choice1 = 5
	elif model1 == 'SP12RTS':
		model_choice1 = 6
	elif model1 == 'SGLOBE':
		model_choice1 = 7
	elif model1 == '3D2016':
		model_choice1 = 8
	elif model1 == 'MYMODEL':
		model_choice1 = 9

	if para1 == 'VS':
		para11 = 1
	elif para1 =='VP':
		para11 = 2
	elif para1 == 'RHO':
		para11 = 3

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth1) + ' ' + str(model_choice1)+ ' ' + str(NSmax) )
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model1 +'_'+ para1+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model1 +'_'+ para1+'_'+ str(depth1)+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_corr_fromfile_bis.csh ' + str(file) +' ' + str(NSmax) + ' ' +str(depth1) + ' ' + str(model_choice1) + ' ' + str(para11))
	filename4 = DIR2 + '/output_files_corr/corr_fromfile_'+ str(model1) +'.xy'
	corr = np.loadtxt(filename4)
	os.chdir(DIR)

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	return corr

def correlation_fromfile2_GUI(file1,file2,NSmax):

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_corr_fromfile.csh ' + str(file1) +' ' + str(file2) + ' ' + str(NSmax))
	filename4 = DIR2 + '/output_files_corr/corr_fromfile.xy'
	corr = np.loadtxt(filename4)

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	return corr

def correlation_fromfile(model1,depth1,para1,file,NSmax):

	if model1 == 'SEISGLOB2':
		model_choice1 = 1
	elif model1 == 'S40RTS':
		model_choice1 = 2
	elif model1 == 'SEMUCBWM1':
		model_choice1 = 3
	elif model1 == 'S362WMANIM':
		model_choice1 = 4
	elif model1 == 'SEISGLOB1':
		model_choice1 = 5
	elif model1 == 'SP12RTS':
		model_choice1 = 6
	elif model1 == 'SGLOBE':
		model_choice1 = 7
	elif model1 == '3D2016':
		model_choice1 = 8
	elif model1 == 'MYMODEL':
		model_choice1 = 9

	if para1 == 'VS':
		para11 = 1
	elif para1 =='VP':
		para11 = 2
	elif para1 == 'RHO':
		para11 = 3

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	shutil.copy2(os.path.abspath(file),DIR2 + '/input_files/map_file1.xyz')
	file_new = DIR2 + '/input_files/map_file1.xyz'

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_map.csh ' +  str(depth1) + ' ' + str(model_choice1)+ ' ' + str(NSmax) )
	filename5 = DIR2 + '/output_files_map/map_NEW_'+ model1 +'_'+ para1+'.xyz'
	filename6 = DIR2 + '/input_files/map_NEW_'+ model1 +'_'+ para1+'_'+ str(depth1)+'.xyz'
	shutil.copy2(filename5,filename6)
	os.system('csh make_corr_fromfile_bis.csh ' + str(file_new) +' ' + str(NSmax) + ' ' +str(depth1) + ' ' + str(model_choice1) + ' ' + str(para11))
	filename4 = DIR2 + '/output_files_corr/corr_fromfile_'+ str(model1) +'.xy'
	corr = np.loadtxt(filename4)

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)
	os.chdir(cwd)

	return corr

def correlation_fromfile2(file1,file2,NSmax):

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	shutil.copy2(os.path.abspath(file1),DIR2 + '/input_files/map_file1.xyz')
	file1_new = DIR2 + '/input_files/map_file1.xyz'
	shutil.copy2(os.path.abspath(file2),DIR2 + '/input_files/map_file2.xyz')
	file2_new = DIR2 + '/input_files/map_file2.xyz'

	os.chdir(DIR + '/fortran_files/')
	os.system('csh make_corr_fromfile.csh ' + str(file1_new) +' ' + str(file2_new) + ' ' + str(NSmax))
	filename4 = DIR2 + '/output_files_corr/corr_fromfile.xy'
	corr = np.loadtxt(filename4)

	path = DIR2 + "/output_files_corr"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_corr')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir('../../')

	path = DIR2 + "/input_files"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/input_files')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)

	path = DIR2 + "/output_files_map"
	dirs = os.listdir(path)
	os.chdir(DIR2 + '/output_files_map')
	for i in range(len(dirs)):
		os.remove(dirs[i])
	# os.chdir(DIR)
	os.chdir(cwd)

	return corr

################################################################
#                      TOOLS FOR TIME
################################################################

def get_travel_time_list(model1,elat,elon,edepth,slat,slon,PHASES,PATH_INTERPOLATION = False):

	if (elon)>180:
		elon = elon-360
	if elon<-180:
		elon=elon+360
	if (slon)>180:
		slon = slon-360
	if slon<-180:
		slon=slon+360

	center_lat2,center_lon2 = midpoint(elat,elon,slat,slon)
	great_circle_dist, azi2, baz2 = gps2dist_azimuth(int(round(elat)),int(round(elon)),int(round(slat)),int(round(slon)))
	great_circle_dist /=  111194.9
	KM_PER_DEG = 111.1949
	width2, az, baz = gps2dist_azimuth(int(round(center_lat2)),int(round(center_lon2)),int(round(slat)),int(round(slon)))

	cross_section_midpoint(model1,'VS',int(center_lat2),int(center_lon2),int(az),2890,int(great_circle_dist),40)

	model_name = 'prem' # select the 1D model
	model = TauPyModel(model='Taup_models/prem/prem_Obspy_all_interpolated.npz')

	# if os.path.exists(modelpath):
	# 	taup.taup_create.build_taup_model(modelpath,output_folder=DIR2 + '/input_files')
	# 	modelpath2 = glob.glob(DIR2 + '/input_files/*.npz')
	# 	model2 = TauPyModel(model=modelpath2[0]) 	
	# else:
		# model2 = TauPyModel(model=modelpath)	

	if model1 == 'SEISGLOB2':
		model_choice1 = 'prem'
	elif model1 == 'S40RTS':
		model_choice1 = 'prem'
	elif model1 == 'SEMUCBWM1':
		model_choice1 = 'SEMUCBWM1'
	elif model1 == 'S362WMANIM':
		model_choice1 = 'S362WMANIM'
	elif model1 == 'SEISGLOB1':
		model_choice1 = 'prem'
	elif model1 == 'SP12RTS':
		model_choice1 = 'prem'
	elif model1 == 'SGLOBE':
		model_choice1 = 'prem'

	phase_name = []

	model2D = np.loadtxt(DIR2 + '/output_files_cross/' + str(model1) + '_PathPy.sph')

	p_phases = ['p','P','Pn','Pdiff','PKP','PKiKP','PKIKP','PcP','pP','pPdiff','pPKP','pPKiKP','PKKP','PKIKKIKP','PP','PKPPKP','PKIKPPKIKP','pPKIKP']
	s_phases = ['s','S','Sn','Sdiff','sS','sSdiff','ScS','SS']
	not_added_converted_phases = ['SKIKS','SKKS','SKIKKIKS','SKSSKS','SKIKSSKIKS','SKiKP','SKP','SKIKP','SKKP','SKIKKIKP','PKS','PKIKS','PKKS','PKIKKIKS','PcS','ScP','SP','PS','sSKS','sSKIKS','pS','pSdiff','pSKS','pSKIKS','sP','sPdiff','sPKP','sPKIKP','sPKiKP']
	added_converted_phases = ['SKS','SKKS']

	# load the interpolated PREM model
	PREM = np.loadtxt('Taup_models/prem/prem_interpolated.txt')
	PREM_depth = PREM[:,0] # depth in km
	PREM_VP = PREM[:,1] # km/s
	PREM_VS = PREM[:,2] # km/s
	PREM_RHO = PREM[:,3]
	PREM_r = 6371.0 - PREM_depth # radius in km
	# load the interpolated PREM mantle  model
	PREM_mantle = np.loadtxt('Taup_models/prem/prem_mantle_interpolated.txt')
	PREM_mantle_depth = PREM_mantle[:,0] # depth in km
	PREM_mantle_VP = PREM_mantle[:,1] # km/s
	PREM_mantle_VS = PREM_mantle[:,2] # km/s
	PREM_mantle_RHO = PREM_mantle[:,3]
	PREM_mantle_r = 6371.0 - PREM_mantle_depth # radius in km
	# load the interpolated PREM outercore  model
	PREM_outercore = np.loadtxt('Taup_models/prem/prem_outercore_interpolated.txt')
	PREM_outercore_depth = PREM_outercore[:,0] # depth in km
	PREM_outercore_VP = PREM_outercore[:,1] # km/s
	PREM_outercore_VS = PREM_outercore[:,2] # km/s
	PREM_outercore_RHO = PREM_outercore[:,3]
	PREM_outercore_r = 6371.0 - PREM_outercore_depth # radius in km
	# load the interpolated PREM innercore model
	PREM_innercore = np.loadtxt('Taup_models/prem/prem_innercore_interpolated.txt')
	PREM_innercore_depth = PREM_innercore[:,0] # depth in km
	PREM_innercore_VP = PREM_innercore[:,1] # km/s
	PREM_innercore_VS = PREM_innercore[:,2] # km/s
	PREM_innercore_RHO = PREM_innercore[:,3]
	PREM_innercore_r = 6371.0 - PREM_innercore_depth # radius in km

	# load the interpolated PREM model
	REF = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_interpolated.txt')
	REF_depth = REF[:,0] # depth in km
	REF_VP = REF[:,1] # km/s
	REF_VS = REF[:,2] # km/s
	REF_RHO = REF[:,3]
	REF_r = 6371.0 - REF_depth # radius in km
	# load the interpolated REF mantle  model
	REF_mantle = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_mantle_interpolated.txt')
	REF_mantle_depth = REF_mantle[:,0] # depth in km
	REF_mantle_VP = REF_mantle[:,1] # km/s
	REF_mantle_VS = REF_mantle[:,2] # km/s
	REF_mantle_RHO = REF_mantle[:,3]
	REF_mantle_r = 6371.0 - REF_mantle_depth # radius in km
	# load the interpolated REF outercore  model
	REF_outercore = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_outercore_interpolated.txt')
	REF_outercore_depth = REF_outercore[:,0] # depth in km
	REF_outercore_VP = REF_outercore[:,1] # km/s
	REF_outercore_VS = REF_outercore[:,2] # km/s
	REF_outercore_RHO = REF_outercore[:,3]
	REF_outercore_r = 6371.0 - REF_outercore_depth # radius in km
	# load the interpolated REF innercore model
	REF_innercore = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_innercore_interpolated.txt')
	REF_innercore_depth = REF_innercore[:,0] # depth in km
	REF_innercore_VP = REF_innercore[:,1] # km/s
	REF_innercore_VS = REF_innercore[:,2] # km/s
	REF_innercore_RHO = REF_innercore[:,3]
	REF_innercore_r = 6371.0 - REF_innercore_depth # radius in km

	MODEL_2D = True

	# compute the distance from the event to the station
	ev_distance = locations2degrees(elat,elon,slat,slon)
	arrivals = model.get_ray_paths(source_depth_in_km=edepth,distance_in_degree=ev_distance,phase_list=PHASES)
	tt2D = np.zeros(len(arrivals))
	dtt2D = np.zeros(len(arrivals))
	ttREF = np.zeros(len(arrivals))

	for i_phase in range(len(arrivals)):

		arrival = arrivals[i_phase]

		if arrival.name not in p_phases and arrival.name not in s_phases and arrival.name not in added_converted_phases:
			phase_name.append('NO PHASE')
			scr_msg = 'The phase ' + arrival.name + ' can not be computed. Sorry!'
			print(scr_msg)
			print('---------------------------------------------------------------')
			tt2D[i_phase] = 0.
			dtt2D[i_phase] = 0.
		else:
			phase_name.append(arrival.name)
			theor_tt = arrival.time # theoretical travel time in seconds
			depth = arrival.path['depth'] # km
			r = 6371.0 - depth # radius - km
			r_deg = r * KM_PER_DEG # radius in deg
			theta_rad = arrival.path['dist'] # distance in radians
			theta_deg = theta_rad*180./np.pi # distance in degrees

			# interpolate r and theta given by taup to refine depending the given 2D model
			if PATH_INTERPOLATION: 
				N = 10000
				theta_coords = np.linspace(theta_rad.min(), theta_rad.max(), N)
				r_inter = np.interp(theta_coords, theta_rad, r)
				check_inter_r = False
				if check_inter_r:
					plt.plot(theta_coords,r_inter)
					plt.show()
				r1 = r_inter[:-1]
				r2 = r_inter[1:]
				t1 = theta_coords[:-1]
				t2 = theta_coords[1:]
				dS = np.zeros(N-1)
				dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
				model_vel = np.zeros(N-1) # velocity model for the ray path
				r_coords = 6370-depth
				for i in range(N-1): # fill the velocity model
				 r_mean = (r_inter[i] + r_inter[i+1])/2 # take the middle of the cell
				 if r_mean < 6371.0 and r_mean > 3480: #  mantle
				  index = np.argmin(np.abs(PREM_mantle_r-r_mean))
				  if arrival.name in p_phases:
				  	model_vel[i] = PREM_mantle_VP[index] 
				  elif arrival.name in s_phases:
				  	model_vel[i] = PREM_mantle_VS[index]
				  elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				  	model_vel[i] = PREM_mantle_VS[index]
				  elif r_mean <= 3480.0 and r_mean > 1222:
				  	index = np.argmin(np.abs(PREM_outercore_r-r_mean))
				  if arrival.name in p_phases: # outercore
				   model_vel[i] = PREM_outercore_VP[index]
				  elif arrival.name in s_phases:
				  	model_vel[i] = PREM_outercore_VS[index]
				  elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				  	model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values 
				  elif r_mean <= 1222.0 and r_mean >= 0: # innercore
				   index = np.argmin(np.abs(PREM_innercore_r-r_mean))
				  if arrival.name in p_phases:
				  	model_vel[i] = PREM_innercore_VP[index] 
				  elif arrival.name in s_phases:
				  	model_vel[i] = PREM_innercore_VS[index] 
			else: # no path interpolation
			 r1 = r[:-1]
			 r2 = r[1:]
			 rtot = r1+r2
			 t1 = theta_rad[:-1]
			 t2 = theta_rad[1:]
			 dS = np.zeros(len(r)-1)
			 dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
			 model_vel = np.zeros(len(r)-1) # velocity model for the ray path
			 for i in range(len(r)-1): # fill the velocity model
			  r_mean = (r[i] + r[i+1])/2 # take the middle of the cell
			  if r_mean < 6371.0 and r_mean >= 3480: #  mantle
			   index = np.argmin(np.abs(PREM_mantle_r-r_mean))
			   if arrival.name in p_phases:
			   	model_vel[i] = PREM_mantle_VP[index] # no interpolation of PREM values
			   elif arrival.name in s_phases:
			   	model_vel[i] = PREM_mantle_VS[index] # no interpolation of PREM values
			   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
			   	model_vel[i] = PREM_mantle_VS[index] # no interpolation of PREM values
			   elif r_mean < 3480.0 and r_mean > 1222:
			   	index = np.argmin(np.abs(PREM_outercore_r-r_mean))
			   if arrival.name in p_phases: # outercore
			    model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values
			   elif arrival.name in s_phases:
			   	model_vel[i] = PREM_outercore_VS[index] # no interpolation of PREM values
			   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
			   	model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values
			   elif r_mean <= 1222.0 and r_mean >= 0: # innercore
			    index = np.argmin(np.abs(PREM_innercore_r-r_mean))
			   if arrival.name in p_phases:
			   	model_vel[i] = PREM_innercore_VP[index] # no interpolation of PREM values
			   elif arrival.name in s_phases:
			   	model_vel[i] = PREM_innercore_VS[index] # no interpolation of PREM values

			PREM_tt = 0.    # store the PREM numerical computed travel time
			for i in range(len(dS)):
				if model_vel[i] !=0:
					PREM_tt = PREM_tt + dS[i]/model_vel[i]
					# PREM_tt = np.sum(dS/model_vel)

			if MODEL_2D:
				# initialize variables
				model2D_r = np.zeros([len(model2D),1])
				model2D_theta = np.zeros([len(model2D),1])
				model2D_VP = np.zeros([len(model2D),1])
				model2D_VS = np.zeros([len(model2D),1])
				model2D_RHO = np.zeros([len(model2D),1])
				model2D_coords = np.zeros([len(model2D),2])
				# fill the matrices
				model2D_r = model2D[:,0] # km
				model2D_theta = model2D[:,1] # deg
				model2D_VP = model2D[:,2]
				model2D_VS = model2D[:,3]
				model2D_RHO = model2D[:,4]
				model2D_coords[:,0] = model2D[:,0] # rad 
				model2D_coords[:,1] = model2D[:,1] # theta
				# interpolate r and theta given by taupy to refine depending the given 2D model
				N = 1500
				theta_coords = np.linspace(theta_rad.min(), theta_rad.max(), N)
				r_inter = np.zeros([N,1])
				r_inter = np.interp(theta_coords, theta_rad, r)
				path_coords = np.zeros([N,2])
				path_coords[:,0] = r_inter
				path_coords[:,1] = theta_coords
				r1 = r_inter[:-1]
				r2 = r_inter[1:]
				t1 = theta_coords[:-1]
				t2 = theta_coords[1:]
				dS = np.zeros(N-1)
				dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
				path_coords[:,1] *= 180./np.pi # convert to deg
				r_coords = 6370-depth
				model_vel = np.zeros(N-1) # velocity model for the ray path
				model_vel_REF = np.zeros(N-1) # velocity model for the ray path
				for i in range(N-1): # fill the velocity model
				 r_mean = (path_coords[i,0] + path_coords[i+1,0])/2 # take the middle of the cell
				 if r_mean < 6371.0 and r_mean >= 3480.: #  mantle
				  index = np.argmin(np.abs(REF_mantle_r-r_mean))
				  if arrival.name in p_phases:
				  	REF_vel = REF_mantle_VP[index]
				  	sample_pt = path_coords[i,:]
				  	min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
				  	index = distance.cdist([min_loc],model2D_coords).argmin()
				  	# model_vel[i] = REF_vel*(1.0+model2D_VP[index]/100.) 
				  	# NEW VERSION JUNE 2017
				  	# S. Durand
				  	model_vel[i] = (model2D_VP[index]/100.)
				  	model_vel_REF[i] = REF_vel
				  elif arrival.name in s_phases:
				  	REF_vel = REF_mantle_VS[index]
				  	sample_pt = path_coords[i,:]
				  	min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
				  	index = distance.cdist([min_loc],model2D_coords).argmin()
				  	# model_vel[i] = REF_vel*(1.0+model2D_VS[index]/100.)
				  	# NEW VERSION JUNE 2017
				  	# S. Durand
				  	model_vel[i] = (model2D_VS[index]/100.)
				  	model_vel_REF[i] = REF_vel
				  elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				  	REF_vel = REF_mantle_VS[index]
				  	sample_pt = path_coords[i,:]
				  	min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
				  	index = distance.cdist([min_loc],model2D_coords).argmin()
				  	# model_vel[i] = REF_vel*(1.0+model2D_VS[index]/100.)
				  	# NEW VERSION JUNE 2017
				  	# S. Durand
				  	model_vel[i] = (model2D_VS[index]/100.)
				  	model_vel_REF[i] = REF_vel
				 elif r_mean < 3480.0 and r_mean > 1221.5:
				 	index = np.argmin(np.abs(REF_outercore_r-r_mean))
				 	# NEW VERSION JUNE 2017
				 	# S. Durand
				 	model_vel[i] = 0.
				 	if arrival.name in p_phases: # outercore
				 	 model_vel_REF[i] = REF_outercore_VP[index]
				 	elif arrival.name in s_phases:
				 		model_vel_REF[i] = REF_outercore_VS[index]
				 	elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				 		model_vel_REF[i] = REF_outercore_VP[index]
				 elif r_mean <= 1221.5 and r_mean >= 0: # innercore
				  # NEW VERSION JUNE 2017
				  # S. Durand
				  model_vel[i] = 0.
				  index = np.argmin(np.abs(REF_innercore_r-r_mean))
				  if arrival.name in p_phases:
				  	model_vel_REF[i] = REF_innercore_VP[index]
				  elif arrival.name in s_phases:
				  	model_vel_REF[i] = REF_innercore_VS[index]

			for i in range(len(dS)):
				if model_vel_REF[i] !=0:
					tt2D[i_phase] = tt2D[i_phase] - dS[i]*model_vel[i]/((1+model_vel[i])*model_vel_REF[i])
					ttREF[i_phase] = ttREF[i_phase] + dS[i]/model_vel_REF[i]
			dtt2D[i_phase] = abs(theor_tt-PREM_tt)

			path = DIR2 + "/input_files"
			dirs = os.listdir(path)
			os.chdir(DIR2 + '/input_files')
			for i in range(len(dirs)):
				os.remove(dirs[i])
			os.chdir(cwd)

	return tt2D, dtt2D, ttREF, phase_name


            # travel time for the 2D model
            # for i in range(len(dS)):
            # 	if model_vel[i] !=0:
            # 		tt2D[i_phase] = tt2D[i_phase] + dS[i]/model_vel[i]
            # 		dtt2D[i_phase] = abs(theor_tt-PREM_tt)
            # tt2D[i_phase] = np.sum(dS/model_vel)
            # dtt2D[i_phase] = abs(theor_tt-PREM_tt)

            # travel time for the 2D model

    # return tt2D, dtt2D, phase_name

def get_travel_time(model1,elat,elon,edepth,slat,slon,PHASES,PATH_INTERPOLATION = False):

				if (elon)>180:
					elon = elon-360
				if elon<-180:
					elon=elon+360
				if (slon)>180:
					slon = slon-360
				if slon<-180:
					slon=slon+360 

				# if os.path.exists(modelpath):
				# 	taup.taup_create.build_taup_model(modelpath,output_folder=DIR2 + '/input_files')
				# 	modelpath2 = glob.glob(DIR2 + '/input_files/*.npz')
				# 	model2 = TauPyModel(model=modelpath2[0]) 	
				# else:
				# 	model2 = TauPyModel(model=modelpath)	

				center_lat2,center_lon2 = midpoint(elat,elon,slat,slon)
				great_circle_dist, azi2, baz2 = gps2dist_azimuth(int(round(elat)),int(round(elon)),int(round(slat)),int(round(slon)))
				great_circle_dist /=  111194.9
				KM_PER_DEG = 111.1949
				width2, az, baz = gps2dist_azimuth(int(round(center_lat2)),int(round(center_lon2)),int(round(slat)),int(round(slon)))

				cross_section_midpoint(model1,'VS',int(center_lat2),int(center_lon2),int(az),2890,int(great_circle_dist),40)

				model_name = 'prem' # select the 1D model
				model = TauPyModel(model='Taup_models/prem/prem_Obspy_all_interpolated.npz')

				phase_name = []

				model2D = np.loadtxt(DIR2 + '/output_files_cross/' + str(model1) + '_PathPy.sph')

				p_phases = ['p','P','Pn','Pdiff','PKP','PKiKP','PKIKP','PcP','pP','pPdiff','pPKP','pPKiKP','PKKP','PKIKKIKP','PP','PKPPKP','PKIKPPKIKP','pPKIKP']
				s_phases = ['s','S','Sn','Sdiff','sS','sSdiff','ScS','SS']
				not_added_converted_phases = ['SKIKS','SKKS','SKIKKIKS','SKSSKS','SKIKSSKIKS','SKiKP','SKP','SKIKP','SKKP','SKIKKIKP','PKS','PKIKS','PKKS','PKIKKIKS','PcS','ScP','SP','PS','sSKS','sSKIKS','pS','pSdiff','pSKS','pSKIKS','sP','sPdiff','sPKP','sPKIKP','sPKiKP']
				added_converted_phases = ['SKS','SKKS']

				if model1 == 'SEISGLOB2':
					model_choice1 = 'prem'
				elif model1 == 'S40RTS':
					model_choice1 = 'prem'
				elif model1 == 'SEMUCBWM1':
					model_choice1 = 'SEMUCBWM1'
				elif model1 == 'S362WMANIM':
					model_choice1 = 'S362WMANIM'
				elif model1 == 'SEISGLOB1':
					model_choice1 = 'prem'
				elif model1 == 'SP12RTS':
					model_choice1 = 'prem'
				elif model1 == 'SGLOBE':
					model_choice1 = 'prem'

				# load the interpolated PREM model
				PREM = np.loadtxt('Taup_models/prem/prem_interpolated.txt')
				PREM_depth = PREM[:,0] # depth in km
				PREM_VP = PREM[:,1] # km/s
				PREM_VS = PREM[:,2] # km/s
				PREM_RHO = PREM[:,3]
				PREM_r = 6371.0 - PREM_depth # radius in km
				# load the interpolated PREM mantle  model
				PREM_mantle = np.loadtxt('Taup_models/prem/prem_mantle_interpolated.txt')
				PREM_mantle_depth = PREM_mantle[:,0] # depth in km
				PREM_mantle_VP = PREM_mantle[:,1] # km/s
				PREM_mantle_VS = PREM_mantle[:,2] # km/s
				PREM_mantle_RHO = PREM_mantle[:,3]
				PREM_mantle_r = 6371.0 - PREM_mantle_depth # radius in km
				# load the interpolated PREM outercore  model
				PREM_outercore = np.loadtxt('Taup_models/prem/prem_outercore_interpolated.txt')
				PREM_outercore_depth = PREM_outercore[:,0] # depth in km
				PREM_outercore_VP = PREM_outercore[:,1] # km/s
				PREM_outercore_VS = PREM_outercore[:,2] # km/s
				PREM_outercore_RHO = PREM_outercore[:,3]
				PREM_outercore_r = 6371.0 - PREM_outercore_depth # radius in km
				# load the interpolated PREM innercore model
				PREM_innercore = np.loadtxt('Taup_models/prem/prem_innercore_interpolated.txt')
				PREM_innercore_depth = PREM_innercore[:,0] # depth in km
				PREM_innercore_VP = PREM_innercore[:,1] # km/s
				PREM_innercore_VS = PREM_innercore[:,2] # km/s
				PREM_innercore_RHO = PREM_innercore[:,3]
				PREM_innercore_r = 6371.0 - PREM_innercore_depth # radius in km

				# load the interpolated PREM model
				REF = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_interpolated.txt')
				REF_depth = REF[:,0] # depth in km
				REF_VP = REF[:,1] # km/s
				REF_VS = REF[:,2] # km/s
				REF_RHO = REF[:,3]
				REF_r = 6371.0 - REF_depth # radius in km
				# load the interpolated REF mantle  model
				REF_mantle = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_mantle_interpolated.txt')
				REF_mantle_depth = REF_mantle[:,0] # depth in km
				REF_mantle_VP = REF_mantle[:,1] # km/s
				REF_mantle_VS = REF_mantle[:,2] # km/s
				REF_mantle_RHO = REF_mantle[:,3]
				REF_mantle_r = 6371.0 - REF_mantle_depth # radius in km
				# load the interpolated REF outercore  model
				REF_outercore = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_outercore_interpolated.txt')
				REF_outercore_depth = REF_outercore[:,0] # depth in km
				REF_outercore_VP = REF_outercore[:,1] # km/s
				REF_outercore_VS = REF_outercore[:,2] # km/s
				REF_outercore_RHO = REF_outercore[:,3]
				REF_outercore_r = 6371.0 - REF_outercore_depth # radius in km
				# load the interpolated REF innercore model
				REF_innercore = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_innercore_interpolated.txt')
				REF_innercore_depth = REF_innercore[:,0] # depth in km
				REF_innercore_VP = REF_innercore[:,1] # km/s
				REF_innercore_VS = REF_innercore[:,2] # km/s
				REF_innercore_RHO = REF_innercore[:,3]
				REF_innercore_r = 6371.0 - REF_innercore_depth # radius in km

				MODEL_2D = True

				# compute the distance from the event to the station
				ev_distance = locations2degrees(elat,elon,slat,slon)
				arrivals = model.get_ray_paths(source_depth_in_km=edepth,distance_in_degree=ev_distance,phase_list=PHASES)
				tt2D = 0.
				dtt2D = 0.
				ttREF = 0.

				arrival = arrivals[0]

				phase_name.append(arrival.name)
				theor_tt = arrival.time # theoretical travel time in seconds
				depth = arrival.path['depth'] # km
				r = 6371.0 - depth # radius - km
				r_deg = r * KM_PER_DEG # radius in deg
				theta_rad = arrival.path['dist'] # distance in radians
				theta_deg = theta_rad*180./np.pi # distance in degrees

				if PATH_INTERPOLATION: # interpolate r and theta given by taup to refine depending the given 2D model
				    N = 10000
				    theta_coords = np.linspace(theta_rad.min(), theta_rad.max(), N)
				    r_inter = np.interp(theta_coords, theta_rad, r)
				    r1 = r_inter[:-1]
				    r2 = r_inter[1:]
				    t1 = theta_coords[:-1]
				    t2 = theta_coords[1:]
				    dS = np.zeros(N-1)
				    dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
				    model_vel = np.zeros(N-1) # velocity model for the ray path
				    r_coords = 6370-depth
				    for i in range(N-1): # fill the velocity model
				        r_mean = (r_inter[i] + r_inter[i+1])/2 # take the middle of the cell
				        if r_mean < 6371.0 and r_mean > 3480.0: #  mantle
				            index = np.argmin(np.abs(PREM_mantle_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_mantle_VP[index] 
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_mantle_VS[index]
				            elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				            	   model_vel[i] = PREM_mantle_VS[index]
				        elif r_mean <= 3480.0 and r_mean > 1221.5:
				        	   index = np.argmin(np.abs(PREM_outercore_r-r_mean))
				        	   if arrival.name in p_phases: # outercore
				        	       model_vel[i] = PREM_outercore_VP[index] 
				        	   elif arrival.name in s_phases:
				        	   	   model_vel[i] = PREM_outercore_VS[index]
				        	   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				        	   	   model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values 
				        elif r_mean <= 1221.5 and r_mean >= 0: # innercore
				            index = np.argmin(np.abs(PREM_innercore_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_innercore_VP[index] 
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_innercore_VS[index] 

				else: # no path interpolation
				    r1 = r[:-1]
				    r2 = r[1:]
				    rtot = r1+r2
				    t1 = theta_rad[:-1]
				    t2 = theta_rad[1:]
				    dS = np.zeros(len(r)-1)
				    dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
				    model_vel = np.zeros(len(r)-1) # velocity model for the ray path
				    for i in range(len(r)-1): # fill the velocity model
				        r_mean = (r[i] + r[i+1])/2 # take the middle of the cell
				        if r_mean < 6371.0 and r_mean >= 3480.: #  mantle
				            index = np.argmin(np.abs(PREM_mantle_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_mantle_VP[index] # no interpolation of PREM values
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_mantle_VS[index] # no interpolation of PREM values
				            elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				            	   model_vel[i] = PREM_mantle_VS[index] # no interpolation of PREM values
				        elif r_mean < 3480.0 and r_mean > 1221.5:
				            index = np.argmin(np.abs(PREM_outercore_r-r_mean))
				            if arrival.name in p_phases: # outercore
				                model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_outercore_VS[index] # no interpolation of PREM values
				            elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				            	   model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values
				        elif r_mean <= 1221.5 and r_mean >= 0: # innercore
				            index = np.argmin(np.abs(PREM_innercore_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_innercore_VP[index] # no interpolation of PREM values
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_innercore_VS[index] # no interpolation of PREM values

				# store the PREM numerical computed travel time
				PREM_tt = 0.
				for i in range(len(dS)):
					if model_vel[i] !=0:
						PREM_tt = PREM_tt + dS[i]/model_vel[i]

				if MODEL_2D:
					   # initialize variables
					   model2D_r = np.zeros([len(model2D),1])
					   model2D_theta = np.zeros([len(model2D),1])
					   model2D_VP = np.zeros([len(model2D),1])
					   model2D_VS = np.zeros([len(model2D),1])
					   model2D_RHO = np.zeros([len(model2D),1])
					   model2D_coords = np.zeros([len(model2D),2])
					   # fill the matrices
					   model2D_r = model2D[:,0] # km
					   model2D_theta = model2D[:,1] # deg
					   model2D_VP = model2D[:,2]
					   model2D_VS = model2D[:,3]
					   model2D_RHO = model2D[:,4]
					   model2D_coords[:,0] = model2D[:,0] # rad 
					   model2D_coords[:,1] = model2D[:,1] # theta
					   # interpolate r and theta given by taupy to refine depending the given 2D model
					   N = 1500
					   theta_coords = np.linspace(theta_rad.min(), theta_rad.max(), N)
					   r_inter = np.zeros([N,1])
					   r_inter = np.interp(theta_coords, theta_rad, r)
					   path_coords = np.zeros([N,2])
					   path_coords[:,0] = r_inter
					   path_coords[:,1] = theta_coords
					   r1 = r_inter[:-1]
					   r2 = r_inter[1:]
					   t1 = theta_coords[:-1]
					   t2 = theta_coords[1:]
					   dS = np.zeros(N-1)
					   dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
					   path_coords[:,1] *= 180./np.pi # convert to deg
					   r_coords = 6370-depth
					   model_vel = np.zeros(N-1) # velocity model for the ray path
					   model_vel_REF = np.zeros(N-1) # velocity model for the ray path
					   for i in range(N-1): # fill the velocity model
					       r_mean = (path_coords[i,0] + path_coords[i+1,0])/2 # take the middle of the cell
					       if r_mean < 6371.0 and r_mean >= 3480.: #  mantle
					           index = np.argmin(np.abs(REF_mantle_r-r_mean))
					           if arrival.name in p_phases:
					           	   REF_vel = REF_mantle_VP[index]
					           	   sample_pt = path_coords[i,:]
					           	   min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
					           	   index = distance.cdist([min_loc],model2D_coords).argmin()
					           	   # model_vel[i] = REF_vel*(1.0+model2D_VP[index]/100.) 
					           	   # NEW VERSION JUNE 2017
					           	   # S. Durand
					           	   model_vel[i] = (model2D_VP[index]/100.)
					           	   model_vel_REF[i] = REF_vel
					           elif arrival.name in s_phases:
					           	   REF_vel = REF_mantle_VS[index]
					           	   sample_pt = path_coords[i,:]
					           	   min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
					           	   index = distance.cdist([min_loc],model2D_coords).argmin()
					           	   # model_vel[i] = REF_vel*(1.0+model2D_VS[index]/100.)
					           	   # NEW VERSION JUNE 2017
					           	   # S. Durand
					           	   model_vel[i] = (model2D_VS[index]/100.)
					           	   model_vel_REF[i] = REF_vel
					           elif arrival.name == 'SKS' or arrival.name == 'SKKS':
					           	   REF_vel = REF_mantle_VS[index]
					           	   sample_pt = path_coords[i,:]
					           	   min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
					           	   index = distance.cdist([min_loc],model2D_coords).argmin()
					           	   # model_vel[i] = REF_vel*(1.0+model2D_VS[index]/100.)
					           	   # NEW VERSION JUNE 2017
					           	   # S. Durand
					           	   model_vel[i] = (model2D_VS[index]/100.)
					           	   model_vel_REF[i] = REF_vel
					       elif r_mean < 3480.0 and r_mean > 1221.5:
					       	   index = np.argmin(np.abs(REF_outercore_r-r_mean))
					       	   # NEW VERSION JUNE 2017
					       	   # S. Durand
					       	   model_vel[i] = 0.
					       	   if arrival.name in p_phases: # outercore
					       	       model_vel_REF[i] = REF_outercore_VP[index]
					       	   elif arrival.name in s_phases:
					       	   	   model_vel_REF[i] = REF_outercore_VS[index]
					       	   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
					       	   	   model_vel_REF[i] = REF_outercore_VP[index]
					       elif r_mean <= 1221.5 and r_mean >= 0: # innercore
					           # NEW VERSION JUNE 2017
					           # S. Durand
					           model_vel[i] = 0.
					           index = np.argmin(np.abs(REF_innercore_r-r_mean))
					           if arrival.name in p_phases:
					           	   model_vel_REF[i] = REF_innercore_VP[index]
					           elif arrival.name in s_phases:
					           	   model_vel_REF[i] = REF_innercore_VS[index]

				# travel time for the 2D model
				for i in range(len(dS)):
					if model_vel_REF[i] !=0:
						# tt2D = tt2D + dS[i]/model_vel[i]
						tt2D = tt2D - dS[i]*model_vel[i]/model_vel_REF[i]
						ttREF = ttREF + dS[i]/model_vel_REF[i]
				dtt2D = abs(theor_tt-PREM_tt)

				path = DIR2 + "/input_files"
				dirs = os.listdir(path)
				os.chdir(DIR2 + '/input_files')
				for i in range(len(dirs)):
					os.remove(dirs[i])
				os.chdir(cwd)

				return tt2D, dtt2D, ttREF, theor_tt, phase_name[0]

def get_travel_time_fromfile(file,elat,elon,edepth,slat,slon,PHASES,PATH_INTERPOLATION = False):

				if (elon)>180:
					elon = elon-360
				if elon<-180:
					elon=elon+360
				if (slon)>180:
					slon = slon-360
				if slon<-180:
					slon=slon+360

				# if os.path.exists(modelpath):
				# 	taup.taup_create.build_taup_model(modelpath,output_folder=DIR2 + '/input_files')
				# 	modelpath2 = glob.glob(DIR2 + '/input_files/*.npz')
				# 	model2 = TauPyModel(model=modelpath2[0])
				# else:
				# 	model2 = TauPyModel(model=modelpath)	

				center_lat2,center_lon2 = midpoint(elat,elon,slat,slon)
				great_circle_dist, azi2, baz2 = gps2dist_azimuth(int(round(elat)),int(round(elon)),int(round(slat)),int(round(slon)))
				great_circle_dist /=  111194.9
				KM_PER_DEG = 111.1949
				width2, az, baz = gps2dist_azimuth(int(round(center_lat2)),int(round(center_lon2)),int(round(slat)),int(round(slon)))

				model_name = 'prem' # select the 1D model
				model = TauPyModel(model='Taup_models/prem/prem_Obspy_all_interpolated.npz')

				phase_name = []

				shutil.copy2(os.path.abspath(file),DIR2 + '/output_files_cross/path.xyz')
				file_new = DIR2 + '/output_files_cross/path.xyz'

				model2D = np.loadtxt(file_new)

				p_phases = ['p','P','Pn','Pdiff','PKP','PKiKP','PKIKP','PcP','pP','pPdiff','pPKP','pPKiKP','PKKP','PKIKKIKP','PP','PKPPKP','PKIKPPKIKP','pPKIKP']
				s_phases = ['s','S','Sn','Sdiff','sS','sSdiff','ScS','SS']
				not_added_converted_phases = ['SKIKS','SKKS','SKIKKIKS','SKSSKS','SKIKSSKIKS','SKiKP','SKP','SKIKP','SKKP','SKIKKIKP','PKS','PKIKS','PKKS','PKIKKIKS','PcS','ScP','SP','PS','sSKS','sSKIKS','pS','pSdiff','pSKS','pSKIKS','sP','sPdiff','sPKP','sPKIKP','sPKiKP']
				added_converted_phases = ['SKS','SKKS']

				if model1 == 'SEISGLOB2':
					model_choice1 = 'prem'
				elif model1 == 'S40RTS':
					model_choice1 = 'prem'
				elif model1 == 'SEMUCBWM1':
					model_choice1 = 'SEMUCBWM1'
				elif model1 == 'S362WMANIM':
					model_choice1 = 'S362WMANIM'
				elif model1 == 'SEISGLOB1':
					model_choice1 = 'prem'
				elif model1 == 'SP12RTS':
					model_choice1 = 'prem'
				elif model1 == 'SGLOBE':
					model_choice1 = 'prem'

				# load the interpolated PREM model
				PREM = np.loadtxt('Taup_models/prem/prem_interpolated.txt')
				PREM_depth = PREM[:,0] # depth in km
				PREM_VP = PREM[:,1] # km/s
				PREM_VS = PREM[:,2] # km/s
				PREM_RHO = PREM[:,3]
				PREM_r = 6371.0 - PREM_depth # radius in km
				# load the interpolated PREM mantle  model
				PREM_mantle = np.loadtxt('Taup_models/prem/prem_mantle_interpolated.txt')
				PREM_mantle_depth = PREM_mantle[:,0] # depth in km
				PREM_mantle_VP = PREM_mantle[:,1] # km/s
				PREM_mantle_VS = PREM_mantle[:,2] # km/s
				PREM_mantle_RHO = PREM_mantle[:,3]
				PREM_mantle_r = 6371.0 - PREM_mantle_depth # radius in km
				# load the interpolated PREM outercore  model
				PREM_outercore = np.loadtxt('Taup_models/prem/prem_outercore_interpolated.txt')
				PREM_outercore_depth = PREM_outercore[:,0] # depth in km
				PREM_outercore_VP = PREM_outercore[:,1] # km/s
				PREM_outercore_VS = PREM_outercore[:,2] # km/s
				PREM_outercore_RHO = PREM_outercore[:,3]
				PREM_outercore_r = 6371.0 - PREM_outercore_depth # radius in km
				# load the interpolated PREM innercore model
				PREM_innercore = np.loadtxt('Taup_models/prem/prem_innercore_interpolated.txt')
				PREM_innercore_depth = PREM_innercore[:,0] # depth in km
				PREM_innercore_VP = PREM_innercore[:,1] # km/s
				PREM_innercore_VS = PREM_innercore[:,2] # km/s
				PREM_innercore_RHO = PREM_innercore[:,3]
				PREM_innercore_r = 6371.0 - PREM_innercore_depth # radius in km

				# load the interpolated PREM model
				REF = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_interpolated.txt')
				REF_depth = REF[:,0] # depth in km
				REF_VP = REF[:,1] # km/s
				REF_VS = REF[:,2] # km/s
				REF_RHO = REF[:,3]
				REF_r = 6371.0 - REF_depth # radius in km
				# load the interpolated REF mantle  model
				REF_mantle = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_mantle_interpolated.txt')
				REF_mantle_depth = REF_mantle[:,0] # depth in km
				REF_mantle_VP = REF_mantle[:,1] # km/s
				REF_mantle_VS = REF_mantle[:,2] # km/s
				REF_mantle_RHO = REF_mantle[:,3]
				REF_mantle_r = 6371.0 - REF_mantle_depth # radius in km
				# load the interpolated REF outercore  model
				REF_outercore = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_outercore_interpolated.txt')
				REF_outercore_depth = REF_outercore[:,0] # depth in km
				REF_outercore_VP = REF_outercore[:,1] # km/s
				REF_outercore_VS = REF_outercore[:,2] # km/s
				REF_outercore_RHO = REF_outercore[:,3]
				REF_outercore_r = 6371.0 - REF_outercore_depth # radius in km
				# load the interpolated REF innercore model
				REF_innercore = np.loadtxt('Taup_models/'  + model_choice1 +'/'+ model_choice1  + '_innercore_interpolated.txt')
				REF_innercore_depth = REF_innercore[:,0] # depth in km
				REF_innercore_VP = REF_innercore[:,1] # km/s
				REF_innercore_VS = REF_innercore[:,2] # km/s
				REF_innercore_RHO = REF_innercore[:,3]
				REF_innercore_r = 6371.0 - REF_innercore_depth # radius in km

				MODEL_2D = True

				# compute the distance from the event to the station
				ev_distance = locations2degrees(elat,elon,slat,slon)
				arrivals = model.get_ray_paths(source_depth_in_km=edepth,distance_in_degree=ev_distance,phase_list=PHASES)
				tt2D = 0.
				dtt2D = 0.
				ttREF = 0.

				arrival = arrivals[0]

				phase_name.append(arrival.name)
				theor_tt = arrival.time # theoretical travel time in seconds
				depth = arrival.path['depth'] # km
				r = 6371.0 - depth # radius - km
				r_deg = r * KM_PER_DEG # radius in deg
				theta_rad = arrival.path['dist'] # distance in radians
				theta_deg = theta_rad*180./np.pi # distance in degrees

				if PATH_INTERPOLATION: # interpolate r and theta given by taup to refine depending the given 2D model
				    N = 10000
				    theta_coords = np.linspace(theta_rad.min(), theta_rad.max(), N)
				    r_inter = np.interp(theta_coords, theta_rad, r)
				    r1 = r_inter[:-1]
				    r2 = r_inter[1:]
				    t1 = theta_coords[:-1]
				    t2 = theta_coords[1:]
				    dS = np.zeros(N-1)
				    dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
				    model_vel = np.zeros(N-1) # velocity model for the ray path
				    r_coords = 6370-depth
				    for i in range(N-1): # fill the velocity model
				        r_mean = (r_inter[i] + r_inter[i+1])/2 # take the middle of the cell
				        if r_mean < 6371.0 and r_mean > 3480.0: #  mantle
				            index = np.argmin(np.abs(PREM_mantle_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_mantle_VP[index] 
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_mantle_VS[index]
				            elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				            	   model_vel[i] = PREM_mantle_VS[index]
				        elif r_mean <= 3480.0 and r_mean > 1221.5:
				        	   index = np.argmin(np.abs(PREM_outercore_r-r_mean))
				        	   if arrival.name in p_phases: # outercore
				        	       model_vel[i] = PREM_outercore_VP[index] 
				        	   elif arrival.name in s_phases:
				        	   	   model_vel[i] = PREM_outercore_VS[index]
				        	   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				        	   	   model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values 
				        elif r_mean <= 1221.5 and r_mean >= 0: # innercore
				            index = np.argmin(np.abs(PREM_innercore_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_innercore_VP[index] 
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_innercore_VS[index] 
				else: # no path interpolation
				    r1 = r[:-1]
				    r2 = r[1:]
				    rtot = r1+r2
				    t1 = theta_rad[:-1]
				    t2 = theta_rad[1:]
				    dS = np.zeros(len(r)-1)
				    dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
				    model_vel = np.zeros(len(r)-1) # velocity model for the ray path
				    for i in range(len(r)-1): # fill the velocity model
				        r_mean = (r[i] + r[i+1])/2 # take the middle of the cell
				        if r_mean < 6371.0 and r_mean >= 3480.: #  mantle
				            index = np.argmin(np.abs(PREM_mantle_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_mantle_VP[index] # no interpolation of PREM values
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_mantle_VS[index] # no interpolation of PREM values
				            elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				            	   model_vel[i] = PREM_mantle_VS[index] # no interpolation of PREM values
				        elif r_mean < 3480.0 and r_mean > 1221.5:
				        	   index = np.argmin(np.abs(PREM_outercore_r-r_mean))
				        	   if arrival.name in p_phases: # outercore
				        	       model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values
				        	   elif arrival.name in s_phases:
				        	   	   model_vel[i] = PREM_outercore_VS[index] # no interpolation of PREM values
				        	   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
				        	   	   model_vel[i] = PREM_outercore_VP[index] # no interpolation of PREM values
				        elif r_mean <= 1221.5 and r_mean >= 0: # innercore
				            index = np.argmin(np.abs(PREM_innercore_r-r_mean))
				            if arrival.name in p_phases:
				            	   model_vel[i] = PREM_innercore_VP[index] # no interpolation of PREM values
				            elif arrival.name in s_phases:
				            	   model_vel[i] = PREM_innercore_VS[index] # no interpolation of PREM values

				# store the PREM numerical computed travel time
				PREM_tt = 0.
				for i in range(len(dS)):
					if model_vel[i] !=0:
						PREM_tt = PREM_tt + dS[i]/model_vel[i]

				if MODEL_2D:
					   # initialize variables
					   model2D_r = np.zeros([len(model2D),1])
					   model2D_theta = np.zeros([len(model2D),1])
					   model2D_VP = np.zeros([len(model2D),1])
					   model2D_VS = np.zeros([len(model2D),1])
					   model2D_RHO = np.zeros([len(model2D),1])
					   model2D_coords = np.zeros([len(model2D),2])
					   # fill the matrices
					   model2D_r = model2D[:,0] # km
					   model2D_theta = model2D[:,1] # deg
					   model2D_VP = model2D[:,2]
					   model2D_VS = model2D[:,3]
					   model2D_RHO = model2D[:,4]
					   model2D_coords[:,0] = model2D[:,0] # rad 
					   model2D_coords[:,1] = model2D[:,1] # theta
					   # interpolate r and theta given by taupy to refine depending the given 2D model
					   N = 1500
					   theta_coords = np.linspace(theta_rad.min(), theta_rad.max(), N)
					   r_inter = np.zeros([N,1])
					   r_inter = np.interp(theta_coords, theta_rad, r)
					   path_coords = np.zeros([N,2])
					   path_coords[:,0] = r_inter
					   path_coords[:,1] = theta_coords
					   r1 = r_inter[:-1]
					   r2 = r_inter[1:]
					   t1 = theta_coords[:-1]
					   t2 = theta_coords[1:]
					   dS = np.zeros(N-1)
					   dS = np.sqrt(r1**2 + r2**2 - 2*r1*r2 * np.cos(t1-t2)) # km
					   path_coords[:,1] *= 180./np.pi # convert to deg
					   r_coords = 6370-depth
					   model_vel = np.zeros(N-1) # velocity model for the ray path
					   model_vel_REF = np.zeros(N-1) # velocity model for the ray path
					   for i in range(N-1): # fill the velocity model
					       r_mean = (path_coords[i,0] + path_coords[i+1,0])/2 # take the middle of the cell
					       if r_mean < 6371.0 and r_mean >= 3480.: #  mantle
					           index = np.argmin(np.abs(REF_mantle_r-r_mean))
					           if arrival.name in p_phases:
					           	   REF_vel = REF_mantle_VP[index]
					           	   sample_pt = path_coords[i,:]
					           	   min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
					           	   index = distance.cdist([min_loc],model2D_coords).argmin()
					           	   # model_vel[i] = REF_vel*(1.0+model2D_VP[index]/100.)
					           	   # NEW VERSION JUNE 2017
					           	   # S. Durand
					           	   model_vel[i] = (model2D_VP[index]/100.)
					           	   model_vel_REF[i] = REF_vel
					           elif arrival.name in s_phases:
					           	   REF_vel = REF_mantle_VS[index]
					           	   sample_pt = path_coords[i,:]
					           	   min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
					           	   index = distance.cdist([min_loc],model2D_coords).argmin()
					           	   # model_vel[i] = REF_vel*(1.0+model2D_VS[index]/100.)
					           	   # NEW VERSION JUNE 2017
					           	   # S. Durand
					           	   model_vel[i] = (model2D_VS[index]/100.)
					           	   model_vel_REF[i] = REF_vel
					           elif arrival.name == 'SKS' or arrival.name == 'SKKS':
					           	   REF_vel = REF_mantle_VS[index]
					           	   sample_pt = path_coords[i,:]
					           	   min_loc = model2D_coords[distance.cdist([sample_pt],model2D_coords).argmin()]
					           	   index = distance.cdist([min_loc],model2D_coords).argmin()
					           	   # model_vel[i] = REF_vel*(1.0+model2D_VS[index]/100.)
					           	   # NEW VERSION JUNE 2017
					           	   # S. Durand
					           	   model_vel[i] = (model2D_VS[index]/100.)
					           	   model_vel_REF[i] = REF_vel
					       elif r_mean < 3480.0 and r_mean > 1221.5:
					       	   index = np.argmin(np.abs(REF_outercore_r-r_mean))
					       	   # NEW VERSION JUNE 2017
					       	   # S. Durand
					       	   model_vel[i] = 0.
					       	   if arrival.name in p_phases: # outercore
					       	       model_vel_REF[i] = REF_outercore_VP[index]
					       	   elif arrival.name in s_phases:
					       	   	   model_vel_REF[i] = REF_outercore_VS[index]
					       	   elif arrival.name == 'SKS' or arrival.name == 'SKKS':
					       	   	   model_vel_REF[i] = REF_outercore_VP[index]
					       elif r_mean <= 1221.5 and r_mean >= 0: # innercore
					           # NEW VERSION JUNE 2017
					           # S. Durand
					           model_vel[i] = 0.
					           index = np.argmin(np.abs(REF_innercore_r-r_mean))
					           if arrival.name in p_phases:
					           	   model_vel_REF[i] = REF_innercore_VP[index]
					           elif arrival.name in s_phases:
					           	   model_vel_REF[i] = REF_innercore_VS[index]

				# travel time for the 2D model
				for i in range(len(dS)):
					if model_vel_REF[i] !=0:
						# tt2D = tt2D + dS[i]/model_vel[i]
						tt2D = tt2D - dS[i]*model_vel[i]/model_vel_REF[i]
						ttREF = ttREF + dS[i]/model_vel_REF[i]
				dtt2D = abs(theor_tt-PREM_tt)

				path = DIR2 + "/input_files"
				dirs = os.listdir(path)
				os.chdir(DIR2 + '/input_files')
				for i in range(len(dirs)):
					os.remove(dirs[i])
				os.chdir(cwd)				

				return tt2D, dtt2D, ttREF, theor_tt, phase_name[0]

def path_tracer_time(lon,lat,rng,az):

    rng = np.radians(rng)/2.

    lat = np.radians(lat)
    lon = np.radians(lon)
    az = np.radians(az)
    phi0 = lat
    lambda0 = lon

    epsilon = 1E-7    # Set tolerance
    if phi0 >= np.pi/2-epsilon:
        az = np.pi
    elif phi0 <= epsilon-np.pi/2:
        az = 0

    # Calculate coordinates of the ending point
    newlat2 = np.arcsin( np.sin(phi0)*np.cos(rng) + np.cos(phi0)*np.sin(rng)*np.cos(az) )
    A = np.sin(rng)*np.sin(az)
    B = np.cos(phi0)*np.cos(rng) - np.sin(phi0)*np.sin(rng)*np.cos(az)
    newlon2 = lambda0 + np.arctan2(A,B)
    ############################################################################
    lon2 = newlon2
    lon2 = np.pi*((np.absolute(lon2)/np.pi) - 2*np.ceil(((np.absolute(lon2)/np.pi)-1)/2)) * np.sign(lon2)
    newlon2 = lon2
    ############################################################################
    latout2 = np.degrees(newlat2) 
    lonout2 = np.degrees(newlon2)

    return [latout2,lonout2]

def path_tracer(lon,lat,rng,az):

    rng = np.radians(rng)

    lat = np.radians(lat)
    lon = np.radians(lon)
    az = np.radians(az)
    phi0 = lat
    lambda0 = lon

    epsilon = 1E-7    # Set tolerance
    if phi0 >= np.pi/2-epsilon:
        az = np.pi
    elif phi0 <= epsilon-np.pi/2:
        az = 0

    az2 = az+np.pi

    # Calculate coordinates of the starting point
    newlat = np.arcsin( np.sin(phi0)*np.cos(rng) + np.cos(phi0)*np.sin(rng)*np.cos(az) )
    A = np.sin(rng)*np.sin(az)
    B = np.cos(phi0)*np.cos(rng) - np.sin(phi0)*np.sin(rng)*np.cos(az)
    newlon = lambda0 + np.arctan2( A,B)
    ############################################################################
    lon = newlon
    lon = np.pi*((np.absolute(lon)/np.pi) - 2*np.ceil(((np.absolute(lon)/np.pi)-1)/2)) * np.sign(lon)
    newlon = lon
    ############################################################################
    latout = np.degrees(newlat) 
    lonout = np.degrees(newlon)

    # Calculate coordinates of the ending point
    newlat2 = np.arcsin( np.sin(phi0)*np.cos(rng) + np.cos(phi0)*np.sin(rng)*np.cos(az2) )
    A = np.sin(rng)*np.sin(az2)
    B = np.cos(phi0)*np.cos(rng) - np.sin(phi0)*np.sin(rng)*np.cos(az2)
    newlon2 = lambda0 + np.arctan2(A,B)
    ############################################################################
    lon2 = newlon2
    lon2 = np.pi*((np.absolute(lon2)/np.pi) - 2*np.ceil(((np.absolute(lon2)/np.pi)-1)/2)) * np.sign(lon2)
    newlon2 = lon2
    ############################################################################
    latout2 = np.degrees(newlat2) 
    lonout2 = np.degrees(newlon2)

    return [latout,lonout,latout2,lonout2]


