import sys
import os
import os.path
import imp
import inspect
import glob
import math
import subprocess, sys
from obspy.core import AttribDict
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from obspy.taup import TauPyModel
from obspy import taup
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, DictFormatter
from obspy.geodetics.base import gps2dist_azimuth
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import  QFileDialog, QMessageBox
import inspect
from glob import iglob
import shutil
from obspy.geodetics import locations2degrees
from scipy.spatial import distance
import re
import SeisTomoPy as SeisTomoPy
from os.path import expanduser
import geopandas
# import conda

# conda_file_dir = conda.__file__
# conda_dir = conda_file_dir.split('lib')[0]
# proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
# os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap


home = expanduser("~")
DIR2 = home + '/SeisTomoPy_files/'
DIR = os.path.dirname(os.path.abspath(__file__))
cwd = os.getcwd()


def compile_ui_files():

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

    filename=DIR + "/qt_window.ui"
    ui_file = filename
    py_ui_file = DIR + "/qt_window.py"
    if not os.path.exists(py_ui_file) or \
                (os.path.getmtime(ui_file) >= os.path.getmtime(py_ui_file)):
        from PyQt5 import uic
        with open(py_ui_file, 'w') as open_file:
            uic.compileUi(ui_file, open_file)
    try:
        import_name = "qt_window"
        globals()[import_name] = imp.load_source(import_name, py_ui_file)
    except ImportError as e:
        print("Error importing %s" % py_ui_file)
        print(e.message)

class Window(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        # Injected by the compile_ui_files() function.
        self.ui = qt_window.Ui_SeisTomoPy() 
        self.ui.setupUi(self)

# Set image for welcome 

        self.fig_intro1 = self.ui.mapfig_intro1.fig
        self.ax_intro1 = self.fig_intro1.add_axes([-0.01, -0.01, 1.03, 1.01])
        image = mpimg.imread(DIR + "/logo.png")
        self.ax_intro1.imshow(image)
        self.ax_intro1.axis("off")
        self.fig_intro1.canvas.draw()

        self.fig_intro2 = self.ui.mapfig_intro2.fig
        self.ax_intro2 = self.fig_intro2.add_axes([-0.02, -0.03, 1.035, 1.03])
        image = mpimg.imread(DIR + "/welcome.png")
        self.ax_intro2.imshow(image)
        self.ax_intro2.axis("off")
        self.fig_intro2.canvas.draw()

# Set some starting values CROSS

        center_lat = 0
        center_lon = 40
        azi_cross = 30
        vmax_cross = 2
        width_cross = 180
        depth_cross = 2890
        tomo_choice_cross = 0
        para_cross = 0
        NSmax_cross = 40
        self.ui.midpt_latitude.setValue(center_lat)
        self.ui.midpt_longitude.setValue(center_lon)
        self.ui.azimuth.setValue(azi_cross)
        self.ui.Vmax_cross.setValue(vmax_cross)
        self.ui.HSdeg_cross.setValue(NSmax_cross)
        self.ui.depth_cross.setValue(depth_cross)
        self.ui.width_cross.setValue(width_cross)
        self.ui.model_cross.setCurrentIndex(tomo_choice_cross)
        self.ui.parameter_cross.setCurrentIndex(para_cross)
        if tomo_choice_cross == 0:
            model_cross = 'SEISGLOB2'
        elif tomo_choice_cross == 1:
            model_cross = 'S40RTS'
        elif tomo_choice_cross == 2:
            model_cross = 'SEMUCB-WM1'
        elif tomo_choice_cross == 3:
            model_cross = 'S362WMANI+M'
        elif tomo_choice_cross == 4:
            model_cross = 'SEISGLOB1'
        elif tomo_choice_cross == 5:
            model_cross = 'SP12RTS'
        elif tomo_choice_cross == 6:
            model_cross = 'SGLOBE-rani'
        elif tomo_choice_cross == 7:
            model_cross = '3D2016-09S'
        elif tomo_choice_cross == 8:
            model_cross = 'MY MODEL'

        if para_cross == 0:
            para1 = 'Vs'
        elif para_cross == 1:
            para1 = 'Vp'
        elif para_cross == 2:
            para1 = 'Density'

        self.ui.label_16.setText('Model: ' + str(model_cross) + ' and Parameter: ' + str(para1))
        rng = float(self.ui.width_cross.value())/2. # distance in degrees
        r = self.center_pt
        lon, lat = r.longitude, r.latitude
        az = float(self.ui.azimuth.value())

        end_lat,end_lon,end_lat2,end_lon2 = SeisTomoPy.path_tracer(lon,lat,rng,az)
        self.ui.label_4.setText(str(int(end_lat2)))
        self.ui.label_9.setText(str(int(end_lon2)))

        self.plot_map()
        self._plot_path()

# Set some starting values MAP

        depth_map = 2800
        NSmax_map = 40
        vmax_map = 2
        para_map = 0
        self.ui.Vmax_map.setValue(vmax_map)
        self.ui.parameter_map.setCurrentIndex(para_map)
        self.ui.depth_map.setValue(depth_map)
        self.ui.HSdeg_map.setValue(NSmax_map)
        self.ui.model_map.setCurrentIndex(4)
        tomo_choice_map = 4
        if tomo_choice_map == 0:
            model_map = 'SEISGLOB2'
        elif tomo_choice_map == 1:
            model_map = 'S40RTS'
        elif tomo_choice_map == 2:
            model_map = 'SEMUCB-WM1'
        elif tomo_choice_map == 3:
            model_map = 'S362WMANI+M'
        elif tomo_choice_map == 4:
            model_map = 'SEISGLOB1'
        elif tomo_choice_map == 5:
            model_map = 'SP12RTS'
        elif tomo_choice_map == 6:
            model_map = 'SGLOBE-rani'
        elif tomo_choice_map == 7:
            model_map = '3D2016'
        elif tomo_choice_map == 8:
            model_map = 'MY MODEL'

        if para_map == 0:
            para1 = 'Vs'
        elif para_map == 1:
            para1 = 'Vp'
        elif para_map == 2:
            para1 = 'Density'

        self.ui.label_20.setText('Model: ' + str(model_map) + ' and Parameter: ' + str(para1))

# Set some starting values SPEC

        dep_spec = 400
        self.ui.depth_spec.setValue(dep_spec)

# Set some starting values CORR

        dep_corr1 = 100
        dep_corr2 = 100
        tomo_choice_corr1 = 0
        tomo_choice_corr2 = 1
        para_corr1 = 0
        para_corr2 = 0
        self.ui.depth1_corr.setValue(dep_corr1)
        self.ui.depth2_corr.setValue(dep_corr2)
        self.ui.model1_corr.setCurrentIndex(tomo_choice_corr1)
        self.ui.model2_corr.setCurrentIndex(tomo_choice_corr2)
        self.ui.parameter1_corr.setCurrentIndex(para_corr1)
        self.ui.parameter2_corr.setCurrentIndex(para_corr2)

################################################################
#                      TOOLS FOR CROSS SECTIONS
################################################################

    def plot_map(self):
        PC=np.loadtxt(DIR + '/hotspots.xy')
        PB=np.loadtxt(DIR + '/plate_boundaries.xy')


        self.fig2_cross = self.ui.mapfig2_cross.fig

        self.ax2_cross = self.fig2_cross.add_axes([0.01, 0.01, .98, .98])

        self.map = Basemap(projection='moll', lon_0=0, resolution="c",ax=self.ax2_cross)

        self.map.drawmapboundary(fill_color='#cccccc')
        self.map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)

        self.fig2_cross.patch.set_alpha(0.0)

        self.fig2_cross.canvas.mpl_connect('button_press_event', self._on_map_mouse_click_event)
        # self.fig2_cross.canvas.mpl_connect('scroll_event',self._zoom_fun)

        xPC, yPC = self.map(PC[:,0],PC[:,1])
        xPB, yPB = self.map(PB[:,0],PB[:,1])
        self.__points_chauds = self.map.scatter(xPC, yPC, s=10, zorder=10,
                                                   color="magenta", marker="o",
                                                   edgecolor="k")
        colors=["7fff00"]
        self.__plate_boundaries = self.map.scatter(xPB, yPB, s=2, zorder=10,
                                                   color="0.3", marker=".")

        EQ=np.loadtxt(DIR + '/catalogue.xy')
        xEQ, yEQ = self.map(EQ[:,0],EQ[:,1])
        self.__eq_map_obj = self.map.scatter(xEQ, yEQ, s=10, zorder=10,
                                                       color="white", marker="o",
                                                       edgecolor="k")

        self.fig2_cross.canvas.draw()

    def _zoom_fun(self, event):
        # get the current x and y limits
        base_scale = 2.
        cur_xlim = self.ax2_cross.get_xlim()
        cur_ylim = self.ax2_cross.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
        cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1./base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1.
        # set new limits
        self.ax2_cross.set_xlim([xdata - cur_xrange*scale_factor,
                     xdata + cur_xrange*scale_factor])
        self.ax2_cross.set_ylim([ydata - cur_yrange*scale_factor,
                     ydata + cur_yrange*scale_factor])
        self.fig2_cross.canvas.draw() # force re-draw


    def _on_map_mouse_click_event(self, event):
        if None in (event.xdata, event.ydata):
            return
        # Get map coordinates by the inverse transform.
        lng, lat = self.map(event.xdata, event.ydata, inverse=True)
        # Left click: set mid point
        if event.button == 1:
            self.ui.midpt_longitude.setValue(lng)
            self.ui.midpt_latitude.setValue(lat)
            rng = float(self.ui.width_cross.value())/2. # distance in degrees
            r = self.center_pt
            lon, lat = r.longitude, r.latitude
            az = float(self.ui.azimuth.value())
            end_lat,end_lon,end_lat2,end_lon2 = SeisTomoPy.path_tracer(lon,lat,rng,az)
            self.ui.label_4.setText(str(int(end_lat2)))
            self.ui.label_9.setText(str(int(end_lon2)))
            self._plot_path()

    def _plot_path(self):
        PC=np.loadtxt(DIR + '/hotspots.xy')

        rng = float(self.ui.width_cross.value())/2. # distance in degrees
        r = self.center_pt
        lon, lat = r.longitude, r.latitude
        az = float(self.ui.azimuth.value())

        try:
            if self.__midpt_map_obj.longitude == lon and \
                    self.__midpt_map_obj.latitude == lat and\
                    self.__azimuth_map_obj == az and\
                    self.__width_map_obj == rng*2.:
                return
        except AttributeError:
            pass

        try:
            self.__midpt_map_obj.remove()
            self.__startpoint_map_obj.remove()
            self.__endpoint_map_obj.remove()
            self.__great_circle_obj.remove()            
            self.__great_circle_obj1.remove()
            self.__PC_map_obj4.remove()
            self.__PC_map_obj5.remove()
            self.__great_circle_obj_path.remove()
            self.__great_circle_obj_path1.remove()
        except AttributeError:
            pass

        x1, y1 = self.map(lon, lat)
        self.__midpt_map_obj = self.map.scatter(x1, y1, s=40, zorder=10,
                                                   color="yellow", marker="o",
                                                   edgecolor="k")
        self.__midpt_map_obj.longitude = lon
        self.__midpt_map_obj.latitude = lat
        self.__azimuth_map_obj = az
        self.__width_map_obj = rng*2.

        #  plot great circle
        end_lat,end_lon,end_lat2,end_lon2 = SeisTomoPy.path_tracer(lon,lat,rng,az)
        x2, y2 = self.map(end_lon,end_lat)
        x3, y3 = self.map(end_lon2,end_lat2)
        self.__great_circle_obj, = self.map.drawgreatcircle(lon,lat,end_lon,end_lat,ls='None',marker='o',markersize=2,color='k')
        self.__great_circle_obj1, = self.map.drawgreatcircle(lon,lat,end_lon2,end_lat2,ls='None',marker='o',markersize=2,color='k')

        # plot starting and end points
        self.__endpoint_map_obj = self.map.scatter(x2, y2, s=40, zorder=10,
                                                   color="red", marker="o",
                                                   edgecolor="k")
        self.__startpoint_map_obj = self.map.scatter(x3, y3, s=40, zorder=10,
                                                   color="lime", marker="o",
                                                   edgecolor="k")
        self.fig2_cross.canvas.draw()

        self.ui.label_4.setText(str(int(end_lat2)))
        self.ui.label_9.setText(str(int(end_lon2)))

    def _plot_eq(self):

        resample = self.ui.earthquake.checkState()
        if resample:
            EQ=np.loadtxt(DIR + '/catalogue.xy')
            xEQ, yEQ = self.map(EQ[:,0],EQ[:,1])
            self.__eq_map_obj = self.map.scatter(xEQ, yEQ, s=20, zorder=10,
                                                       color="white", marker="o",
                                                       edgecolor="k")
            self.fig2_cross.canvas.draw()
        else:
            self.__eq_map_obj.remove()
            self.fig2_cross.canvas.draw()

    def _plot_hotspots(self):

        pcval = self.ui.hotspots.checkState()
        if pcval:
            PC=np.loadtxt(DIR + '/hotspots.xy')
            xPC, yPC = self.map(PC[:,0],PC[:,1])
            self.__points_chauds = self.map.scatter(xPC, yPC, s=20, zorder=10,
                                                   color="magenta", marker="o",
                                                   edgecolor="k")
            self.fig2_cross.canvas.draw()
        else:
            self.__points_chauds.remove()
            self.fig2_cross.canvas.draw()

    def on_plot_cross_released(self):
        lati = int(self.ui.midpt_latitude.value())
        longi = int(self.ui.midpt_longitude.value())
        azi_cross = int(self.ui.azimuth.value())
        vmax_cross = int(self.ui.Vmax_cross.value())
        dep_cross = int(self.ui.depth_cross.value())
        width_cross = int(self.ui.width_cross.value())
        tomo_choice_cross = int(self.ui.model_cross.currentIndex())+1
        NSmax_cross = int(self.ui.HSdeg_cross.value())

        if tomo_choice_cross-1 == 0:
            model_cross = 'SEISGLOB2'
        elif tomo_choice_cross-1 == 1:
            model_cross = 'S40RTS'
        elif tomo_choice_cross-1 == 2:
            model_cross = 'SEMUCBWM1'
        elif tomo_choice_cross-1 == 3:
            model_cross = 'S362WMANIM'
        elif tomo_choice_cross-1 == 4:
            model_cross = 'SEISGLOB1'
        elif tomo_choice_cross-1 == 5:
            model_cross = 'SP12RTS'
        elif tomo_choice_cross-1 == 6:
            model_cross = 'SGLOBE'
        elif tomo_choice_cross-1 == 7:
            model_cross = '3D2016'
        elif tomo_choice_cross-1 == 8:
            model_cross = 'MYMODEL'

        para_cross = int(self.ui.parameter_cross.currentIndex())
        if para_cross == 0:
            para = 'VS'
        elif para_cross == 1:
            para = 'VP'
        elif para_cross == 2:
            para = 'RHO'

        Z, th, r = SeisTomoPy.cross_section_midpoint(model_cross,para,lati,longi,azi_cross,dep_cross,width_cross,NSmax_cross)

        try:
            self.fig1_cross.clf()
        except AttributeError:
            pass

        model = Z[:,2]
        model_new = np.reshape(model, (len(th),len(r)))
        model_new2=np.transpose(model_new)
        model_new2 = np.fliplr(model_new2)
        X, Y = np.meshgrid(th, r)
        self.fig1_cross = self.ui.mapfig1_cross.fig

        th0, th1 = (90-(width_cross/2), 90+(width_cross/2))
        r0, r1 = (6371-dep_cross, 6281)
        thstep, rstep = (10, 500)

        # scale degrees to radians:
        tr_scale = Affine2D().scale(np.pi/180., 1.)
        tr = tr_scale + PolarAxes.PolarTransform()
        theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
        theta_tick_formatter = angle_helper.FormatterDMS()
        r_ticks = [(5871.,r"$500$"),(5371.,r"$1000$"),(4871.,r"$1500$"),\
         (4371,r"$2000$"),(3891,r"$2500$")]
        r_grid_locator  = FixedLocator([v for v, s in r_ticks])
        r_tick_formatter = DictFormatter(dict(r_ticks))
        grid_helper = GridHelperCurveLinear(tr,
                                            extremes=(th0, th1, r0, r1),
                                            grid_locator1=theta_grid_locator,
                                            grid_locator2=r_grid_locator,
                                            tick_formatter1=None,
                                            tick_formatter2=r_tick_formatter)

        self.ax1_cross = FloatingSubplot(self.fig1_cross, 111, grid_helper=grid_helper)
        self.fig1_cross.add_subplot(self.ax1_cross)

        # adjust x axis (theta):
        self.ax1_cross.axis["bottom"].set_visible(True)
        self.ax1_cross.axis["bottom"].toggle(all=False)
        self.ax1_cross.axis["top"].set_visible(True)
        self.ax1_cross.axis["top"].toggle(all=False)
        # adjust y axis (r):
        self.ax1_cross.axis["left"].set_visible(True)
        self.ax1_cross.axis["left"].set_axis_direction("bottom") # tick direction
        self.ax1_cross.axis["left"].major_ticklabels.set_axis_direction("right")

        # create a parasite axes whose transData is theta, r:
        self.auxa = self.ax1_cross.get_aux_axes(tr)
        # make aux_ax to have a clip path as in a?:
        self.auxa.patch = self.ax1_cross.patch 
        # this has a side effect that the patch is drawn twice, and possibly over some other
        # artists. So, we decrease the zorder a bit to prevent this:
        self.ax1_cross.patch.zorder = -2
        self.__s = self.auxa.pcolormesh(X, Y, model_new2, shading='gouraud', cmap='RdYlBu', vmin=-vmax_cross, vmax=vmax_cross)
        eqval = self.ui.earthquake.checkState()
        EQ1=np.loadtxt(DIR2+ 'output_files_cross/eq_profile.xy')
        pcval = self.ui.hotspots.checkState()
        PC3=np.loadtxt(DIR2+'output_files_cross/hotspots_profile.xy')

        if bool(eqval):
        	if np.size(EQ1)/2 != 1:
        		self.__eq = self.auxa.scatter(EQ1[:,0],EQ1[:,1], s=30, zorder=3,\
                                               color="white", marker="o",\
                                               edgecolor="k",clip_on=False)
        if bool(pcval):
        	if np.size(PC3)/2 != 1:
        		self.__points_chauds2 = self.auxa.scatter(PC3[:,0],PC3[:,1], s=30, zorder=3,\
                                               color="magenta", marker="o",\
                                               edgecolor="k",clip_on=False)

        PL=np.loadtxt(DIR2+'output_files_cross/dots_profile.xy')
        self.__points_ligne = self.auxa.scatter(float(PL[0,0]),6281., s=80,zorder=3,\
                                               color="red", marker="o",\
                                               edgecolor="k",clip_on=False)
        self.__points_ligne2 = self.auxa.scatter(float(PL[1,0]),6281., s=80, zorder=3,\
                                       color="yellow", marker="o",\
                                       edgecolor="k",clip_on=False)
        self.__points_ligne3 = self.auxa.scatter(float(PL[2,0]),6281., s=80, zorder=3,\
                                       color="lime", marker="o",\
                                       edgecolor="k",clip_on=False)


        self.fig1_cross.subplots_adjust(bottom=0.2, right=0.9, left=0.1, top=0.95)
        cb = self.fig1_cross.add_axes([0.2, 0.1, 0.6, 0.05]) 
        step2 = float(vmax_cross)/4
        self.__cbar = self.fig1_cross.colorbar(self.__s,orientation="horizontal", cax = cb, ticks=np.arange(-vmax_cross,vmax_cross+step2,step2))
        self.auxa.add_artist(plt.Circle([0, 0], radius=5961., ls='--', lw=1, color='k', fill=False,
                                transform=self.ax1_cross.transData._b, zorder=2))
        self.auxa.add_artist(plt.Circle([0, 0], radius=5701., ls='--', lw=1, color='k', fill=False,
                                transform=self.ax1_cross.transData._b, zorder=2))

        self.fig1_cross.canvas.draw()

    def on_savefig_cross_released(self):
        lati = int(self.ui.midpt_latitude.value())
        longi = int(self.ui.midpt_longitude.value())
        azi_cross = int(self.ui.azimuth.value())
        selected_directory = QFileDialog.getExistingDirectory()
        tomo_choice_cross = int(self.ui.model_cross.currentIndex())
        if tomo_choice_cross == 0:
            model_cross = 'SEISGLOB2'
        elif tomo_choice_cross == 1:
            model_cross = 'S40RTS'
        elif tomo_choice_cross == 2:
            model_cross = 'SEMUCBWM1'
        elif tomo_choice_cross == 3:
            model_cross = 'S362WMANIM'
        elif tomo_choice_cross == 4:
            model_cross = 'SEISGLOB1'
        elif tomo_choice_cross == 5:
            model_cross = 'SP12RTS'
        elif tomo_choice_cross == 6:
            model_cross = 'SGLOBE'
        elif tomo_choice_cross == 7:
            model_cross = '3D2016'
        elif tomo_choice_cross == 8:
            model_cross = 'MYMODEL'

        para1 = int(self.ui.parameter_cross.currentIndex())
        if para1 == 0:
            para = 'VS'
        elif para1 == 1:
            para = 'VP'
        elif para1 == 2:
            para = 'RHO'

        if not selected_directory:
            return

        if os.path.isfile(str(selected_directory) + '/crossec_' + str(model_cross) \
            + '_' + str(para)  + '_'  + str(lati) + '_' + str(longi) + '_' + str(azi_cross) + '.pdf'):

            QMessageBox.critical(None, "Message", "This cross section figure already exists in this directory")
        else:
            self.fig1_cross.savefig(str(selected_directory) + '/crossec_' + str(model_cross) \
            + '_' + str(para)  + '_'  + str(lati) + '_' + str(longi) + '_' + str(azi_cross) + '.pdf',bbox_inches='tight')

        if os.path.isfile(str(selected_directory) + '/map_'+ str(model_cross) \
            + '_' + str(para)  + '_'  + str(lati) + '_' + str(longi) + '_' + str(azi_cross) + '.pdf'):

            QMessageBox.critical(None, "Message", "This map figure already exists in this directory")
        else:
            self.fig2_cross.savefig(str(selected_directory) + '/map_'+ str(model_cross) \
            + '_' + str(para)  + '_'  + str(lati) + '_' + str(longi) + '_' + str(azi_cross) + '.pdf',bbox_inches='tight')

    def on_saveoutput_cross_released(self):
        lati = int(self.ui.midpt_latitude.value())
        longi = int(self.ui.midpt_longitude.value())
        azi_cross = int(self.ui.azimuth.value())
        selected_directory2 = QFileDialog.getExistingDirectory()

        if not selected_directory2:
            return

        if os.path.isdir(str(selected_directory2) + '/output_crosssection_' +\
         str(lati) + '_' + str(longi) + '_' + str(azi_cross)):

            QMessageBox.critical(None, "Message", "This directory already exists in this directory")
        else:
            os.mkdir(str(selected_directory2) + '/output_crosssection_' + str(lati) + '_' + str(longi) + '_' + str(azi_cross))
            files_to_copy1 = set(glob.glob(DIR2 + 'output_files_cross/output_cross_*'))
            for _files in files_to_copy1:
                shutil.copyfile( _files, \
                    str(selected_directory2) + '/output_crosssection_' + str(lati) + '_' + str(longi) + '_' + str(azi_cross)+ '/' + os.path.splitext(os.path.basename(_files))[0]+'.out')
            files_to_copy2 = set(glob.glob(DIR2 + 'output_files_cross/*.sph'))
            for _files in files_to_copy2:
                shutil.copyfile( _files, \
                    str(selected_directory2) + '/output_crosssection_' + str(lati) + '_' + str(longi) + '_' + str(azi_cross) + '/' + os.path.splitext(os.path.basename(_files))[0] +'.sph')
                
            path = DIR2 + "/output_files_cross"
            dirs = os.listdir(path)
            os.chdir(DIR2 + '/output_files_cross')
            for i in range(len(dirs)):
                os.remove(dirs[i])
            os.chdir(DIR)

################################################################
#                      TOOLS FOR MAP
################################################################

    def on_savefig_map_released(self):

        dep_map = int(self.ui.depth_map.value())
        selected_directory4 = QFileDialog.getExistingDirectory()
        tomo_choice_map = int(self.ui.model_map.currentIndex())
        if tomo_choice_map == 0:
            model_map = 'SEISGLOB2'
        elif tomo_choice_map == 1:
            model_map = 'S40RTS'
        elif tomo_choice_map == 2:
            model_map = 'SEMUCBWM1'
        elif tomo_choice_map == 3:
            model_map = 'S362WMANIM'
        elif tomo_choice_map == 4:
            model_map = 'SEISGLOB1'
        elif tomo_choice_map == 5:
            model_map = 'SP12RTS'
        elif tomo_choice_map == 6:
            model_map = 'SGLOBE'
        elif tomo_choice_map == 7:
            model_map = '3D2016'
        elif tomo_choice_map == 8:
            model_map = 'YOUR MODEL'

        para1 = int(self.ui.parameter_map.currentIndex())
        if para1 == 0:
            para = 'VS'
        elif para1 == 1:
            para = 'VP'
        elif para1 == 2:
            para = 'RHO'

        if not selected_directory4:
            return

        if os.path.isfile(str(selected_directory4) + '/map_' + str(model_map) \
            + '_' + str(para)  + '_'  + str(dep_map) + 'km.pdf'):

            QMessageBox.critical(None, "Message", "This figure already exists in this directory")
        else:
            self.fig_map.savefig(str(selected_directory4) + '/map_' + str(model_map) \
            + '_' + str(para)  + '_'  + str(dep_map) + 'km.pdf',bbox_inches='tight')

    def on_saveoutput_map_released(self):

        dep_map = int(self.ui.depth_map.value())
        selected_directory5 = QFileDialog.getExistingDirectory()
        tomo_choice_map = int(self.ui.model_map.currentIndex())
        if tomo_choice_map == 0:
            model_map = 'SEISGLOB2'
        elif tomo_choice_map == 1:
            model_map = 'S40RTS'
        elif tomo_choice_map == 2:
            model_map = 'SEMUCBWM1'
        elif tomo_choice_map == 3:
            model_map = 'S362WMANIM'
        elif tomo_choice_map == 4:
            model_map = 'SEISGLOB1'
        elif tomo_choice_map == 5:
            model_map = 'SP12RTS'
        elif tomo_choice_map == 6:
            model_map = 'SGLOBE'
        elif tomo_choice_map == 7:
            model_map = '3D2016'
        elif tomo_choice_map == 8:
            model_map = 'MYMODEL'

        para1 = int(self.ui.parameter_map.currentIndex())
        if para1 == 0:
            para = 'VS'
        elif para1 == 1:
            para = 'VP'
        elif para1 == 2:
            para = 'RHO'

        if not selected_directory5:
            return

        if os.path.isdir(str(selected_directory5) + '/output_map_' + \
             str(dep_map) + 'km'):

            QMessageBox.critical(None, "Message", "This directory already exists in this directory")
        else:
            shutil.copytree(DIR2 + '/output_files_map', str(selected_directory5) + '/output_map_' + \
             str(dep_map) + 'km')

            path = DIR2 + "/output_files_map"
            dirs = os.listdir(path)
            os.chdir(DIR2 + '/output_files_map')
            for i in range(len(dirs)):
                os.remove(dirs[i])
            os.chdir(DIR)

    def on_plot_map_released(self):
        vmax_map = int(self.ui.Vmax_map.value())
        dep_map = int(self.ui.depth_map.value())
        tomo_choice_map = int(self.ui.model_map.currentIndex())+1
        NSmax_map = int(self.ui.HSdeg_map.value())
        HP2 = self.ui.hotspots_map.checkState()
        PBs = self.ui.plateboundaries.checkState()
        lon_map = int(self.ui.lon0_map.value())
        if tomo_choice_map-1 == 0:
            model_map = 'SEISGLOB2'
        elif tomo_choice_map-1 == 1:
            model_map = 'S40RTS'
        elif tomo_choice_map-1 == 2:
            model_map = 'SEMUCBWM1'
        elif tomo_choice_map-1 == 3:
            model_map = 'S362WMANIM'
        elif tomo_choice_map-1 == 4:
            model_map = 'SEISGLOB1'
        elif tomo_choice_map-1 == 5:
            model_map = 'SP12RTS'
        elif tomo_choice_map-1 == 6:
            model_map = 'SGLOBE'
        elif tomo_choice_map-1 == 7:
            model_map = '3D2016'
        elif tomo_choice_map-1 == 8:
            model_map = 'MYMODEL'

        para_map = int(self.ui.parameter_map.currentIndex())
        if para_map == 0:
            para = 'VS'
        elif para_map == 1:
            para = 'VP'
        elif para_map == 2:
            para = 'RHO'

        Z2, lat, lon = SeisTomoPy.tomomap_GUI(model_map,para,dep_map,NSmax_map)

        try:
            self.fig_map.clf()
        except AttributeError:
            pass

        self.fig_map  = self.ui.mapfig_map.fig
        model2 = Z2[:,2]
        X2, Y2 = np.meshgrid(lon, lat)
        model_new3 = np.reshape(model2, (len(lon),len(lat)))
        model_new4=np.transpose(model_new3)

        self.ax_map  = self.fig_map.add_axes([0.05, 0.01, .9, .9])
        self.map3 = Basemap(projection='moll', lon_0=lon_map, resolution="c",ax=self.ax_map)
        self.__s2 = self.map3.pcolormesh(X2, Y2, model_new4, cmap='RdYlBu', shading='gouraud',latlon=True, vmin=-vmax_map, vmax=vmax_map)
        self.__coast = self.map3.drawcoastlines(linewidth=1)
        step = float(vmax_map)/4
        self.__cbar = self.fig_map.colorbar(self.__s2,orientation="horizontal", ticks=np.arange(-vmax_map,vmax_map+step,step))

        self.fig_map.canvas.draw()

        self._plot_hotspots_map()
        self._plot_plateboundaries()

    def _plot_hotspots_map(self):

        resample2 = self.ui.hotspots_map.checkState()
        if resample2:
            PC2=np.loadtxt('hotspots.xy')
            xPC2, yPC2 = self.map3(PC2[:,0],PC2[:,1])
            self.__points_chauds3 = self.map3.scatter(xPC2, yPC2, s=20, zorder=10,
                                                   color="magenta", marker="o",
                                                   edgecolor="k")
            self.fig_map.canvas.draw()
        else:
            try:
                self.__points_chauds3.remove()
                self.fig_map.canvas.draw()
            except AttributeError:
                pass

    def _plot_plateboundaries(self):

        resample2 = self.ui.plateboundaries.checkState()
        if resample2:
            PB=np.loadtxt('plate_boundaries.xy')
            xPB, yPB = self.map3(PB[:,0],PB[:,1])
            self.__plate_boundaries2 = self.map3.scatter(xPB, yPB, s=2, zorder=10,
                                                   color="0.3", marker=".")
            self.fig_map.canvas.draw()
        else:
            try:
                self.__plate_boundaries2.remove()
                self.fig_map.canvas.draw()
            except AttributeError:
                pass

################################################################
#                      TOOLS FOR PATH
################################################################

    def on_path_plot_released(self):
        vmax_path = int(self.ui.Vmax_3.value())
        depmin = 6371-int(self.ui.depthmin_path.value())
        depmax = 6371-int(self.ui.depthmax_path.value())
        depstep = (self.ui.depthstep_path.value())
        width_path = int(self.ui.width_path.value())
        widthstep = (self.ui.widthstep_path.value())
        para_path = int(self.ui.parameter_6.currentIndex())
        dmodel_path = (self.ui.parameter_7.currentIndex())

        if para_path == 0:
            para6 = 1
        elif para_path == 1:
            para6 = 0
        elif para_path == 2:
            para6 = 2

        if dmodel_path == 0:
            dmodel = "prem"
            modelpath = TauPyModel(model=dmodel)
        elif dmodel_path == 1:
            dmodel = "iasp91"
            modelpath = TauPyModel(model=dmodel)
        elif dmodel_path == 2:
            dmodel = "ak135"
            modelpath = TauPyModel(model=dmodel)
        elif dmodel_path == 3:  
            dmodel = "pwdk"
            modelpath = TauPyModel(model=dmodel)
        elif dmodel_path == 4:  
            modelfile = DIR2 + '/input_files/model_file.nd'      
            taup.taup_create.build_taup_model(modelfile,output_folder=DIR2 + '/input_files')
            modelpath = TauPyModel(model=DIR2 + '/input_files/model_file.npz')

        try:
            self.fig_path.clf()
        except AttributeError:
            pass

        Z = np.loadtxt(DIR2 + 'input_files/cross_section_path.xyz')
        th = np.arange(90-(width_path/2), 90+(width_path/2)+1, widthstep) #*np.pi/180. # deg
        r = np.arange(depmax, depmin+depstep, depstep)

        model = Z[:,para6+2]

        if (len(th)*len(r))!=len(model):

        	QMessageBox.critical(None, "Message", "Either your depth limits or widths limits are wrong")

        model_new = np.reshape(model, (len(th),len(r)))
        model_new6=np.transpose(model_new)
        model_new6 = np.fliplr(model_new6)
        X, Y = np.meshgrid(th, r)
        self.fig_path = self.ui.mapfig_path.fig

        th0, th1 = (90-(width_path/2), 90+(width_path/2))
        r0, r1 = (depmin, 0)
        # r0, r1 = (depmax, depmin)
        thstep, rstep = (widthstep, depstep)
 
        # scale degrees to radians:
        tr_scale = Affine2D().scale(np.pi/180., 1.)
        tr = tr_scale + PolarAxes.PolarTransform()
        theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
        theta_tick_formatter = angle_helper.FormatterDMS()
        r_ticks = [(5871.,r"$500$"),(5371.,r"$1000$"),(4871.,r"$1500$"),\
         (4371,r"$2000$"),(3891,r"$2500$"),(3480,r"$CMB$"),(1221,r"$ICB$"),(6371,r"$0$")]
        r_grid_locator  = FixedLocator([v for v, s in r_ticks])
        r_tick_formatter = DictFormatter(dict(r_ticks))
        grid_helper = GridHelperCurveLinear(tr,
                                            extremes=(th0, th1, r0, r1),
                                            grid_locator1=theta_grid_locator,
                                            grid_locator2=r_grid_locator,
                                            tick_formatter1=None,
                                            tick_formatter2=r_tick_formatter)

        self.ax_path = FloatingSubplot(self.fig_path, 111, grid_helper=grid_helper)
        self.fig_path.add_subplot(self.ax_path)

        # adjust x axis (theta):
        self.ax_path.axis["bottom"].set_visible(True)
        self.ax_path.axis["bottom"].toggle(all=False)
        self.ax_path.axis["top"].set_visible(True)
        self.ax_path.axis["top"].toggle(all=False)
        # adjust y axis (r):
        self.ax_path.axis["left"].set_visible(True)
        self.ax_path.axis["left"].set_axis_direction("bottom") # tick direction
        self.ax_path.axis["left"].major_ticklabels.set_axis_direction("right")

        # create a parasite axes whose transData is theta, r:
        self.auxa = self.ax_path.get_aux_axes(tr)
        # make aux_ax to have a clip path as in a?:
        self.auxa.patch = self.ax_path.patch 
        # this has a side effect that the patch is drawn twice, and possibly over some other
        # artists. So, we decrease the zorder a bit to prevent this:
        self.ax_path.patch.zorder = -2
        self.__s = self.auxa.pcolormesh(X, Y, model_new6, shading='gouraud', cmap='RdYlBu', vmin=-vmax_path, vmax=vmax_path)

        self.fig_path.subplots_adjust(bottom=0.2, right=0.9, left=0.1, top=0.9)
        cb = self.fig_path.add_axes([0.2, 0.07, 0.6, 0.05]) 
        step2 = float(vmax_path)/4
        self.__cbar = self.fig_path.colorbar(self.__s,orientation="horizontal", cax = cb, ticks=np.arange(-vmax_path,vmax_path+step2,step2))
        self.auxa.add_artist(plt.Circle([0, 0], radius=1221, ls='--', lw=1, color='k', fill=False,
                                transform=self.ax_path.transData._b, zorder=2))
        self.auxa.add_artist(plt.Circle([0, 0], radius=3480., ls='--', lw=1, color='k', fill=False,
                                transform=self.ax_path.transData._b, zorder=2))
        # # self.auxa.add_artist(plt.Circle([0, 0], radius=5371., ls='--', lw=1, color='k', fill=False,
        # #                         transform=self.ax1_cross.transData._b, zorder=2))
        stafile = DIR2 + '/input_files/station_file_path.xy'
        evtfile = DIR2 + '/input_files/event_file_path.xy'
        if os.path.exists(stafile):
            STA = np.loadtxt(DIR2 + '/input_files/station_file_path.xy')
            if STA.ndim == 1:
                self.__station = self.auxa.scatter(-(width_path-180)/2.+width_path-(STA[0]),6371.-(STA[1]), s=80,zorder=3,\
                                               color="black", marker="^",\
                                               edgecolor="k",clip_on=False)
                lensta = 1
            else:
                self.__station = self.auxa.scatter(-(width_path-180)/2.+width_path-(STA[:,0]),6371.-(STA[:,1]), s=80,zorder=3,\
                                               color="black", marker="^",\
                                               edgecolor="k",clip_on=False)
                lensta = len(STA)
            os.remove(DIR2 + '/input_files/station_file_path.xy')
        if os.path.exists(evtfile):
            EVT = np.loadtxt(DIR2 + '/input_files/event_file_path.xy')
            if EVT.ndim == 1:
                lenevt = 1
                self.__event = self.auxa.scatter(-(width_path-180)/2.+width_path-(EVT[0]),6371.-(EVT[1]), s=80,zorder=3,\
                                               color="red", marker="*",\
                                               edgecolor="k",clip_on=False)
            else:
                lenevt = len(EVT)
                self.__event = self.auxa.scatter(-(width_path-180)/2.+width_path-(EVT[:,0]),6371.-(EVT[:,1]), s=80,zorder=3,\
                                               color="red", marker="*",\
                                               edgecolor="k",clip_on=False)
            os.remove(DIR2 + '/input_files/event_file_path.xy')
        os.remove(DIR2 + '/input_files/cross_section_path.xyz')

        if str(self.ui.lineEdit.text()):
            liste = str(self.ui.lineEdit.text())
            liste_ph = liste.split()
            for i in range(lenevt):
                for j in range(lensta):
                    if lenevt == 1:
                        if lensta ==1:
                            arrivals = modelpath.get_ray_paths(source_depth_in_km=EVT[1],distance_in_degree=abs(STA[0]-EVT[0]),phase_list=eval(str(liste_ph)))
                        else:
                            arrivals = modelpath.get_ray_paths(source_depth_in_km=EVT[1],distance_in_degree=abs(STA[j,0]-EVT[0]),phase_list=eval(str(liste_ph)))
                    else:
                        if lensta ==1:
                            arrivals = modelpath.get_ray_paths(source_depth_in_km=EVT[i,1],distance_in_degree=abs(STA[0]-EVT[i,0]),phase_list=eval(str(liste_ph)))
                        else:
                            arrivals = modelpath.get_ray_paths(source_depth_in_km=EVT[i,1],distance_in_degree=abs(STA[j,0]-EVT[i,0]),phase_list=eval(str(liste_ph)))
                    nber_phase=np.size(arrivals)
                    for k in range(nber_phase):
                        arrival = arrivals[k]
                        path_depth = arrival.path['depth']
                        if lenevt == 1:
                            path_deg = arrival.path['dist']*180./np.pi+EVT[0]
                            self.__path = self.auxa.plot(-(width_path-180)/2.+width_path-path_deg,6371.-path_depth,linewidth=1.0,color="black",linestyle="-",clip_on=False)
                        else:
                            path_deg = arrival.path['dist']*180./np.pi+EVT[i,0]
                            self.__path = self.auxa.plot(-(width_path-180)/2.+width_path-path_deg,6371.-path_depth,linewidth=1.0,color="black",linestyle="-",clip_on=False)

        path = DIR2 + "/input_files"
        dirs = os.listdir(path)
        os.chdir(DIR2 + '/input_files')
        for i in range(len(dirs)):
            os.remove(dirs[i])
        os.chdir(DIR)

        self.fig_path.canvas.draw()
        self.ui.label_45.setText('')
        self.ui.label_48.setText('')
        self.ui.label_49.setText('')
        self.ui.label_58.setText('')
        self.ui.lineEdit.setText('')

    def on_savefig_path_released(self):

        selected_directory6 = QFileDialog.getExistingDirectory()

        if not selected_directory6:
            return

        if os.path.isfile(str(selected_directory6) + '/path.pdf'):
            QMessageBox.critical(None, "Message", "This figure already exists in this directory")
        else:
            self.fig_path.savefig(str(selected_directory6) + '/path.pdf',bbox_inches='tight')

    def on_loadfile_path_released(self):
        pwd = os.getcwd()
        self.corr_file3, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.sph", pwd,"Cross-section file (*.sph);;"))
        if not self.corr_file3:
            return
        if self.corr_file3.endswith('.sph'):
            self.ui.label_45.setText(os.path.basename(self.corr_file3))
            filename7 = DIR2 + '/input_files/cross_section_path.xyz'
            shutil.copy2(self.corr_file3,filename7)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file3.split('.')[-1])

    def on_loadfile_eq_released(self):
        pwd = os.getcwd()
        self.corr_file4, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.xy", pwd,"Event file (*.xy);;"))
        if not self.corr_file4:
            return
        if self.corr_file4.endswith('.xy'):
            self.ui.label_48.setText(os.path.basename(self.corr_file4))
            filename7 = DIR2 + '/input_files/event_file_path.xy'
            shutil.copy2(self.corr_file4,filename7)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file4.split('.')[-1])

    def on_loadfile_st_released(self):
        pwd = os.getcwd()
        self.corr_file5, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.xy", pwd,"Station file (*.xy);;"))
        if not self.corr_file5:
            return
        if self.corr_file5.endswith('.xy'):
            self.ui.label_49.setText(os.path.basename(self.corr_file5))
            filename7 = DIR2 + '/input_files/station_file_path.xy'
            shutil.copy2(self.corr_file5,filename7)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file5.split('.')[-1])


    def on_loadfile_model_released(self):
        pwd = os.getcwd()
        self.corr_file6, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.nd", pwd,"Model file (*.nd);;"))
        if not self.corr_file6:
            return
        if self.corr_file6.endswith('.nd'):
            self.ui.label_58.setText(os.path.basename(self.corr_file6))
            filename7 = DIR2 + '/input_files/model_file.nd'
            shutil.copy2(self.corr_file6,filename7)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file6.split('.')[-1])

################################################################
#                      TOOLS FOR TIME
################################################################

    def on_calculate_time_released(self):
        # self.progressBar.start()
        tomo_choice_time = int(self.ui.model_time.currentIndex())+1

        if tomo_choice_time-1 == 0:
            model_path = 'SEISGLOB2'
        elif tomo_choice_time-1 == 1:
            model_path = 'S40RTS'
        elif tomo_choice_time-1 == 2:
            model_path = 'SEMUCBWM1'
        elif tomo_choice_time-1 == 3:
            model_path = 'S362WMANIM'
        elif tomo_choice_time-1 == 4:
            model_path = 'SEISGLOB1'
        elif tomo_choice_time-1 == 5:
            model_path = 'SP12RTS'
        elif tomo_choice_time-1 == 6:
            model_path = 'SGLOBE'

        # dmodel_path = (self.ui.parameter_8.currentIndex())
        # if dmodel_path == 0:
        #     dmodel = "prem"
        #     modelpath = TauPyModel(model=dmodel)
        # elif dmodel_path == 1:
        #     dmodel = "iasp91"
        #     modelpath = TauPyModel(model=dmodel)
        # elif dmodel_path == 2:
        #     dmodel = "ak135"
        #     modelpath = TauPyModel(model=dmodel)
        # elif dmodel_path == 3:  
        #     dmodel = "pwdk"
        #     modelpath = TauPyModel(model=dmodel)
        # elif dmodel_path == 4:  
        #     dmodel = DIR2 + '/input_files/model_file.nd'

        G = np.loadtxt(DIR2 + 'output_files_time/timepy_data_file.xy')
        if G.ndim ==1:
            ELAT = G[0]
            ELON = G[1]
            EDEPTH = G[2]
            SLAT = G[3]
            SLON = G[4]
            stalen = 1
            evtlen = 1
        else:
            ELAT = G[:,0]
            ELON = G[:,1]
            EDEPTH = G[:,2]
            SLAT = G[:,3]
            SLON = G[:,4]
            stalen = len(ELAT)
            evtlen = len(ELAT)

        if str(self.ui.lineEdit_4.text()):
        	liste = str(self.ui.lineEdit_4.text())
        	liste_ph = liste.split()
        	PHASES = eval(str(liste_ph))
        else:
        	PHASES = ['ttall']

        sol_filename1 = DIR2 + '/output_files_time/travel_times.out'
        solution_file1 = open(sol_filename1, 'w') # open travel time file
        file_str = 'elat & elon & edepth & slat & slon & Seismic phase & 2D travel time (s) & Travel time precision (s) & Theoretical travel time (s) \n'
        solution_file1.write(file_str)

        big_string=''
        model_name = 'prem' # select the 1D model
        model = TauPyModel(model=DIR + '/../Taup_models/prem/prem_Obspy_all_interpolated.npz')

        p_phases = ['p','P','Pn','Pdiff','PKP','PKiKP','PKIKP','PcP','pP','pPdiff','pPKP','pPKiKP','PKKP','PKIKKIKP',\
        'PP','PKPPKP','PKIKPPKIKP','pPKIKP']
        s_phases = ['s','S','Sn','Sdiff','sS','sSdiff','ScS','SS']
        not_added_converted_phases = ['SKIKS','SKKS','SKIKKIKS','SKSSKS','SKIKSSKIKS','SKiKP','SKP','SKIKP','SKKP'\
        ,'SKIKKIKP','PKS','PKIKS','PKKS','PKIKKIKS','PcS','ScP','SP','PS','sSKS','sSKIKS','pS',\
        'pSdiff','pSKS','pSKIKS','sP','sPdiff','sPKP','sPKIKP','sPKiKP']
        added_converted_phases = ['SKS','SKKS']

        for k in range(stalen):
            if stalen == 1:
                ev_distance = locations2degrees(ELAT,ELON,SLAT,SLON)
            else:
                ev_distance = locations2degrees(ELAT[k],ELON[k],SLAT[k],SLON[k])

            for i_phase in range(len(PHASES)):
                if stalen == 1:
                    arrivals = model.get_ray_paths(source_depth_in_km=EDEPTH,distance_in_degree=ev_distance,phase_list=[PHASES[i_phase]])
                else:
                    arrivals = model.get_ray_paths(source_depth_in_km=EDEPTH[k],distance_in_degree=ev_distance,phase_list=[PHASES[i_phase]])

            	if len(arrivals) == 0:
                    big_string = big_string + PHASES[i_phase] + ' does not exist \n' + \
                    '----------------------------------------------------------------------------'+ '\n'         
            	else:
                    arrival = arrivals[0]
                    if arrival.name != PHASES[i_phase]:
                            big_string = big_string + PHASES[i_phase] + ' does not exist \n' +\
                            '----------------------------------------------------------------------------'+ '\n' 
                    else:
                        if arrival.name not in p_phases and arrival.name not in s_phases and arrival.name not in added_converted_phases:
                            big_string = big_string + 'The phase ' + arrival.name + ' can not be computed yet. Sorry! \n' + \
                            '----------------------------------------------------------------------------'+ '\n' 
                        else:
                            if stalen == 1:
                                tt, dtt, tref, ttheor, phase_name = SeisTomoPy.get_travel_time(model_path,ELAT,ELON,EDEPTH,SLAT,SLON,[PHASES[i_phase]],dmodel,PATH_INTERPOLATION = False)
                            else:
                                tt, dtt, tref, ttheor, phase_name = SeisTomoPy.get_travel_time(model_path,ELAT[k],ELON[k],EDEPTH[k],SLAT[k],SLON[k],[PHASES[i_phase]],PATH_INTERPOLATION = False)
                            val1 = np.round((tt+tref)*10.)/10.
                            val2 = np.round((dtt)*100.)/100.
                            val3 = np.round((ttheor)*100.)/100.
                            if stalen == 1:
                                file_str = str(ELAT) + ' ' + str(ELON) + ' ' +  str(EDEPTH) + ' ' + str(SLAT) + ' ' + str(SLON) + ' ' + str(phase_name) + ' ' + str(val1) + ' ' + str(val2) + ' ' + str(val3) + '\n'
                            else:
                                file_str = str(ELAT[k]) + ' ' + str(ELON[k]) + ' ' +  str(EDEPTH[k]) + ' ' + str(SLAT[k]) + ' ' + str(SLON[k]) + ' ' + str(phase_name) + ' ' + str(val1) + ' ' + str(val2) + ' ' + str(val3) + '\n'
                            solution_file1.write(file_str)

                            big_string = big_string + 'Phase name: '+ phase_name + '\n' + '2D travel time: ' + str(val1)+ ' s' + '\n' + 'Error +/-: ' + str(val2)+ ' s' + '\n' + '----------------------------------------------------------------------------'+ '\n' 

        path = DIR2 + "output_files_cross"
        dirs = os.listdir(path)
        os.chdir(DIR2 + 'output_files_cross')
        for i in range(len(dirs)):
            os.remove(dirs[i])
        os.chdir(DIR)

        solution_file1.close()
        self.ui.textEdit_time.setText(big_string)
        self.ui.lineEdit_4.setText('')
        self.ui.label_115.setText('')
        # self.ui.label_60.setText('')

    def on_saveoutput_time_released(self):

        selected_directory6 = QFileDialog.getExistingDirectory()
        if not selected_directory6:
            return

        if os.path.isdir(str(selected_directory6) + '/travel_times_output'):
            QMessageBox.critical(None, "Message", "This directory already exists")
        else:
            shutil.copytree(DIR2 + 'output_files_time', str(selected_directory6) + '/travel_times_output')
            path = DIR2 + "output_files_time"
            dirs = os.listdir(path)
            os.chdir(DIR2 + 'output_files_time')
            for i in range(len(dirs)):
                os.remove(dirs[i])
            os.chdir(DIR)

    def on_loadfile_time_released(self):
        pwd = os.getcwd()
        self.corr_file6, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.xy", pwd,"Data file (*.xy);;"))
        if not self.corr_file6:
            return
        if self.corr_file6.endswith('.xy'):
            self.ui.label_115.setText(os.path.basename(self.corr_file6))
            filename7 = DIR2 + 'output_files_time/timepy_data_file.xy'
            shutil.copy2(self.corr_file6,filename7)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file6.split('.')[-1])

    def on_loadfile_path_2_released(self):
        pwd = os.getcwd()
        self.corr_file6, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.nd", pwd,"Model file (*.nd);;"))
        if not self.corr_file6:
            return
        if self.corr_file6.endswith('.nd'):
            self.ui.label_60.setText(os.path.basename(self.corr_file6))
            filename7 = DIR2 + '/input_files/model_file.nd'
            shutil.copy2(self.corr_file6,filename7)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file6.split('.')[-1])

################################################################
#                      TOOLS FOR SPECTRUM
################################################################

    def on_plot_spec_released(self):
        dep_spec = int(self.ui.depth_spec.value())
        NSmax_spec = int(self.ui.HSdeg_spec.value())
        para_spec = int(self.ui.parameter_spec.currentIndex())+1
        if para_spec == 1:
            para = 'VS'
        elif para_spec == 2:
            para = 'VP'
        elif para_spec == 3:
            para = 'RHO'

        try:
            self.fig_spec.clf()
        except AttributeError:
            pass

        self.fig_spec  = self.ui.mapfig_spec.fig
        self.ax_spec  = self.fig_spec.add_axes([0.1, 0.15, .6, .8])


        if bool(self.ui.SEISGLOB2.checkState()):
       	    model_spec = 'SEISGLOB2'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp1 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="red",marker="d",markeredgecolor="k", label='SEISGLOB2')
            else:
                sp1 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="red",marker="d", markeredgecolor="k",label='SEISGLOB2')
        if bool(self.ui.SEISGLOB1.checkState()):
            model_spec = 'SEISGLOB1'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp1 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="orange",marker="d",markeredgecolor="k", label='SEISGLOB1')
            else:
                sp1 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="orange",marker="d", markeredgecolor="k",label='SEISGLOB1')
        if bool(self.ui.SP12RTS.checkState()):
            model_spec = 'SP12RTS'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp2 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="navy",marker="d",markeredgecolor="k", label='SP12RTS')
            else:
                sp2 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="navy",marker="d", markeredgecolor="k",label='SP12RTS')
        if bool(self.ui.S40RTS.checkState()):
            model_spec = 'S40RTS'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp2 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="blue",marker="d",markeredgecolor="k", label='S40RTS')
            else:
                sp2 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="blue",marker="d", markeredgecolor="k",label='S40RTS')
        if bool(self.ui.SEMUCB.checkState()):
            model_spec = 'SEMUCBWM1'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp3 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="dodgerblue",marker="d",markeredgecolor="k", label='SEMUCB-WM1')
            else:
                sp3 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="dodgerblue",marker="d",markeredgecolor="k", label='SEMUCB-WM1')
        if bool(self.ui.S362ANI.checkState()):
            model_spec = 'S362WMANIM'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp4 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="cyan",marker="d",markeredgecolor="k", label='S362WMANI+M')
            else:
                sp4 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="cyan",marker="d", markeredgecolor="k",label='S362WMANI+M')
        if bool(self.ui.SGLOBE.checkState()):
            model_spec = 'SGLOBE'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp4 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="k",marker="d",markeredgecolor="k", label='SGLOBE-rani')
            else:
                sp4 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="k",marker="d", markeredgecolor="k",label='SGLOBE-rani')
        if bool(self.ui.D2016.checkState()):
            model_spec = '3D2016'
            sp = SeisTomoPy.spectrum_GUI(model_spec,para,dep_spec,NSmax_spec)
            if bool(self.ui.norm.checkState()):
                sp4 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="m",marker="d",markeredgecolor="k", label='3D2016-09S')
            else:
                sp4 =self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="m",marker="d", markeredgecolor="k",label='3D2016-09S')
        if os.path.isfile(DIR2 + '/input_files/map_file.xyz'):
        	sp = SeisTomoPy.spectrum_fromfile_GUI(DIR2 + '/input_files/map_file.xyz',NSmax_spec)
        	self.ui.label_31.setText('')
        	if bool(self.ui.norm.checkState()):
        		sp4 =self.ax_spec.plot(sp[:,0], sp[:,1]/np.amax(sp[:,1]), linewidth=1.0,color="lime",marker="d", markeredgecolor="k",label='from file')
        	else:
        		sp4 = self.ax_spec.plot(sp[:,0], sp[:,1], linewidth=1.0,color="lime",marker="d",markeredgecolor="k", label='from file')

        path = DIR2 + "/input_files"
        dirs = os.listdir(path)
        os.chdir(DIR2 + '/input_files')
        for i in range(len(dirs)):
            os.remove(dirs[i])
        os.chdir(DIR)

        path = DIR2 + "output_files_map"
        dirs = os.listdir(path)
        os.chdir(DIR2 + 'output_files_map')
        for i in range(len(dirs)):
            os.remove(dirs[i])
        os.chdir(DIR)

        self.ax_spec.legend(bbox_to_anchor=(1.45, 1))
        self.ax_spec.set_xlabel("Harmonic degree l")
        self.ax_spec.set_ylabel("Spectrum amplitude")
        self.ax_spec.set_xlim([0.5, NSmax_spec+0.5])
        self.ax_spec.grid()
        self.fig_spec.canvas.draw()

    def on_savefig_spec_released(self):

        dep_spec = int(self.ui.depth_spec.value())

        selected_directory6 = QFileDialog.getExistingDirectory()

        if not selected_directory6:
            return

        if os.path.isfile(str(selected_directory6) + '/spectre_' + str(dep_spec) \
            + 'km.pdf'):
            QMessageBox.critical(None, "Message", "This figure already exists in this directory")
        else:
            self.fig_spec.savefig(str(selected_directory6) + '/spectre_' + str(dep_spec) \
            + 'km.pdf',bbox_inches='tight')

    def on_saveoutput_spec_released(self):

        dep_spec = int(self.ui.depth_spec.value())

        selected_directory7 = QFileDialog.getExistingDirectory()

        if not selected_directory7:
            return

        if os.path.isdir(str(selected_directory7) + '/output_spectre_' + \
             str(dep_spec) + 'km'):
            QMessageBox.critical(None, "Message", "This directory already exists")
        else:
            shutil.copytree(DIR2 + '/output_files_spectre', str(selected_directory7) + '/output_spectre_' + \
             str(dep_spec) + 'km')

            path = DIR2 + "/output_files_spectre"
            dirs = os.listdir(path)
            os.chdir(DIR2 + '/output_files_spectre')
            for i in range(len(dirs)):
                os.remove(dirs[i])
            os.chdir(DIR)

    def on_loadfile_spec_released(self):
        pwd = os.getcwd()
        self.spectre_file, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.xyz", pwd,"Map file (*.xyz);;"))
        if not self.spectre_file:
            return
        if self.spectre_file.endswith('.xyz'):
            self.ui.label_31.setText(os.path.basename(self.spectre_file))
            filename6 = DIR2 + '/input_files/map_file.xyz'
            shutil.copy2(self.spectre_file,filename6)
        else:
            raise IOError('unknown file type *.%s' %
                          self.spectre_file.split('.')[-1])

################################################################
#                      TOOLS FOR CORRELATION
################################################################

    def on_plot_corr_released(self):
        dep_corr1 = int(self.ui.depth1_corr.value())
        NSmax_corr = int(self.ui.HSdeg_corr.value())
        para_corr1 = int(self.ui.parameter1_corr.currentIndex())+1
        dep_corr2 = int(self.ui.depth2_corr.value())
        para_corr2 = int(self.ui.parameter2_corr.currentIndex())+1
        if para_corr1 == 1:
            para1 = 'VS'
        elif para_corr1 == 2:
            para1 = 'VP'
        elif para_corr1 == 3:
            para1 = 'RHO'

        if para_corr2 == 1:
            para2 = 'VS'
        elif para_corr2 == 2:
            para2 = 'VP'
        elif para_corr2 == 3:
            para2 = 'RHO'

        tomo_choice_corr1 = int(self.ui.model1_corr.currentIndex())+1
        tomo_choice_corr2 = int(self.ui.model2_corr.currentIndex())+1
        if tomo_choice_corr1-1 == 0:
            model_corr1 = 'SEISGLOB2'
        elif tomo_choice_corr1-1 == 1:
            model_corr1 = 'S40RTS'
        elif tomo_choice_corr1-1 == 2:
            model_corr1 = 'SEMUCBWM1'
        elif tomo_choice_corr1-1 == 3:
            model_corr1 = 'S362WMANIM'
        elif tomo_choice_corr1-1 == 4:
            model_corr1 = 'SEISGLOB1'
        elif tomo_choice_corr1-1 == 5:
            model_corr1 = 'SP12RTS'
        elif tomo_choice_corr1-1 == 6:
            model_corr1 = 'SGLOBE'
        elif tomo_choice_corr1-1 == 7:
            model_corr1 = '3D2016'
        elif tomo_choice_corr1-1 == 8:
            model_corr1 = 'MYMODEL'
        elif tomo_choice_corr1-1 == 9:
            model_corr1 = 'None'
        if tomo_choice_corr2-1 == 0:
            model_corr2 = 'SEISGLOB2'
        elif tomo_choice_corr2-1 == 1:
            model_corr2 = 'S40RTS'
        elif tomo_choice_corr2-1 == 2:
            model_corr2 = 'SEMUCBWM1'
        elif tomo_choice_corr2-1 == 3:
            model_corr2 = 'S362WMANIM'
        elif tomo_choice_corr2-1 == 4:
            model_corr2 = 'SEISGLOB1'
        elif tomo_choice_corr2-1 == 5:
            model_corr2 = 'SP12RTS'
        elif tomo_choice_corr2-1 == 6:
            model_corr2 = 'SGLOBE'
        elif tomo_choice_corr2-1 == 7:
            model_corr2 = '3D2016'
        elif tomo_choice_corr2-1 == 8:
            model_corr2 = 'MYMODEL'
        elif tomo_choice_corr2-1 == 9:
            model_corr2 = 'None'

        try:
            self.fig_corr.clf()
        except AttributeError:
            pass

        self.fig_corr  = self.ui.mapfig_corr.fig
        self.ax_corr  = self.fig_corr.add_axes([0.1, 0.15, .5, .8])

        if (model_corr1 != 'None') and (model_corr2 != 'None'):
            sp = SeisTomoPy.correlation_GUI(model_corr1,model_corr2,dep_corr1,dep_corr2,para1,para2,NSmax_corr)
            sp1 =self.ax_corr.plot(sp[:,0], sp[:,1], linewidth=1.0,color="grey",marker="d", label=str(model_corr1) + '-' +str(model_corr2))
        elif ((model_corr1 != 'None') and (model_corr2 == 'None')) or ((model_corr1 == 'None') and (model_corr2 != 'None')):
            if os.path.isfile(DIR2 + '/input_files/map_file1.xyz'):
                self.ui.label_33.setText('')
                sp = SeisTomoPy.correlation_fromfile_GUI(model_corr2,dep_corr2,para2,DIR2 + '/input_files/map_file1.xyz',NSmax_corr)
                sp1 =self.ax_corr.plot(sp[:,0], sp[:,1], linewidth=1.0,color="grey",marker="d", label='from file-' +str(model_corr2))
            elif os.path.isfile(DIR2 + '/input_files/map_file2.xyz'):
                self.ui.label_37.setText('')
                sp = SeisTomoPy.correlation_fromfile_GUI(model_corr1,dep_corr1,para1,DIR2 + '/input_files/map_file2.xyz',NSmax_corr)
                sp1 =self.ax_corr.plot(sp[:,0], sp[:,1], linewidth=1.0,color="grey",marker="d", label='from file-' +str(model_corr1))
        elif (model_corr1 == 'None') and (model_corr2 == 'None'):
            if os.path.isfile(DIR2 + '/input_files/map_file1.xyz') and os.path.isfile(DIR2 + '/input_files/map_file2.xyz'):
                sp = SeisTomoPy.correlation_fromfile2_GUI(DIR2 + '/input_files/map_file1.xyz',DIR2 + '/input_files/map_file2.xyz',NSmax_corr)
                self.ui.label_33.setText('')
                self.ui.label_37.setText('')
                sp1 =self.ax_corr.plot(sp[:,0], sp[:,1], linewidth=1.0,color="grey",marker="d", label='from file-from file')

        path = DIR2 + "/input_files"
        dirs = os.listdir(path)
        os.chdir(DIR2 + '/input_files')
        for i in range(len(dirs)):
            os.remove(dirs[i])
        os.chdir(DIR)

        path = DIR2 + "output_files_map"
        dirs = os.listdir(path)
        os.chdir(DIR2 + 'output_files_map')
        for i in range(len(dirs)):
            os.remove(dirs[i])
        os.chdir(DIR)

        conf66 = np.loadtxt('conf66.dat')
        conf90 = np.loadtxt('conf90.dat')
        conf95 = np.loadtxt('conf95.dat')
        conf =self.ax_corr.plot(conf66[:,0], conf66[:,1], linewidth=1.0,color="silver", ls='--',label='Confidence level 66%')
        conf =self.ax_corr.plot(conf90[:,0], conf90[:,1], linewidth=1.0,color="gray", ls='--',label='Confidence level 90%')
        conf =self.ax_corr.plot(conf95[:,0], conf95[:,1], linewidth=1.0,color="black", ls='--',label='Confidence level 95%')
        conf =self.ax_corr.plot(conf66[:,0], -conf66[:,1], linewidth=1.0,color="silver", ls='--')
        conf =self.ax_corr.plot(conf90[:,0], -conf90[:,1], linewidth=1.0,color="gray", ls='--')
        conf =self.ax_corr.plot(conf95[:,0], -conf95[:,1], linewidth=1.0,color="black", ls='--')
        self.ax_corr.set_xlabel("Harmonic degree l")
        self.ax_corr.set_ylabel("Correlation")
        self.ax_corr.set_xlim([0.5, NSmax_corr+0.5])
        self.ax_corr.set_ylim([-1, 1])
        self.ax_corr.legend(bbox_to_anchor=(1.6, 1))
        self.ax_corr.grid()
        self.fig_corr.canvas.draw()

    def on_savefig_corr_released(self):

        dep_corr1 = int(self.ui.depth1_corr.value())
        dep_corr2 = int(self.ui.depth2_corr.value())

        selected_directory6 = QFileDialog.getExistingDirectory()

        if not selected_directory6:
            return

        if os.path.isfile(str(selected_directory6) + '/correlation.pdf'):
        	QMessageBox.critical(None, "Message", "This figure already exists in this directory")
        else:
        	self.fig_corr.savefig(str(selected_directory6) + '/correlation.pdf',bbox_inches='tight')

    def on_loadfile1_corr_released(self):
        pwd = os.getcwd()
        self.corr_file1, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.xyz", pwd,"Map file (*.xyz);;"))
        if not self.corr_file1:
            return
        if self.corr_file1.endswith('.xyz'):
            self.ui.label_33.setText(os.path.basename(self.corr_file1))
            filename6 = DIR2 + '/input_files/map_file1.xyz'
            shutil.copy2(self.corr_file1,filename6)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file1.split('.')[-1])

    def on_loadfile2_corr_released(self):
        pwd = os.getcwd()
        self.corr_file2, _ = (QFileDialog.getOpenFileName(
            self, "Choose *.xyz", pwd,"Map file (*.xyz);;"))
        if not self.corr_file2:
            return
        if self.corr_file2.endswith('.xyz'):
            self.ui.label_37.setText(os.path.basename(self.corr_file2))
            filename6 = DIR2 + '/input_files/map_file2.xyz'
            shutil.copy2(self.corr_file2,filename6)
        else:
            raise IOError('unknown file type *.%s' %
                          self.corr_file2.split('.')[-1])

    def on_saveoutput_corr_released(self):

        selected_directory8 = QFileDialog.getExistingDirectory()

        if not selected_directory8:
            return

        if os.path.isdir(str(selected_directory8) + '/output_corr'):

        	QMessageBox.critical(None, "Message", "This directory already exists in this directory")
        else:
        	shutil.copytree(DIR2 + '/output_files_corr', str(selected_directory8) + '/output_corr')
           	path = DIR2 + '/output_files_corr'
            dirs = os.listdir(path)
            os.chdir(DIR2 + '/output_files_corr')
            for i in range(len(dirs)):
                os.remove(dirs[i])
            os.chdir(DIR)

################################################################
#                         UPDATES
################################################################

    def update_cross(self, force=False):

        rng = float(self.ui.width_cross.value())/2. # distance in degrees
        r = self.center_pt
        lon, lat = r.longitude, r.latitude
        az = float(self.ui.azimuth.value())

        try:
            end_latf,end_lonf,end_lati,end_loni = SeisTomoPy.path_tracer(lon,lat,rng,az)
            self._plot_path()
            self.ui.label_4.setText(str(int(end_lati)))
            self.ui.label_9.setText(str(int(end_loni)))
            tomo_choice_cross = int(self.ui.model_cross.currentIndex())
            para_cross = int(self.ui.parameter_cross.currentIndex())
            if tomo_choice_cross == 0:
                model_cross = 'SEISGLOB2'
            elif tomo_choice_cross == 1:
                model_cross = 'S40RTS'
            elif tomo_choice_cross == 2:
                model_cross = 'SEMUCB-WM1'
            elif tomo_choice_cross == 3:
                model_cross = 'S362WMANI+M'
            elif tomo_choice_cross == 4:
                model_cross = 'SEISGLOB1'
            elif tomo_choice_cross == 5:
                model_cross = 'SP12RTS'
            elif tomo_choice_cross == 6:
                model_cross = 'SGLOBE-rani'
            elif tomo_choice_cross == 7:
                model_cross = '3D2016-09S'
            elif tomo_choice_cross == 8:
                model_cross = 'YOUR MODEL'


            if para_cross == 0:
                para = 'Vs'
            elif para_cross == 1:
                para = 'Vp'
            elif para_cross == 2:
                para = 'Density'

            self.ui.label_16.setText('Model: ' + str(model_cross) + ' and Parameter: ' + str(para))

        except AttributeError:
            return

    def update_map(self, force=False):

        try:
            tomo_choice_map = int(self.ui.model_map.currentIndex())
            para_map = int(self.ui.parameter_map.currentIndex())
            if tomo_choice_map == 0:
                model_map = 'SEISGLOB2'
            elif tomo_choice_map == 1:
                model_map = 'S40RTS'
            elif tomo_choice_map == 2:
                model_map = 'SEMUCB-WM1'
            elif tomo_choice_map == 3:
                model_map = 'S362WMANI+M'
            elif tomo_choice_map == 4:
                model_map = 'SEISGLOB1'
            elif tomo_choice_map == 5:
                model_map = 'SP12RTS'
            elif tomo_choice_map == 6:
                model_map = 'SGLOBE-rani'
            elif tomo_choice_map == 7:
                model_map = '3D2016-09S'
            elif tomo_choice_map == 8:
                model_map = 'YOUR MODEL'

            if para_map == 0:
                para = 'Vs'
            elif para_map == 1:
                para = 'Vp'
            elif para_map == 2:
                para = 'Density'

            self.ui.label_20.setText('Model: ' + str(model_map) + ' and Parameter: ' + str(para))

        except AttributeError:
            return

    @property
    def center_pt(self):
        center_pt = AttribDict({'latitude': float(self.ui.midpt_latitude.value()),
                'longitude': float(self.ui.midpt_longitude.value())})
        return center_pt

    def on_midpt_latitude_valueChanged(self, *args):
        self.update_cross()

    def on_midpt_longitude_valueChanged(self, *args):
        self.update_cross()

    def on_azimuth_valueChanged(self, *args):
        self.update_cross()

    def on_depth_cross_valueChanged(self, *args):
        self.update_cross()

    def on_width_cross_valueChanged(self, *args):
        self.update_cross()

    def on_model_cross_currentIndexChanged(self):
        self.update_cross()

    def on_parameter_cross_currentIndexChanged(self):
        self.update_cross()

    def on_earthquake_stateChanged(self): 
        self._plot_eq()

    def on_hotspots_stateChanged(self): 
        self._plot_hotspots()

    def on_hotspots_map_stateChanged(self): 
        self._plot_hotspots_map()

    def on_plateboundaries_stateChanged(self): 
        self._plot_plateboundaries()

    def on_model_map_currentIndexChanged(self):
        self.update_map()

    def on_parameter_map_currentIndexChanged(self):
        self.update_map()

    def on_actionDocumentation_triggered(self):
        subprocess.call(["open",DIR  + "/../../Documentation/SeisTomoPy_DOC.pdf"])

compile_ui_files()

app = QtWidgets.QApplication(sys.argv)
# app = QtGui.QApplication(sys.argv, QtGui.QApplication.GuiClient)
# Show and bring window to foreground.
window = Window()
window.show()
app.installEventFilter(window)
window.raise_()
os._exit(app.exec_())
