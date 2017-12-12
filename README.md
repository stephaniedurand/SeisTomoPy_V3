# SeisTomoPy

Fast visualization, comparison and calculations in global tomographic models

# What is SeisTomoPy for ?

SeisTomoPy is a new Python tool that facilitates the use of a suite of tomographic models available to the public, with a single programme. SeisTomoPy provides six tools that allow to visualize tomographic models, compare them and extract information for further scientific purposes. The tool comes with a graphical interface with intuitive buttons and simple parameters but the same information can also be gained by using the Python classes that can be run routinely in Python scripts.

# Requirements

SeisTomoPy has a number of dependencies listed below.

    gfortran : GNU Fortran (MacPorts gcc48 4.8.5_0) 4.8.5

    Python 2.7, 3.5, 3.6

    iPython 5.3.0

    matplotlib 2.0.2

    numpy 1.11.0

    obspy 1.0.1

    pyqt 5.6.0

    scipy 0.17.0

For installing the Python dependencies, please run :

conda install -c conda-forge matplotlib numpy obspy pyqt scipy pyqtgraph basemap

# Installing SeisTomoPy

$ git clone https://github.com/stephaniedurand/SeisTomoPy.git

$ cd SeisTomoPy

$ pip install -e .

# User Manual

A Documentation is provided in SeisTomoPy/Documentation/ directory.

# Tutorial

A tutorial is available as a jupyter notebook. You can run it going into the folder SeisTomoPy/Documentation/TomoPy_notebook and running:

$ jupyter notebook

A web window should open. Then you just need to click on SeisTomoPy_GUI_Notebook.ipynb to have an introduction to the GUI interface or on SeisTomoPy_class_Notebook.ipynb for an introduction to the use of SeisTomoPy class.

# Running the GUI

$ ipython -m SeisTomoPy.gui

The graphical interface should pop-up.

# Copyright

Durand S., R. Abreu, C. Thomas, 2017, SeisTomoPy: Fast visualization, comparison and calculations in global tomographic models, Seis. Res. Lett., ACCEPTED