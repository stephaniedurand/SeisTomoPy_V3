# SeisTomoPy

Fast visualization, comparison and calculations in global tomographic models

# What is SeisTomoPy for ?

SeisTomoPy is a new Python tool that facilitates the use of a suite of tomographic models available to the public, with a single programme. SeisTomoPy provides six tools that allow to visualize tomographic models, compare them and extract information for further scientific purposes. The tool comes with a graphical interface with intuitive buttons and simple parameters but the same information can also be gained by using the Python classes that can be run routinely in Python scripts.

# Requirements

We recommend you to start with the installation of anaconda3:

https://www.anaconda.com/download/#macos

And please check that at the end of the installation of anaconda3 there is a ~/.bash_profile file that looks like:

```
# added by Anaconda3 2018.12 installer
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$(CONDA_REPORT_ERRORS=false '/Users/YOURNAME/anaconda3/bin/conda' shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    \eval "$__conda_setup"
else
    if [ -f "/Users/YOURNAME/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/YOURNAME/anaconda3/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        \export PATH="/Users/YOURNAME/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<
```

Then SeisTomoPy has a number of dependencies listed below.

    gfortran : GNU Fortran (MacPorts gcc48 4.8.5_0) 4.8.5

    Python 2.7, 3.5, 3.6, 3.7

    iPython 5.8.0

    matplotlib 3.0.2

    basemap 1.2.0

    pyproj 1.9.6
  
    proj4 5.0.2

    numpy 1.15.4

    obspy 1.1.0

    pyqt 5.9.2

    scipy 1.1.0

For installing the Python dependencies, please run :

conda install -c conda-forge obspy h5py future requests tornado flake8 pytest mock basemap pyqt pip jsonschema responses pyqtgraph pytest-xdist

If there is any problem with the compilation of fortran source files, we recommand to install fortran using:
https://gcc.gnu.org/wiki/GFortranBinaries

# Installing SeisTomoPy

$ git clone https://github.com/stephaniedurand/SeisTomoPy_V3.git

$ cd SeisTomoPy_V3

$ pip install -e .

# User Manual

A Documentation is provided in SeisTomoPy_V3/Documentation/ directory.

# Tutorial

A tutorial is available as a jupyter notebook. You can run it going into the folder SeisTomoPy_V3/Documentation/TomoPy_notebook and running:

$ jupyter notebook

A web window should open. Then you just need to click on SeisTomoPy_GUI_Notebook.ipynb to have an introduction to the GUI interface or on SeisTomoPy_class_Notebook.ipynb for an introduction to the use of SeisTomoPy class.

# Running the GUI

$ ipython -m SeisTomoPy.gui

The graphical interface should pop-up.

# Copyright

Durand S., R. Abreu, C. Thomas, 2017, SeisTomoPy: Fast visualization, comparison and calculations in global tomographic models, Seis. Res. Lett., 89(2A), 658-667
