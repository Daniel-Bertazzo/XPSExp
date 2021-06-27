.. XPS Experiment class documentation master file, created by
   sphinx-quickstart on Wed Jun  2 15:58:24 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##############
XPSExp module
##############

Module for manipulation of XPS experiment data in python.

.. toctree::
   :maxdepth: 2
   :caption: Contents
   
   Installation <installation>
   Class documentation <XPS_Experiment>
   Example 1 <example1>
   Example 2 <example2>


Information
============

This module offers a python class capable of handling data obtained from a X-ray Photoelectron
Spectroscopy (XPS) experiment. The main features include:

* Reading and storing information from a text file (``.txt``) containing data and metadata about a XPS experiment;
* Plotting and visualizing the data in various ways;
* Processing the data to obtain the sample spectrum;
* Correction of distortion (linear and quadratic) that may occur during the experiment;
* Saving the integrated data (spectrum) in VAMAS file format to be opened by other surface analysis softwares.


Indices and tables
==================

* :ref:`genindex`