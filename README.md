# CrossSection_Calc_INCL
This file is a python2.7 class that can calculate cross sections using INCL/ABLA, 
reads them back in, and creates a textfile as the output. Note that INCL/ABLA is
not provided here, neither is root. Please see the docstring for more information.

This program is for example used in order to calculate cross sections for my 
cosmogenic nuclide calculation. See, e.g.,

Trappitsch and Leya (2016), The Astrophysical Journal
http://iopscience.iop.org/article/10.3847/0004-637X/823/1/12/meta

Please have a look at the docstring in the code, especially the doc string for the 
class should be pretty instructive on how to use the code. Please not that you will
need the files

- incl_gridrun.py
- gdata/*
- cs/

The files in gdata are important inputs for your run. The cs folder is the one where
cross sections will be written to, so just keep it there. 

Contact me if you have any questions, inquiries, etc.
