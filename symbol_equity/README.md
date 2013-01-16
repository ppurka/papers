# About

These are the files used to generate the data for the figures in
Y. M. Chee, H. M. Kiah, P. Purkayastha, and C. Wang, "Importance of symbol
equity in coded modulation for power line communications", preprint,
8 pages.

The subdirectory `sage` contains the implementation in Sage.


# License

Copyright (C) 2013. `P. Purkayastha <ppurka _at_ gmail _dot_ com>`,
`H. M. Kiah <kiah0001 _at_ ntu _dot_ edu _dot_ sg>`.

This program is released under GPL-3 or later. If the license is not
mentioned in any version of the file, then the default version is GPL-3.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.


# Installation and Execution

Copy the `sage` directory somewhere. Then compile and run the files.
Takes a couple of days!

    $ cd sage
    $ sage -python setup.py build_ext --inplace
    $ time sage simulations.sage | tee -i -a output.txt


# The files

1. `codes.py` - contains all the codes used to generate the data
2. `config.pyx` - a file to contain global variables
3. `powerline_code.pyx` - the main file, written in cython, containing
   all the functions and classes
4. `setup.py` - used to compile the cython files
5. `simulations.sage` - the main Sage file that is run. This contains
   all the setup for the simulations; for instance, the probabilities,
   number of transmissions, etc.

