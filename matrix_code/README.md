# About

These are the files used to generate the data for the figures in

1. Y. M. Chee, H. M. Kiah, P. Purkayastha, and C. Wang, "Importance of
   symbol equity in coded modulation for power line communications",
   preprint, 8 pages.

2. Y. M. Chee, H. M. Kiah, and P. Purkayastha, "Matrix codes and multitone
   frequency shift keying in power line communications", to appear in IEEE
   ISIT 2013, 5 pages.

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
Takes a couple of days!::

    $ cd sage
    $ sage -python setup.py build_ext --inplace
    $ time sage simulations.sage | tee -i -a output.txt

You may modify the files `codes.py`, `codes_linear.py` and/or
`simulations.sage` to comment out certain portions depending on the
specific codes you want to simulate.



# The files

1. `powerline_code.pyx` - has the following three classes for concatenated
   codes (see the examples in `codes_linear.py`):
    1. *PowerlineReedSolomonCode* - this class is for the Reed Solomon
       code. It takes in the constant weight code (as a NonlinearCode
       class), the alphabet size q of the Reed Solomon code, and the
       dimension of the Reed Solomon code. Default length of the Reed
       Solomon code is q-1.
    2. *PowerlineReedSolomonSubcode* - this class is for the Reed Solomon
       subcode which contains all words except the all-one and its
       scalings. It takes in the constant weight code (as a NonlinearCode
       class), the alphabet size q of the Reed Solomon code, and the
       dimension of the Reed Solomon code. Default length of the Reed
       Solomon code is q-1.
       The parameters it takes are the same as the PowerlineReedSolomonCode
       class.
    3. *PowerlineReedSolomonCosetCode* - this class is for the Reed Solomon
       coset. The coset leader is generated from the evaluations of the
       generator polynomial of the RS code of one dimension higher.
       The parameters it takes are the same as the PowerlineReedSolomonCode
       class.

    The class *PowerlineCode* can be used for nonlinear codes (see codes.py
    for examples), while the class *PowerlineLinearCode* can be used for
    large linear codes where it is not desired to generate all the
    codewords of the code.

2. `codes_linear.py` generates the concatenated codes. It contains the
   constant weight codes corresponding to A(9,4,4) and A(13,6,5) and the
   Reed Solomon codes.

3. `codes.py` contains the equitable symbol weight codes, the minimum symbol
   weight codes, and the cosets of Reed Solomon codes, that are used in the
   paper titled "Importance of symbol equity in coded modulation for power
   line communications".

4. `mfsk.py` contains the implementation of the M-ary FSK channel.

5. `noises.py` contains the implementation of AWGN noise and Cyclostationary
   noise.

6. `config.pyx` is a file for holding global variables

6. `simulations.sage` is the main file handling all the simulation parameters.

7. `*.pxd` are header files needed for cython.

