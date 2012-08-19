About
=====

These are files used to generate the figure in
Y. M. Chee and P. Purkayastha, "Efficient decoding of permutation codes
obtained from distance preserving maps," IEEE International Symposium on
Information Theory 2012, Boston, MA, U.S.A., 641--645.

#. The subdirectory ``sage`` contains the implementation via Sage.
#. The subdirectory ``matlab`` contains the implementation on MATLAB or Octave.
#. The file ``check_distances.c`` contains the implementation to check whether
the mapping remains a DPM as conjectured at the end of the paper.

Installation
============

#. Sage - copy the ``sage`` directory somewhere::

    $ cd sage
    $ sage -python setup.py build_ext --inplace

   Other requirements are sage-5.2, along with the patches in
   http://trac.sagemath.org/10512 and http://trac.sagemath.org/13340


#. MATLAB - simply copy the ``matlab`` directory somewhere


Execution
---------

#. Sage::

    $ sage run_simulation.sage

#. MATLAB or Octave - from within the matlab directory, run the file
    ``simulation``. Make sure to set ``OCTAVE=1`` in ``simulation.m`` if run
    via Octave.
