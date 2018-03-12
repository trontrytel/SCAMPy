Building and running scampy
===============

$ cd scampy

generate the simulation specific parameters (accepted keywords: Soares, Bomex, DYCOMS_RF01)
$ python generate_namelist.py Soares

generate the turbulence parameters (accepted keywords: defaults, Soares, DYCOMS_RF01)
$ python generate_paramlist.py defaults

compile the source code by running setup.py
$ CC=mpicc python setup.py build_ext --inplace

serial execution of scampy (both turbulence and case specific parameters need to be passed)
$ python main.py Soares.in paramlist_Soares.in
