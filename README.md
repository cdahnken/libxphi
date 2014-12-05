LIBXPHI
=======

The LIBXPHI library allows to collect hotspot functions that are beneficial to be offloaded to Intel Xeon Phi coprocessors (an instance of the Intel Many Integrated Core Architecture "MIC"). The library already comes with support for some BLAS level 3 functions (xGEMM) but is not limited LAPACK or BLAS.

The library leverages the offload usage model while intercepting calls to a function that is thought to run on the host. Calling the coprocessor's code version may be protected by a threshold that determines whether an offload is beneficial or not.

Instructions
============
Make sure you do this using a Shell where you sourced the Intel Compilers as well as the MPI.

```sh
source /opt/intel/composer_xe_2013_sp1.4.211/bin/compilervars.sh intel64
source /opt/intel/impi/5.0.1.035/intel64/bin/mpivars.sh
```

Change into the directory where you cloned LIBXPHI into and run the build script. 

```sh
./build-library.sh
```

You should now find two libraries 'libxphi.so' and 'libmkl_proxy.so'.

MPIRUN WRAPPER
==============
The MPIRUN wrapper script is used to generate and execute an MPIRUN command line. The wrapper can be found under https://github.com/hfp/mpirun and cloned running:

```sh
git clone https://github.com/hfp/mpirun.git
```
