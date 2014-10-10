rm *.o *.so *.a *.mod *.modmic

COMMOMOPTIONS="-xAVX -fpp -D_BUILD_LIBRARY_  -openmp -fPIC -heap-arrays -offload-attribute-target=mic"

#source  /opt/intel/composer_xe_2013.2.146/bin/compilervars.sh intel64
#ifort -c /opt/intel/composer_xe_2013.2.146/mkl/include/mkl_dfti.f90
ifort $COMMOMOPTIONS   -c xphilib.f90 
ifort $COMMOMOPTIONS   -c xphilib_proxy.f90 
#icc $COMMOMOPTIONS   -c mkl_proxy.c 
icc -shared -fPIC   mkl_proxy.c  -o libmkl_proxy.so -mkl=parallel
ifort  -shared  $COMMOMOPTIONS   xphilib.o xphilib_proxy.o -L. -lmkl_proxy -o libxphi.so

#ifort $COMMOMOPTIONS  main-dynamic.f90 -L. -lmkl_proxy -mkl=parallel -o main-dynamic.x
