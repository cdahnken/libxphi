#!/bin/bash

XPHIZGEMMPATH=$(cd $(dirname $0); pwd -P)


export LD_LIBRARY_PATH=$XPHIZGEMMPATH:$LD_LIBRARY_PATH
export LD_PRELOAD=$XPHIZGEMMPATH/libxphi.so
$*

