/******************************************************************************
 ** Copyright (c) 2013-2014, Intel Corporation                                **
 ** All rights reserved.                                                      **
 **                                                                           **
 ** Redistribution and use in source and binary forms, with or without        **
 ** modification, are permitted provided that the following conditions        **
 ** are met:                                                                  **
 ** 1. Redistributions of source code must retain the above copyright         **
 **    notice, this list of conditions and the following disclaimer.          **
 ** 2. Redistributions in binary form must reproduce the above copyright      **
 **    notice, this list of conditions and the following disclaimer in the    **
 **    documentation and/or other materials provided with the distribution.   **
 ** 3. Neither the name of the copyright holder nor the names of its          **
 **    contributors may be used to endorse or promote products derived        **
 **    from this software without specific prior written permission.          **
 **                                                                           **
 ** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
 ** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
 ** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
 ** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
 ** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
 ** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
 ** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
 ** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
 ** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
 ** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
 ** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
 ******************************************************************************/
/* Christopher Dahnken (Intel Corp.), Hans Pabst (Intel Corp.)
 ******************************************************************************/
#define _GNU_SOURCE
#include<dlfcn.h>
#include<stdio.h>
#include<sched.h>

#pragma offload_attribute (push, target (mic))
void mkl_zgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc){
#ifdef __MIC__
  zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}
#pragma offload_attribute (pop) 

#pragma offload_attribute (push, target (mic))
void mkl_dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc){
#ifdef __MIC__
  dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}
#pragma offload_attribute (pop) 

#pragma offload_attribute (push, target (mic))
void mkl_sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc){
#ifdef __MIC__
  sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}
#pragma offload_attribute (pop) 

