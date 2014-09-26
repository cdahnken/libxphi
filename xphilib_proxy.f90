 ! ******************************************************************************
 ! ** Copyright (c) 2013-2014, Intel Corporation                                **
 ! ** All rights reserved.                                                      **
 ! **                                                                           **
 ! ** Redistribution and use in source and binary forms, with or without        **
 ! ** modification, are permitted provided that the following conditions        **
 ! ** are met:                                                                  **
 ! ** 1. Redistributions of source code must retain the above copyright         **
 ! **    notice, this list of conditions and the following disclaimer.          **
 ! ** 2. Redistributions in binary form must reproduce the above copyright      **
 ! **    notice, this list of conditions and the following disclaimer in the    **
 ! **    documentation and/or other materials provided with the distribution.   **
 ! ** 3. Neither the name of the copyright holder nor the names of its          **
 ! **    contributors may be used to endorse or promote products derived        **
 ! **    from this software without specific prior written permission.          **
 ! **                                                                           **
 ! ** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
 ! ** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
 ! ** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
 ! ** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
 ! ** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
 ! ** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
 ! ** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
 ! ** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
 ! ** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
 ! ** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
 ! ** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
 ! ******************************************************************************
 ! * Christopher Dahnken (Intel Corp.), Hans Pabst (Intel Corp.)
 ! ******************************************************************************

subroutine zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
  use xphilibmod
  implicit none

  character :: transa,transb
  integer :: lda,ldb,ldc,ms,ks,ns
  complex(kind=kind(0.0D0)), dimension(lda,*) :: a
  complex(kind=kind(0.0D0)), dimension(ldb,*) :: b
  complex(kind=kind(0.0D0)), dimension(ldc,*) :: c
  complex(kind=kind(0.0D0)) :: alpha, beta
  
  call xphizgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine zgemm

subroutine dgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
  use xphilibmod
  implicit none

  character :: transa,transb
  integer :: lda,ldb,ldc,ms,ks,ns
  real(kind=kind(0.0D0)), dimension(lda,*) :: a
  real(kind=kind(0.0D0)), dimension(ldb,*) :: b
  real(kind=kind(0.0D0)), dimension(ldc,*) :: c
  real(kind=kind(0.0D0)) :: alpha, beta
  
  call xphidgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine dgemm

subroutine sgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
  use xphilibmod
  implicit none

  character :: transa,transb
  integer :: lda,ldb,ldc,ms,ks,ns
  real(kind=kind(0.0E0)), dimension(lda,*) :: a
  real(kind=kind(0.0E0)), dimension(ldb,*) :: b
  real(kind=kind(0.0E0)), dimension(ldc,*) :: c
  real(kind=kind(0.0E0)) :: alpha, beta
  
  call xphisgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine sgemm




