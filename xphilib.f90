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
 ! ******************************************************************************/
 ! * Christopher Dahnken (Intel Corp.), Hans Pabst (Intel Corp.)
 ! ******************************************************************************

! ==============================================
! -- Module xphilibmod
! -- Chris Dahnken
! -- Intel Corporation 2014
! -- 
! -- This module provide the driver routines
! -- for a replacement of ZGEMM and DGEMM
! -- called from a dynamically linked MKL binary
! ==============================================
module xphilibmod
  use ifport, only : dclock
  implicit none
  !dir$ attributes offload:mic                :: abuff1
  !dir$ attributes align: 64                  :: abuff1
  !dir$ attributes offload:mic                :: bbuff1
  !dir$ attributes align: 64                  :: bbuff1
  !dir$ attributes offload:mic                :: abuff2
  !dir$ attributes align: 64                  :: abuff2
  !dir$ attributes offload:mic                :: bbuff2
  !dir$ attributes align: 64                  :: bbuff2
  !dir$ attributes offload:mic                :: cbuff1
  !dir$ attributes align: 64                  :: cbuff1
  !dir$ attributes offload:mic                :: cbuff2
  !dir$ attributes align: 64                  :: cbuff2
  !dir$ attributes offload:mic                :: cbuff3
  !dir$ attributes align: 64                  :: cbuff3

  private
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: abuff1
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: abuff2
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: bbuff1
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: bbuff2
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: cbuff1
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: cbuff2
  complex(kind=kind(0.0D0)), allocatable, dimension(:),target :: cbuff3


  integer :: signalabuff1=101
  integer :: signalabuff2=102
  integer :: signalbbuff1=103
  integer :: signalbbuff2=104
  integer :: signalcbuff1=105
  integer :: signalcbuff2=106
  integer :: signalcbuff3=107
  logical :: isallocated=.false.
  integer :: computesignal=108


  integer :: blocksize_m
  integer :: blocksize_n
  integer :: blocksize_k

  integer :: offload_device
  logical :: offload_verbose
  real(kind=kind(0.0D0)) :: offload_threshold
  
  public :: xphizgemm,xphidgemm,xphigemm_allocate, xphigemm_deallocate !, xphisgemm, xphizhegvx, xphizgemv,xphidgemv

contains
  ! -- this function is just included since the inclusion 
  ! -- of mkl.fi stretches the compile time a factor of 100x
  logical          function lsame( ca, cb )
    character          ca, cb
    intrinsic          ichar
    integer            inta, intb, zcode
    lsame = ca.eq.cb
    if( lsame ) then 
       return
    endif
    zcode = ichar( 'z' )

    inta = ichar( ca )
    intb = ichar( cb )

    if( zcode.eq.90 .or. zcode.eq.122 ) then
       if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
       if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32

    else if( zcode.eq.233 .or. zcode.eq.169 ) then

       if( inta.ge.129 .and. inta.le.137 .or. &
            inta.ge.145 .and. inta.le.153 .or. &
            inta.ge.162 .and. inta.le.169 ) inta = inta + 64
       if( intb.ge.129 .and. intb.le.137 .or. &
            intb.ge.145 .and. intb.le.153 .or. &
            intb.ge.162 .and. intb.le.169 ) intb = intb + 64

    else if( zcode.eq.218 .or. zcode.eq.250 ) then
       if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
       if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
    end if
    lsame = inta.eq.intb
  end function lsame

  integer function numblocks(l,b)
    implicit none
    integer :: l,b
    numblocks=floor(dble(l)/dble(b))+1
    if(mod(l,b).eq.0) then
       numblocks=numblocks-1
    end if
  end function numblocks

  integer function buffersize(bs,nb,i)
    implicit none
    integer :: bs,nb,i
    integer :: ret
    buffersize=min(nb*i,bs)-min(nb*(i-1),bs)
  end function buffersize

  ! -- midx_of_blocknumber
  ! -- returns the m-dimension index of a give block number i
  integer function midx_of_block(l,blocknum_n,blocknum_k)
    implicit none
    integer :: l
    integer :: blocknum_n,blocknum_k
    midx_of_block=floor(dble(l-1)/dble(blocknum_n*blocknum_k))+1
  end function midx_of_block

  ! -- nidx_of_blocknumber
  ! -- returns the n-dimension index of a give block number i
  integer function nidx_of_block(l,blocknum_m,blocknum_k)
    implicit none
    integer :: l
    integer :: blocknum_m,blocknum_k
    nidx_of_block=mod(floor(dble(l-1)/dble(blocknum_k)),blocknum_m)+1
  end function nidx_of_block

  ! -- kidx_of_blocknumber
  ! -- returns the k-dimension index of a give block number i
  integer function kidx_of_block(l,blocknum_k)
    implicit none
    integer :: l
    integer :: blocknum_k
    kidx_of_block=mod(l-1,blocknum_k)+1
  end function kidx_of_block


  ! -- shapelinear
  ! -- shapelinear is a helper function that reshapes parts of the matrix
  ! -- A into the buffer ABUFF depending on the parameter TRANSA
  ! -- if(TRANSA.eq.'N') ABUFF <- A(r1:r2,c1,c2) else ABUFF <- A(c1:c2,r1,r2)
  subroutine shapelinear(transa,abuff,a,r1,r2,c1,c2,buffsize1,buffsize2,ldim)
    implicit none
    character :: transa    
    integer :: r1,r2,c1,c2,buffsize1,buffsize2
    integer :: ldim
    complex(kind=kind(0.0D0)), dimension(:) :: abuff
    complex(kind=kind(0.0D0)), dimension(ldim,*) :: a
    if(lsame(transa,'n')) then 
       abuff(1:buffsize1*buffsize2)=reshape(a(r1:r2,c1:c2),(/buffsize1*buffsize2/))
    else 
       abuff(1:buffsize1*buffsize2)=reshape(a(c1:c2,r1:r2),(/buffsize1*buffsize2/))
    end if
  end subroutine shapelinear

  ! -- shapelinear_d - a double precision version of shapelinear
  ! -- shapelinear is a helper function that reshapes parts of the matrix
  ! -- A into the buffer ABUFF depending on the parameter TRANSA
  ! -- if(TRANSA.eq.'N') ABUFF <- A(r1:r2,c1,c2) else ABUFF <- A(c1:c2,r1,r2)
  subroutine shapelinear_d(transa,abuff,a,r1,r2,c1,c2,buffsize1,buffsize2,ldim)
    implicit none
    character :: transa    
    integer :: r1,r2,c1,c2,buffsize1,buffsize2
    integer :: ldim
    real(kind=kind(0.0D0)), dimension(:) :: abuff
    real(kind=kind(0.0D0)), dimension(ldim,*) :: a
    if(lsame(transa,'n')) then 
       abuff(1:buffsize1*buffsize2)=reshape(a(r1:r2,c1:c2),(/buffsize1*buffsize2/))
    else 
       abuff(1:buffsize1*buffsize2)=reshape(a(c1:c2,r1:r2),(/buffsize1*buffsize2/))
    end if
  end subroutine shapelinear_d

  subroutine shapelinear_s(transa,abuff,a,r1,r2,c1,c2,buffsize1,buffsize2,ldim)
    implicit none
    character :: transa    
    integer :: r1,r2,c1,c2,buffsize1,buffsize2
    integer :: ldim
    real(kind=kind(0.0E0)), dimension(:) :: abuff
    real(kind=kind(0.0E0)), dimension(ldim,*) :: a
    if(lsame(transa,'n')) then 
       abuff(1:buffsize1*buffsize2)=reshape(a(r1:r2,c1:c2),(/buffsize1*buffsize2/))
    else 
       abuff(1:buffsize1*buffsize2)=reshape(a(c1:c2,r1:r2),(/buffsize1*buffsize2/))
    end if
  end subroutine shapelinear_s


  ! -- xphigemm_allocate
  ! -- sets the block sizes blocksize_m, blocksize_n and blocksize_k and 
  ! -- allocates the buffers of size blocksize_m*blocksize_n*blocksize_k.
  ! -- Parameters passed here can be overwritten with environment 
  ! -- variables, if they are set.
  subroutine xphigemm_allocate(ibm,ibn,ibk)
    use ifport
    implicit none
    integer :: ibm,ibn,ibk
    integer :: bm,bn,bk
    double precision :: t1
    character(len=32) :: carg

    ! -- the offload device should default to 0
    offload_device=0
    ! -- offload threshold should default to 25GFLOP (~the perf of a single Xeon core)
    offload_threshold=25.0

    bm=ibm
    bn=ibn
    bk=ibk

    ! -- override the presets and subroutine arguments
    ! -- with the values set in the env variables.
    ! -- WARNING: currently the env variables must be set
    ! -- since I seem not to be a ble to detect when they are
    ! -- defined and when not. 
    ! -- NEEDSFIX
    call getenv("QE_MIC_VERBOSE",carg)
    if(len(trim(carg)).ne.0) then
!       read(carg,*) offload_device
       read(carg,*) offload_verbose
    end if


    call getenv("QE_MIC_DEVICE",carg)
    if(len(trim(carg)).ne.0) then
       read(carg,*) offload_device
    end if

    call getenv("QE_MIC_BLOCKSIZE_M",carg)
    if(len(trim(carg)).ne.0) then
       read(carg,*) bm
    end if

    call getenv("QE_MIC_BLOCKSIZE_N",carg)
    if(len(trim(carg)).ne.0) then
       read(carg,*) bn
    end if

    call getenv("QE_MIC_BLOCKSIZE_K",carg)
    if(len(trim(carg)).ne.0) then
       read(carg,*) bk
    end if

    call getenv("QE_MIC_OFFLOAD_THRESHOLD",carg)
    if(len(trim(carg)).ne.0) then
       read(carg,*) offload_threshold
       offload_threshold=offload_threshold*dble(1000)*dble(1000)*dble(1000)
    end if

    write(*,*) "allocating buffers", bm,bn,bk 
    write(*,*) "on device ", offload_device
    write(*,*) "threshold ", offload_threshold
    t1=dclock()

    ! -- set the signal numbers
    ! -- this could be done better, so that
    ! -- the signals are unique, e.g. with 
    ! -- the address of the respective fields.
    signalabuff1=111
    signalabuff2=222
    signalbbuff1=333
    signalbbuff2=444
    signalcbuff1=555

    ! -- write the values intot he global fields
    blocksize_m=bm
    blocksize_n=bn
    blocksize_k=bk

    ! -- allocate the host-side buffers
    allocate(abuff1(bm*bk))
    allocate(abuff2(bm*bk))
    allocate(bbuff1(bk*bn))
    allocate(bbuff2(bk*bn))
    allocate(cbuff1(bm*bn))
    allocate(cbuff2(bm*bn))
    allocate(cbuff3(bm*bn))

    ! -- allocate all necessary fields on the device. No data transfer needed at this time
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(abuff1:length(bm*bk) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(abuff2:length(bm*bk) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(bbuff1:length(bk*bn) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(bbuff2:length(bk*bn) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(cbuff1:length(bm*bn) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(cbuff2:length(bm*bn) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(cbuff3:length(bm*bn) ALIGN(64) ALLOC_IF(.true.) FREE_IF(.false.))

    write(*,*) "buffer allocation",dclock()-t1,"s" 
    isallocated=.true.

  end subroutine xphigemm_allocate



  subroutine xphigemm_deallocate()
    implicit none
    integer :: bm,bn,bk

    bm=blocksize_m
    bn=blocksize_n
    bk=blocksize_k

    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(cbuff1:length(bm*bn) ALLOC_IF(.false.) FREE_IF(.true.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(cbuff2:length(bm*bn) ALLOC_IF(.false.) FREE_IF(.true.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(cbuff3:length(bm*bn) ALLOC_IF(.false.) FREE_IF(.true.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(abuff1:length(bm*bk) ALLOC_IF(.false.) FREE_IF(.true.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(abuff2:length(bm*bk) ALLOC_IF(.false.) FREE_IF(.true.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(bbuff1:length(bk*bn) ALLOC_IF(.false.) FREE_IF(.true.))
    !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) NOCOPY(bbuff2:length(bk*bn) ALLOC_IF(.false.) FREE_IF(.true.))

    deallocate(abuff1)
    deallocate(bbuff1)
    deallocate(abuff2)
    deallocate(bbuff2)
    deallocate(cbuff1)
    deallocate(cbuff2)
    deallocate(cbuff3)

  end subroutine xphigemm_deallocate


  subroutine xphidgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
    use ifport, only : dclock
    use iso_c_binding
    implicit none
#ifdef _BUILD_LIBRARY_
    !dir$ attributes offload:mic                :: MKL_DGEMM
#else
    !dir$ attributes offload:mic                :: DGEMM
#endif
    character :: transa,transb
    integer :: lda,ldb,ldc,ms,ks,ns,bm,bn,bk
    real(kind=kind(0.0D0)), dimension(lda,*) :: a
    real(kind=kind(0.0D0)), dimension(ldb,*) :: b
    real(kind=kind(0.0D0)), dimension(ldc,*) :: c
    real(kind=kind(0.0D0)) :: alpha, beta
    real(kind=kind(0.0D0)), dimension(:), pointer :: aptrcur_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: aptrnxt_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: bptrcur_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: bptrnxt_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: cptrcur_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: cptrnxt_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: cptrprv_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: ptrtmp_d

    !dir$ attributes offload:mic :: cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, cur_index_m, cur_index_n, cur_index_k, nxt_buffsize_k, nxt_index_k

    integer :: i,j,k,l
    integer :: nxti,nxtj,nxtk
    integer :: prvi,prvj,prvk

    integer :: nkb,nmb,nnb
    integer :: nblocks
    integer :: iblock
    integer :: cur_buffsize_m
    integer :: cur_buffsize_n
    integer :: cur_buffsize_k
    integer :: prv_buffsize_m
    integer :: prv_buffsize_n
    integer :: prv_buffsize_k
    integer :: nxt_buffsize_m
    integer :: nxt_buffsize_n
    integer :: nxt_buffsize_k
    integer :: cur_index_m
    integer :: cur_index_n
    integer :: cur_index_k
    integer :: prv_index_m
    integer :: prv_index_n
    integer :: prv_index_k
    integer :: nxt_index_m
    integer :: nxt_index_n
    integer :: nxt_index_k
    real(kind=kind(0.0D0)) :: ONE=dble(1.0D0)
    real(kind=kind(0.0D0)) :: ZERO=dble(0.0D0)
    integer, dimension(:), allocatable :: m_of_i 
    integer, dimension(:), allocatable :: n_of_i 
    integer, dimension(:), allocatable :: k_of_i 
    integer :: count
    double precision :: t1
    double precision :: tsend
    double precision :: tcomp
    double precision :: trecv
    double precision :: talloc
    double precision :: tdeall
    logical :: catchcomputesignal=.false.

    IF ((ms.EQ.0) .OR. (ns.EQ.0) .OR.(ks.EQ.0)) RETURN
    
    if(.not.isallocated) then
       call xphigemm_allocate(1000,1000,1000)
       isallocated=.true.
    end if

    t1=dclock()
    ! -- check if the matrix computation is actually less than the 
    ! -- threshold performance. In this case, just do it on the host.
    if((dble(2)*dble(ms)*dble(ns)*dble(ks)).lt.(offload_threshold)) then
!       write(*,*) "about to call mkl_zgemm on host"
#ifdef _BUILD_LIBRARY_
       call mkl_dgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#else
       call dgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#endif
    else
       ! -- if the total number of FLOPs in the matrix is bigger than the threshold
       ! -- throw the matrix into the blocking algorithm.
       bm=blocksize_m
       bn=blocksize_n
       bk=blocksize_k

       ! -- compute the total number of blocks that must run though
       nmb=numblocks(ms,bm)
       nnb=numblocks(ns,bn)
       nkb=numblocks(ks,bk)
       nblocks=nmb*nnb*nkb


       allocate(m_of_i(nblocks))
       allocate(n_of_i(nblocks))
       allocate(k_of_i(nblocks))

       count=0
       ! -- this is recording the plan for the computation 
       ! -- it seems more efficient to record and allocate than
       ! -- to do it in the compute loop. Overhead is in the order 
       ! -- of usecs
       do i=1,nmb
          do j=1,nnb
             do l=1,nkb
                count=count+1
                m_of_i(count)=i
                n_of_i(count)=j
                k_of_i(count)=l
             end do
          end do
       end do


       call c_f_pointer(c_loc(abuff1),aptrnxt_d,[blocksize_m*blocksize_k])
       call c_f_pointer(c_loc(abuff2),aptrcur_d,[blocksize_m*blocksize_k])
       call c_f_pointer(c_loc(bbuff1),bptrnxt_d,[blocksize_k*blocksize_n])
       call c_f_pointer(c_loc(bbuff2),bptrcur_d,[blocksize_k*blocksize_n])
       call c_f_pointer(c_loc(cbuff1),cptrnxt_d,[blocksize_m*blocksize_n])
       call c_f_pointer(c_loc(cbuff2),cptrcur_d,[blocksize_m*blocksize_n])
       call c_f_pointer(c_loc(cbuff3),cptrprv_d,[blocksize_m*blocksize_n])
       

       do iblock=0,nblocks+1
          if( (iblock.eq.1).or.(iblock.eq.(nblocks+1)).or.(m_of_i(iblock).ne.m_of_i(iblock-1)).or.(n_of_i(iblock).ne.n_of_i(iblock-1)) ) then
             ptrtmp_d=>cptrnxt_d
             cptrnxt_d=>cptrprv_d
             cptrprv_d=>cptrcur_d
             cptrcur_d=>ptrtmp_d
          end if
          
          ptrtmp_d=>aptrnxt_d
          aptrnxt_d=>aptrcur_d
          aptrcur_d=>ptrtmp_d
          
          ptrtmp_d=>bptrnxt_d
          bptrnxt_d=>bptrcur_d
          bptrcur_d=>ptrtmp_d

          if((iblock.ge.1).and.(iblock.le.nblocks)) then
             i=m_of_i(iblock) 
             j=n_of_i(iblock) 
             k=k_of_i(iblock) 
             
             cur_buffsize_m=buffersize(ms,bm,i)
             cur_index_m=(i-1)*bm+1
             
             cur_buffsize_n=buffersize(ns,bn,j)
             cur_index_n=(j-1)*bn+1
             
             cur_buffsize_k=buffersize(ks,bk,k)
             cur_index_k=(k-1)*bk+1

             if(lsame(transa,'n').and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,ONE,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else if((.not.lsame(transa,'n')).and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#endif
             else if(lsame(transa,'n').and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else if((.not.lsame(transa,'n')).and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else
                write(*,*) "don't know what to do"
                stop
             end if
          end if

          if(iblock.lt.(nblocks)) then
             nxti=m_of_i(iblock+1) 
             nxtj=n_of_i(iblock+1) 
             nxtk=k_of_i(iblock+1) 

             nxt_buffsize_m=buffersize(ms,bm,nxti)
             nxt_index_m=(nxti-1)*bm+1
             
             nxt_buffsize_n=buffersize(ns,bn,nxtj)
             nxt_index_n=(nxtj-1)*bn+1
             
             nxt_buffsize_k=buffersize(ks,bk,nxtk)
             nxt_index_k=(nxtk-1)*bk+1
             
             call shapelinear_d(transa,aptrnxt_d,a,nxt_index_m,nxt_index_m+nxt_buffsize_m-1,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_buffsize_m,nxt_buffsize_k,lda)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(aptrnxt_d:length(nxt_buffsize_m*nxt_buffsize_k) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalabuff1)
             call shapelinear_d(transb,bptrnxt_d,b,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_index_n,nxt_index_n+nxt_buffsize_n-1,nxt_buffsize_k,nxt_buffsize_n,ldb)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(bptrnxt_d:length(nxt_buffsize_k*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalbbuff2)
             if(mod(iblock,nkb).eq.0) then
                cptrnxt_d(1:nxt_buffsize_m*nxt_buffsize_n)=reshape(c(nxt_index_m:nxt_index_m+nxt_buffsize_m-1,nxt_index_n:nxt_index_n+nxt_buffsize_n-1),(/nxt_buffsize_m*nxt_buffsize_n/))
                cptrnxt_d=beta*cptrnxt_d
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(cptrnxt_d:length(nxt_buffsize_m*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
             end if
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalabuff1) 
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalbbuff2)
          end if
          if(iblock.gt.1) then
             if(mod(iblock-1,nkb).eq.0) then
                prvi=m_of_i(iblock-1) 
                prvj=n_of_i(iblock-1) 
                
                prv_buffsize_m=buffersize(ms,bm,prvi)
                prv_index_m=(prvi-1)*bm+1
                
                prv_buffsize_n=buffersize(ns,bn,prvj)
                prv_index_n=(prvj-1)*bn+1
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) OUT(cptrprv_d:length(prv_buffsize_m*prv_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
                c(prv_index_m:prv_index_m+prv_buffsize_m-1,prv_index_n:prv_index_n+prv_buffsize_n-1)=reshape(cptrprv_d(1:prv_buffsize_m*prv_buffsize_n),(/prv_buffsize_m,prv_buffsize_n/))
             end if
          end if
          if(catchcomputesignal) then
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(computesignal)
                catchcomputesignal=.false.
          end if
       end do

       ! -- deallocate the index data
       deallocate(m_of_i)
       deallocate(n_of_i)
       deallocate(k_of_i)
    end if
    if(offload_verbose) then
       write(*,*) "xphidgemm ",ms,ns,ks,dble(2)*dble(ms)*dble(ns)*dble(ks)/(dclock()-t1)/dble(1000)/dble(1000)/dble(1000),"GFLP"
    endif
  end subroutine xphidgemm



  subroutine xphizgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
    use ifport, only : dclock
    implicit none
#ifdef _BUILD_LIBRARY_
    !dir$ attributes offload:mic                :: MKL_ZGEMM
#else
    !dir$ attributes offload:mic                :: ZGEMM
#endif
    character :: transa,transb
    integer :: lda,ldb,ldc,ms,ks,ns,bm,bn,bk
    complex(kind=kind(0.0D0)), dimension(lda,*) :: a
    complex(kind=kind(0.0D0)), dimension(ldb,*) :: b
    complex(kind=kind(0.0D0)), dimension(ldc,*) :: c
    complex(kind=kind(0.0D0)) :: alpha, beta
    complex(kind=kind(0.0D0)), dimension(:), pointer :: aptrcur
    complex(kind=kind(0.0D0)), dimension(:), pointer :: aptrnxt
    complex(kind=kind(0.0D0)), dimension(:), pointer :: bptrcur
    complex(kind=kind(0.0D0)), dimension(:), pointer :: bptrnxt
    complex(kind=kind(0.0D0)), dimension(:), pointer :: cptrcur
    complex(kind=kind(0.0D0)), dimension(:), pointer :: cptrnxt
    complex(kind=kind(0.0D0)), dimension(:), pointer :: cptrprv
    complex(kind=kind(0.0D0)), dimension(:), pointer :: ptrtmp


    !dir$ attributes offload:mic :: cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, cur_index_m, cur_index_n, cur_index_k, nxt_buffsize_k, nxt_index_k
    integer :: i,j,k,l
    integer :: nxti,nxtj,nxtk
    integer :: prvi,prvj,prvk

    integer :: nkb,nmb,nnb
    integer :: nblocks
    integer :: iblock
    integer :: cur_buffsize_m
    integer :: cur_buffsize_n
    integer :: cur_buffsize_k
    integer :: prv_buffsize_m
    integer :: prv_buffsize_n
    integer :: prv_buffsize_k
    integer :: nxt_buffsize_m
    integer :: nxt_buffsize_n
    integer :: nxt_buffsize_k
    integer :: cur_index_m
    integer :: cur_index_n
    integer :: cur_index_k
    integer :: prv_index_m
    integer :: prv_index_n
    integer :: prv_index_k
    integer :: nxt_index_m
    integer :: nxt_index_n
    integer :: nxt_index_k
    complex(kind=kind(0.0D0)) :: ONE=cmplx(1.0D0,0.0D0)
    complex(kind=kind(0.0D0)) :: ZERO=cmplx(0.0D0,0.0D0)
    integer, dimension(:), allocatable :: m_of_i 
    integer, dimension(:), allocatable :: n_of_i 
    integer, dimension(:), allocatable :: k_of_i 
    integer :: count
    double precision :: t1
    double precision :: tsend
    double precision :: tcomp
    double precision :: trecv
    double precision :: talloc
    double precision :: tdeall
    logical :: catchcomputesignal=.false.

    
    if(.not.isallocated) then
       call xphigemm_allocate(1000,1000,1000)
       isallocated=.true.
    end if

    t1=dclock()
    ! -- check if the matrix computation is actually less than the 
    ! -- threshold performance. In this case, just do it on the host.
    if((dble(8)*dble(ms)*dble(ns)*dble(ks)).lt.(offload_threshold)) then
!       write(*,*) "about to call mkl_zgemm on host"
#ifdef _BUILD_LIBRARY_
       call mkl_zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#else
       call zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#endif
!       call mkl_blas_zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
    else
       ! -- if the total number of FLOPs in the matrix is bigger than the threshold
       ! -- throw the matrix into the blocking algorithm.
       bm=blocksize_m
       bn=blocksize_n
       bk=blocksize_k

       ! -- compute the total number of blocks that must run though
       nmb=numblocks(ms,bm)
       nnb=numblocks(ns,bn)
       nkb=numblocks(ks,bk)
       nblocks=nmb*nnb*nkb


       !       talloc=dclock()
       allocate(m_of_i(nblocks))
       allocate(n_of_i(nblocks))
       allocate(k_of_i(nblocks))

       count=0
       do i=1,nmb
          do j=1,nnb
             do l=1,nkb
                count=count+1
                m_of_i(count)=i
                n_of_i(count)=j
                k_of_i(count)=l
             end do
          end do
       end do


       cptrnxt=>cbuff1
       cptrcur=>cbuff2
       cptrprv=>cbuff3

       aptrnxt=>abuff1
       aptrcur=>abuff2

       bptrnxt=>bbuff1
       bptrcur=>bbuff2


       do iblock=0,nblocks+1
          if( (iblock.eq.1).or.(iblock.eq.(nblocks+1)).or.(m_of_i(iblock).ne.m_of_i(iblock-1)).or.(n_of_i(iblock).ne.n_of_i(iblock-1)) ) then
             ptrtmp=>cptrnxt
             cptrnxt=>cptrprv
             cptrprv=>cptrcur
             cptrcur=>ptrtmp
          end if
          
          ptrtmp=>aptrnxt
          aptrnxt=>aptrcur
          aptrcur=>ptrtmp
          
          ptrtmp=>bptrnxt
          bptrnxt=>bptrcur
          bptrcur=>ptrtmp

          if((iblock.ge.1).and.(iblock.le.nblocks)) then
             i=m_of_i(iblock) 
             j=n_of_i(iblock) 
             k=k_of_i(iblock) 
             
             cur_buffsize_m=buffersize(ms,bm,i)
             cur_index_m=(i-1)*bm+1
             
             cur_buffsize_n=buffersize(ns,bn,j)
             cur_index_n=(j-1)*bn+1
             
             cur_buffsize_k=buffersize(ks,bk,k)
             cur_index_k=(k-1)*bk+1

             if(lsame(transa,'n').and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,ONE,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else if((.not.lsame(transa,'n')).and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else if(lsame(transa,'n').and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else if((.not.lsame(transa,'n')).and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta) signal(computesignal)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#endif
                catchcomputesignal=.true.
             else
                write(*,*) "don't know what to do"
                stop
             end if
          end if

          ! --------------- section 1: SEND ------------------
          if(iblock.lt.(nblocks)) then
             nxti=m_of_i(iblock+1) 
             nxtj=n_of_i(iblock+1) 
             nxtk=k_of_i(iblock+1) 

             nxt_buffsize_m=buffersize(ms,bm,nxti)
             nxt_index_m=(nxti-1)*bm+1
             
             nxt_buffsize_n=buffersize(ns,bn,nxtj)
             nxt_index_n=(nxtj-1)*bn+1
             
             nxt_buffsize_k=buffersize(ks,bk,nxtk)
             nxt_index_k=(nxtk-1)*bk+1
             
             call shapelinear(transa,aptrnxt,a,nxt_index_m,nxt_index_m+nxt_buffsize_m-1,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_buffsize_m,nxt_buffsize_k,lda)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(aptrnxt:length(nxt_buffsize_m*nxt_buffsize_k) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalabuff1)
             call shapelinear(transb,bptrnxt,b,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_index_n,nxt_index_n+nxt_buffsize_n-1,nxt_buffsize_k,nxt_buffsize_n,ldb)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(bptrnxt:length(nxt_buffsize_k*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalbbuff2)
             if(mod(iblock,nkb).eq.0) then
                cptrnxt(1:nxt_buffsize_m*nxt_buffsize_n)=reshape(c(nxt_index_m:nxt_index_m+nxt_buffsize_m-1,nxt_index_n:nxt_index_n+nxt_buffsize_n-1),(/nxt_buffsize_m*nxt_buffsize_n/))
                cptrnxt=beta*cptrnxt
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(cptrnxt:length(nxt_buffsize_m*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
             end if
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalabuff1) 
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalbbuff2)
          end if
          if(iblock.gt.1) then
             if(mod(iblock-1,nkb).eq.0) then
                prvi=m_of_i(iblock-1) 
                prvj=n_of_i(iblock-1) 
                
                prv_buffsize_m=buffersize(ms,bm,prvi)
                prv_index_m=(prvi-1)*bm+1
                
                prv_buffsize_n=buffersize(ns,bn,prvj)
                prv_index_n=(prvj-1)*bn+1
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) OUT(cptrprv:length(prv_buffsize_m*prv_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
                c(prv_index_m:prv_index_m+prv_buffsize_m-1,prv_index_n:prv_index_n+prv_buffsize_n-1)=reshape(cptrprv(1:prv_buffsize_m*prv_buffsize_n),(/prv_buffsize_m,prv_buffsize_n/))
             end if
          end if
          ! --------------- section 2: compute ------------------
          if(catchcomputesignal) then
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(computesignal)
                catchcomputesignal=.false.
          end if
       end do
       deallocate(m_of_i)
       deallocate(n_of_i)
       deallocate(k_of_i)

    end if
    if(offload_verbose) then
       write(*,*) "xphizgemm ",ms,ns,ks,dble(8)*dble(ms)*dble(ns)*dble(ks)/(dclock()-t1)/dble(1000)/dble(1000)/dble(1000),"GFLP"
    end if
  end subroutine xphizgemm

  subroutine xphidgemm_omp(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
    use ifport, only : dclock
    use iso_c_binding
    implicit none
#ifdef _BUILD_LIBRARY_
    !dir$ attributes offload:mic                :: MKL_DGEMM
#else
    !dir$ attributes offload:mic                :: DGEMM
#endif
    character :: transa,transb
    integer :: lda,ldb,ldc,ms,ks,ns,bm,bn,bk
    real(kind=kind(0.0D0)), dimension(lda,*) :: a
    real(kind=kind(0.0D0)), dimension(ldb,*) :: b
    real(kind=kind(0.0D0)), dimension(ldc,*) :: c
    real(kind=kind(0.0D0)) :: alpha, beta
    real(kind=kind(0.0D0)), dimension(:), pointer :: aptrcur_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: aptrnxt_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: bptrcur_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: bptrnxt_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: cptrcur_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: cptrnxt_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: cptrprv_d
    real(kind=kind(0.0D0)), dimension(:), pointer :: ptrtmp_d


    !dir$ attributes offload:mic :: cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, cur_index_m, cur_index_n, cur_index_k, nxt_buffsize_k, nxt_index_k
    integer :: i,j,k,l
    integer :: nxti,nxtj,nxtk
    integer :: prvi,prvj,prvk

    integer :: nkb,nmb,nnb
    integer :: nblocks
    integer :: iblock
    integer :: cur_buffsize_m
    integer :: cur_buffsize_n
    integer :: cur_buffsize_k
    integer :: prv_buffsize_m
    integer :: prv_buffsize_n
    integer :: prv_buffsize_k
    integer :: nxt_buffsize_m
    integer :: nxt_buffsize_n
    integer :: nxt_buffsize_k
    integer :: cur_index_m
    integer :: cur_index_n
    integer :: cur_index_k
    integer :: prv_index_m
    integer :: prv_index_n
    integer :: prv_index_k
    integer :: nxt_index_m
    integer :: nxt_index_n
    integer :: nxt_index_k
    real(kind=kind(0.0D0)) :: ONE=dble(1.0D0)
    real(kind=kind(0.0D0)) :: ZERO=dble(0.0D0)
    integer, dimension(:), allocatable :: m_of_i 
    integer, dimension(:), allocatable :: n_of_i 
    integer, dimension(:), allocatable :: k_of_i 
    integer :: count
    double precision :: t1
    double precision :: tsend
    double precision :: tcomp
    double precision :: trecv
    double precision :: talloc
    double precision :: tdeall

    IF ((ms.EQ.0) .OR. (ns.EQ.0) .OR.(ks.EQ.0)) RETURN
    
    if(.not.isallocated) then
       call xphigemm_allocate(1000,1000,1000)
       isallocated=.true.
    end if

!    write(*,*) "I am here", dble(8)*dble(ms)*dble(ns)*dble(ks), offload_threshold
    t1=dclock()
    ! -- check if the matrix computation is actually less than the 
    ! -- threshold performance. In this case, just do it on the host.
    if((dble(2)*dble(ms)*dble(ns)*dble(ks)).lt.(offload_threshold)) then
!       write(*,*) "about to call mkl_zgemm on host"
#ifdef _BUILD_LIBRARY_
       call mkl_dgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#else
       call dgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#endif
    else
       ! -- if the total number of FLOPs in the matrix is bigger than the threshold
       ! -- throw the matrix into the blocking algorithm.
       bm=blocksize_m
       bn=blocksize_n
       bk=blocksize_k

       ! -- compute the total number of blocks that must run though
       nmb=numblocks(ms,bm)
       nnb=numblocks(ns,bn)
       nkb=numblocks(ks,bk)
       nblocks=nmb*nnb*nkb


       allocate(m_of_i(nblocks))
       allocate(n_of_i(nblocks))
       allocate(k_of_i(nblocks))

       count=0
       ! -- this is recording the plan for the computation 
       ! -- it seems more efficient to record and allocate than
       ! -- to do it in the compute loop. Overhead is in the order 
       ! -- of usecs
       do i=1,nmb
          do j=1,nnb
             do l=1,nkb
                count=count+1
                m_of_i(count)=i
                n_of_i(count)=j
                k_of_i(count)=l
             end do
          end do
       end do


       call c_f_pointer(c_loc(abuff1),aptrnxt_d,[blocksize_m*blocksize_k])
       call c_f_pointer(c_loc(abuff2),aptrcur_d,[blocksize_m*blocksize_k])
       call c_f_pointer(c_loc(bbuff1),bptrnxt_d,[blocksize_k*blocksize_n])
       call c_f_pointer(c_loc(bbuff2),bptrcur_d,[blocksize_k*blocksize_n])
       call c_f_pointer(c_loc(cbuff1),cptrnxt_d,[blocksize_m*blocksize_n])
       call c_f_pointer(c_loc(cbuff2),cptrcur_d,[blocksize_m*blocksize_n])
       call c_f_pointer(c_loc(cbuff3),cptrprv_d,[blocksize_m*blocksize_n])
       

       !$OMP PARALLEL NUM_THREADS(2)
       do iblock=0,nblocks+1
          !$OMP SINGLE
          if( (iblock.eq.1).or.(iblock.eq.(nblocks+1)).or.(m_of_i(iblock).ne.m_of_i(iblock-1)).or.(n_of_i(iblock).ne.n_of_i(iblock-1)) ) then
             ptrtmp_d=>cptrnxt_d
             cptrnxt_d=>cptrprv_d
             cptrprv_d=>cptrcur_d
             cptrcur_d=>ptrtmp_d
          end if
          
          ptrtmp_d=>aptrnxt_d
          aptrnxt_d=>aptrcur_d
          aptrcur_d=>ptrtmp_d
          
          ptrtmp_d=>bptrnxt_d
          bptrnxt_d=>bptrcur_d
          bptrcur_d=>ptrtmp_d
          !$OMP END SINGLE
          !$OMP SECTIONS 
          ! --------------- section 1: SEND ------------------
          !$OMP SECTION
          !   tsend=dclock()
          if(iblock.lt.(nblocks)) then
             nxti=m_of_i(iblock+1) 
             nxtj=n_of_i(iblock+1) 
             nxtk=k_of_i(iblock+1) 

             nxt_buffsize_m=buffersize(ms,bm,nxti)
             nxt_index_m=(nxti-1)*bm+1
             
             nxt_buffsize_n=buffersize(ns,bn,nxtj)
             nxt_index_n=(nxtj-1)*bn+1
             
             nxt_buffsize_k=buffersize(ks,bk,nxtk)
             nxt_index_k=(nxtk-1)*bk+1
             
             call shapelinear_d(transa,aptrnxt_d,a,nxt_index_m,nxt_index_m+nxt_buffsize_m-1,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_buffsize_m,nxt_buffsize_k,lda)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(aptrnxt_d:length(nxt_buffsize_m*nxt_buffsize_k) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalabuff1)
             call shapelinear_d(transb,bptrnxt_d,b,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_index_n,nxt_index_n+nxt_buffsize_n-1,nxt_buffsize_k,nxt_buffsize_n,ldb)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(bptrnxt_d:length(nxt_buffsize_k*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalbbuff1)
             if(mod(iblock,nkb).eq.0) then
                cptrnxt_d(1:nxt_buffsize_m*nxt_buffsize_n)=reshape(c(nxt_index_m:nxt_index_m+nxt_buffsize_m-1,nxt_index_n:nxt_index_n+nxt_buffsize_n-1),(/nxt_buffsize_m*nxt_buffsize_n/))
                cptrnxt_d=beta*cptrnxt_d
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(cptrnxt_d:length(nxt_buffsize_m*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
             end if
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalabuff1) 
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalbbuff1)
          end if
          if(iblock.gt.1) then
             if(mod(iblock-1,nkb).eq.0) then
                prvi=m_of_i(iblock-1) 
                prvj=n_of_i(iblock-1) 
                
                prv_buffsize_m=buffersize(ms,bm,prvi)
                prv_index_m=(prvi-1)*bm+1
                
                prv_buffsize_n=buffersize(ns,bn,prvj)
                prv_index_n=(prvj-1)*bn+1
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) OUT(cptrprv_d:length(prv_buffsize_m*prv_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
                c(prv_index_m:prv_index_m+prv_buffsize_m-1,prv_index_n:prv_index_n+prv_buffsize_n-1)=reshape(cptrprv_d(1:prv_buffsize_m*prv_buffsize_n),(/prv_buffsize_m,prv_buffsize_n/))
             end if
          end if
!          tsend=dclock() -tsend
!          write(*,*) "tsend ",iblock,tsend
          ! --------------- section 2: compute ------------------
          !$OMP SECTION 
!          tcomp=dclock()
          if((iblock.ge.1).and.(iblock.le.nblocks)) then
             i=m_of_i(iblock) 
             j=n_of_i(iblock) 
             k=k_of_i(iblock) 
             
             cur_buffsize_m=buffersize(ms,bm,i)
             cur_index_m=(i-1)*bm+1
             
             cur_buffsize_n=buffersize(ns,bn,j)
             cur_index_n=(j-1)*bn+1
             
             cur_buffsize_k=buffersize(ks,bk,k)
             cur_index_k=(k-1)*bk+1

             if(lsame(transa,'n').and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,ONE,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#endif
             else if((.not.lsame(transa,'n')).and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_k, ONE, cptrcur_d, cur_buffsize_m)
#endif
             else if(lsame(transa,'n').and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_m, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#endif
             else if((.not.lsame(transa,'n')).and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur_d:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta)
#ifdef _BUILD_LIBRARY_
                call mkl_dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#else
                call dgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur_d, cur_buffsize_k, bptrcur_d, cur_buffsize_n, ONE, cptrcur_d, cur_buffsize_m)
#endif
             else
                write(*,*) "don't know what to do"
                stop
             end if
          end if
!          tcomp=dclock()-tcomp
!          write(*,*) "tcomp ",iblock, tcomp
          !$OMP END SECTIONS 
       end do
       !$OMP END PARALLEL

       ! -- deallocate the index data
       deallocate(m_of_i)
       deallocate(n_of_i)
       deallocate(k_of_i)
    end if
    if(offload_verbose) then
       write(*,*) "xphidgemm ",ms,ns,ks,dble(2)*dble(ms)*dble(ns)*dble(ks)/(dclock()-t1)/dble(1000)/dble(1000)/dble(1000),"GFLP"
    endif
  end subroutine xphidgemm_omp



  subroutine xphizgemm_omp(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
    use ifport, only : dclock
    implicit none
#ifdef _BUILD_LIBRARY_
    !dir$ attributes offload:mic                :: MKL_ZGEMM
#else
    !dir$ attributes offload:mic                :: ZGEMM
#endif
    character :: transa,transb
    integer :: lda,ldb,ldc,ms,ks,ns,bm,bn,bk
    complex(kind=kind(0.0D0)), dimension(lda,*) :: a
    complex(kind=kind(0.0D0)), dimension(ldb,*) :: b
    complex(kind=kind(0.0D0)), dimension(ldc,*) :: c
    complex(kind=kind(0.0D0)) :: alpha, beta
    complex(kind=kind(0.0D0)), dimension(:), pointer :: aptrcur
    complex(kind=kind(0.0D0)), dimension(:), pointer :: aptrnxt
    complex(kind=kind(0.0D0)), dimension(:), pointer :: bptrcur
    complex(kind=kind(0.0D0)), dimension(:), pointer :: bptrnxt
    complex(kind=kind(0.0D0)), dimension(:), pointer :: cptrcur
    complex(kind=kind(0.0D0)), dimension(:), pointer :: cptrnxt
    complex(kind=kind(0.0D0)), dimension(:), pointer :: cptrprv
    complex(kind=kind(0.0D0)), dimension(:), pointer :: ptrtmp


    !dir$ attributes offload:mic :: cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, cur_index_m, cur_index_n, cur_index_k, nxt_buffsize_k, nxt_index_k
    integer :: i,j,k,l
    integer :: nxti,nxtj,nxtk
    integer :: prvi,prvj,prvk

    integer :: nkb,nmb,nnb
    integer :: nblocks
    integer :: iblock
    integer :: cur_buffsize_m
    integer :: cur_buffsize_n
    integer :: cur_buffsize_k
    integer :: prv_buffsize_m
    integer :: prv_buffsize_n
    integer :: prv_buffsize_k
    integer :: nxt_buffsize_m
    integer :: nxt_buffsize_n
    integer :: nxt_buffsize_k
    integer :: cur_index_m
    integer :: cur_index_n
    integer :: cur_index_k
    integer :: prv_index_m
    integer :: prv_index_n
    integer :: prv_index_k
    integer :: nxt_index_m
    integer :: nxt_index_n
    integer :: nxt_index_k
    complex(kind=kind(0.0D0)) :: ONE=cmplx(1.0D0,0.0D0)
    complex(kind=kind(0.0D0)) :: ZERO=cmplx(0.0D0,0.0D0)
    integer, dimension(:), allocatable :: m_of_i 
    integer, dimension(:), allocatable :: n_of_i 
    integer, dimension(:), allocatable :: k_of_i 
    integer :: count
    double precision :: t1
    double precision :: tsend
    double precision :: tcomp
    double precision :: trecv
    double precision :: talloc
    double precision :: tdeall


    
    if(.not.isallocated) then
       call xphigemm_allocate(1000,1000,1000)
       isallocated=.true.
    end if

    t1=dclock()
    ! -- check if the matrix computation is actually less than the 
    ! -- threshold performance. In this case, just do it on the host.
    if((dble(8)*dble(ms)*dble(ns)*dble(ks)).lt.(offload_threshold)) then
!       write(*,*) "about to call mkl_zgemm on host"
#ifdef _BUILD_LIBRARY_
       call mkl_zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#else
       call zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
#endif
!       call mkl_blas_zgemm(transa, transb, ms, ns, ks, alpha, a, lda, b, ldb, beta, c, ldc)
    else
       ! -- if the total number of FLOPs in the matrix is bigger than the threshold
       ! -- throw the matrix into the blocking algorithm.
       bm=blocksize_m
       bn=blocksize_n
       bk=blocksize_k

       ! -- compute the total number of blocks that must run though
       nmb=numblocks(ms,bm)
       nnb=numblocks(ns,bn)
       nkb=numblocks(ks,bk)
       nblocks=nmb*nnb*nkb


       !       talloc=dclock()
       allocate(m_of_i(nblocks))
       allocate(n_of_i(nblocks))
       allocate(k_of_i(nblocks))

       count=0
       do i=1,nmb
          do j=1,nnb
             do l=1,nkb
                count=count+1
                m_of_i(count)=i
                n_of_i(count)=j
                k_of_i(count)=l
             end do
          end do
       end do


       cptrnxt=>cbuff1
       cptrcur=>cbuff2
       cptrprv=>cbuff3

       aptrnxt=>abuff1
       aptrcur=>abuff2

       bptrnxt=>bbuff1
       bptrcur=>bbuff2


       !$OMP PARALLEL NUM_THREADS(2)
       do iblock=0,nblocks+1
          !$OMP SINGLE
          if( (iblock.eq.1).or.(iblock.eq.(nblocks+1)).or.(m_of_i(iblock).ne.m_of_i(iblock-1)).or.(n_of_i(iblock).ne.n_of_i(iblock-1)) ) then
             ptrtmp=>cptrnxt
             cptrnxt=>cptrprv
             cptrprv=>cptrcur
             cptrcur=>ptrtmp
          end if
          
          ptrtmp=>aptrnxt
          aptrnxt=>aptrcur
          aptrcur=>ptrtmp
          
          ptrtmp=>bptrnxt
          bptrnxt=>bptrcur
          bptrcur=>ptrtmp
          !$OMP END SINGLE
          !$OMP SECTIONS 
          ! --------------- section 1: SEND ------------------
          !$OMP SECTION
!          tsend=dclock()
          if(iblock.lt.(nblocks)) then
             nxti=m_of_i(iblock+1) 
             nxtj=n_of_i(iblock+1) 
             nxtk=k_of_i(iblock+1) 

             nxt_buffsize_m=buffersize(ms,bm,nxti)
             nxt_index_m=(nxti-1)*bm+1
             
             nxt_buffsize_n=buffersize(ns,bn,nxtj)
             nxt_index_n=(nxtj-1)*bn+1
             
             nxt_buffsize_k=buffersize(ks,bk,nxtk)
             nxt_index_k=(nxtk-1)*bk+1
             
             call shapelinear(transa,aptrnxt,a,nxt_index_m,nxt_index_m+nxt_buffsize_m-1,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_buffsize_m,nxt_buffsize_k,lda)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(aptrnxt:length(nxt_buffsize_m*nxt_buffsize_k) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalabuff1)
             call shapelinear(transb,bptrnxt,b,nxt_index_k,nxt_index_k+nxt_buffsize_k-1,nxt_index_n,nxt_index_n+nxt_buffsize_n-1,nxt_buffsize_k,nxt_buffsize_n,ldb)
             !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(bptrnxt:length(nxt_buffsize_k*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) signal(signalbbuff1)
             if(mod(iblock,nkb).eq.0) then
                cptrnxt(1:nxt_buffsize_m*nxt_buffsize_n)=reshape(c(nxt_index_m:nxt_index_m+nxt_buffsize_m-1,nxt_index_n:nxt_index_n+nxt_buffsize_n-1),(/nxt_buffsize_m*nxt_buffsize_n/))
                cptrnxt=beta*cptrnxt
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) IN(cptrnxt:length(nxt_buffsize_m*nxt_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
             end if
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalabuff1) 
             !DIR$ OFFLOAD_WAIT TARGET(MIC:offload_device) wait(signalbbuff1)
          end if
          if(iblock.gt.1) then
             if(mod(iblock-1,nkb).eq.0) then
                prvi=m_of_i(iblock-1) 
                prvj=n_of_i(iblock-1) 
                
                prv_buffsize_m=buffersize(ms,bm,prvi)
                prv_index_m=(prvi-1)*bm+1
                
                prv_buffsize_n=buffersize(ns,bn,prvj)
                prv_index_n=(prvj-1)*bn+1
                !DIR$ OFFLOAD_TRANSFER TARGET(MIC:offload_device) OUT(cptrprv:length(prv_buffsize_m*prv_buffsize_n) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) 
                c(prv_index_m:prv_index_m+prv_buffsize_m-1,prv_index_n:prv_index_n+prv_buffsize_n-1)=reshape(cptrprv(1:prv_buffsize_m*prv_buffsize_n),(/prv_buffsize_m,prv_buffsize_n/))
             end if
          end if
          ! --------------- section 2: compute ------------------
          !$OMP SECTION 
          if((iblock.ge.1).and.(iblock.le.nblocks)) then
             i=m_of_i(iblock) 
             j=n_of_i(iblock) 
             k=k_of_i(iblock) 
             
             cur_buffsize_m=buffersize(ms,bm,i)
             cur_index_m=(i-1)*bm+1
             
             cur_buffsize_n=buffersize(ns,bn,j)
             cur_index_n=(j-1)*bn+1
             
             cur_buffsize_k=buffersize(ks,bk,k)
             cur_index_k=(k-1)*bk+1

             if(lsame(transa,'n').and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALIGN(64) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,ONE,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#endif
             else if((.not.lsame(transa,'n')).and.lsame(transb,'n')) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_k, ONE, cptrcur, cur_buffsize_m)
#endif
             else if(lsame(transa,'n').and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_m, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#endif
             else if((.not.lsame(transa,'n')).and.(.not.lsame(transb,'n'))) then 
                !DIR$ OFFLOAD TARGET(MIC:offload_device) & 
                !DIR$ & IN(aptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(bptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) &
                !DIR$ & IN(cptrcur:length(0) ALLOC_IF(.false.) FREE_IF(.false.)) IN(transa,transb,cur_buffsize_m,cur_buffsize_n,cur_buffsize_k,alpha,beta)
#ifdef _BUILD_LIBRARY_
                call mkl_zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#else
                call zgemm(transa, transb, cur_buffsize_m, cur_buffsize_n, cur_buffsize_k, alpha, aptrcur, cur_buffsize_k, bptrcur, cur_buffsize_n, ONE, cptrcur, cur_buffsize_m)
#endif
             else
                write(*,*) "don't know what to do"
                stop
             end if
          end if
          !$OMP END SECTIONS 
       end do
       !$OMP END PARALLEL
       deallocate(m_of_i)
       deallocate(n_of_i)
       deallocate(k_of_i)

    end if
    if(offload_verbose) then
       write(*,*) "xphizgemm ",ms,ns,ks,dble(8)*dble(ms)*dble(ns)*dble(ks)/(dclock()-t1)/dble(1000)/dble(1000)/dble(1000),"GFLP"
    end if
  end subroutine xphizgemm_omp
end module xphilibmod

