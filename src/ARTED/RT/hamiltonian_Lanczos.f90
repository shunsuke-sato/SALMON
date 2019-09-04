!
!  Copyright 2017-2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
#define TIMER_BEG(id) call timer_thread_begin(id)
#define TIMER_END(id) call timer_thread_end(id)

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

subroutine hamiltonian_Lanczos(zu,flag_current)
  use Global_Variables
  use timer
  use omp_lib
  use opt_variables
  use hpsi, only: hpsi_omp_KB_RT
  implicit none
  integer    :: tid
  integer    :: ikb,ik,ib,i
  integer, parameter :: ndim_Lanczos = 4
  complex(8) :: zpsi_t(0:PNL-1,ndim_Lanczos,0:NUMBER_THREADS-1)
  real(8) :: alpha_L(ndim_Lanczos,0:NUMBER_THREADS-1)
  real(8) :: beta_L(ndim_Lanczos,0:NUMBER_THREADS-1)
  real(8) :: ham_L(ndim_Lanczos,ndim_Lanczos,0:NUMBER_THREADS-1)
  real(8) :: ss
  complex(8) :: zvec(ndim_Lanczos,0:NUMBER_THREADS-1)
  integer :: j
  complex(8), intent(inout) :: zu(NL,NBoccmax,NK_s:NK_e)
  logical, intent(in)       :: flag_current
!LAPACK ==
      integer :: lwork
      real(8),allocatable :: work_lp(:,:)
      real(8),allocatable :: rwork(:,:),w(:,:)
      integer :: info
!LAPACK ==

#ifndef ARTED_CURRENT_PREPROCESSING
#define UNUSED_VARIABLE(VAR) if(.false.) call salmon_unusedvar(VAR)
  UNUSED_VARIABLE(flag_current)
#endif

  lwork=6*ndim_Lanczos**2
  allocate(work_lp(lwork,0:NUMBER_THREADS-1))
  allocate(rwork(3*ndim_Lanczos-2,0:NUMBER_THREADS-1))
  allocate(w(ndim_Lanczos,0:NUMBER_THREADS-1))

  call timer_begin(LOG_HPSI)

  tid = 0
!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

!$omp do private(ik,ib,ss,info)
  do ikb=1,NKB
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    call init(zpsi_t(:,1,tid),zu(:,ib,ik))
    ss = sum(abs(zpsi_t(:,1,tid))**2)*Hxyz
    zpsi_t(:,1,tid) = zpsi_t(:,1,tid)/sqrt(ss)
    call hpsi_omp_KB_RT(ik,zpsi_t(:,1,tid),zhtpsi(:,1,tid))
    alpha_L(1,tid) = sum(conjg(zhtpsi(:,1,tid))*zpsi_t(:,1,tid))*Hxyz
    zhtpsi(:,1,tid) = zhtpsi(:,1,tid) - alpha_L(1,tid)*zpsi_t(:,1,tid)

    do j = 2, ndim_Lanczos

      beta_L(j,tid) = sqrt( sum(abs(zhtpsi(:,1,tid))**2)*Hxyz )
      zpsi_t(:,j,tid) = zhtpsi(:,1,tid)/beta_L(j,tid)
      call hpsi_omp_KB_RT(ik,zpsi_t(:,j,tid),zhtpsi(:,1,tid))
      alpha_L(j,tid) = sum(conjg(zhtpsi(:,1,tid))*zpsi_t(:,j,tid))*Hxyz
      zhtpsi(:,1,tid) = zhtpsi(:,1,tid) - alpha_L(j,tid)*zpsi_t(:,j,tid) &
                                        -  beta_L(j,tid)*zpsi_t(:,j-1,tid)

    end do

    ham_L(:,:,tid) = 0d0
    ham_L(1,1,tid) = alpha_L(1,tid)
    do j = 2, ndim_Lanczos
      ham_L(j,j,tid) = alpha_L(j,tid)
      ham_L(j,j-1,tid) = beta_L(j,tid)
      ham_L(j-1,j,tid) = beta_L(j,tid)
    end do

    call dsyev('V', 'U', ndim_Lanczos, ham_L(:,:,tid), ndim_Lanczos &
      , w(:,tid), work_lp(:,tid), lwork, info) 

    zvec(:,tid) = 0d0; zvec(1,tid) = 1d0
    zvec(:,tid) = matmul(transpose(ham_L(:,:,tid)),zvec(:,tid)) 
    zvec(:,tid) = exp(-zI*dt*w(:,tid))*zvec(:,tid)
    zvec(:,tid) = matmul(ham_L(:,:,tid),zvec(:,tid)) 
    zhtpsi(:,1,tid) = matmul(zpsi_t(:,:,tid),zvec(:,tid))

    call update(zhtpsi(:,1,tid),zu(:,ib,ik))

#ifdef ARTED_CURRENT_PREPROCESSING
    if(flag_current) call current_omp_KB_ST(ib,ik,zu(:,ib,ik))
#endif
  end do
!$omp end do
!$omp end parallel

  call timer_end(LOG_HPSI)

contains
  subroutine init(tpsi,zu)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timer
    implicit none
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
    integer :: ix,iy,iz

    TIMER_BEG(LOG_HPSI_INIT)
!dir$ vector aligned
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      tpsi(iz,iy,ix)=zu(iz,iy,ix)
    end do
    end do
    end do
    TIMER_END(LOG_HPSI_INIT)
  end subroutine

  subroutine update(tpsi,zu)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timer
    implicit none
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
    integer :: ix,iy,iz

    TIMER_BEG(LOG_HPSI_UPDATE)
!dir$ vector aligned
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      zu(iz,iy,ix)=tpsi(iz,iy,ix)
    end do
    end do
    end do
    TIMER_END(LOG_HPSI_UPDATE)
  end subroutine
end subroutine
