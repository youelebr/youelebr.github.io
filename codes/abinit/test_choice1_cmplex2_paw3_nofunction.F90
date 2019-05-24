!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   choice 1    paw_opt 3 cplex 2   v1       !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  define ABI_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE) 

#  define ABI_DEALLOCATE(ARR) \
   deallocate(ARR) 

Program testpaw3

  implicit none

  INTERFACE 
     FUNCTION wallclock()
       integer, parameter :: dp=kind(1.0d0)
       real(dp) :: wallclock
     END FUNCTION wallclock
subroutine opernlb_ylm(choice,cplex,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
& ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
& nspinor,paw_opt,ph3d,svect,ucvol,vect)
 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
 integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,idir,matblk,ndgxdtfac,nincat
 integer,intent(in) :: nkpg,nlmn,npw,nspinor,paw_opt
 real(dp),intent(in) :: ucvol
 integer,intent(in) :: cplex_dgxdt(ndgxdtfac),indlmn(6,nlmn),nloalg(5)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3)),vect(2,npw*nspinor)
end subroutine opernlb_ylm

  END INTERFACE

  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp)) 

  !scalars
  integer :: cplex
  integer :: cplex_fac
  integer :: dimffnl
  integer :: ia3
  integer :: nincat
  integer :: matblk
  integer :: nlmn,npw,nspinor
  integer :: paw_opt
  real(dp) :: ucvol

  !arrays
!!$  integer,intent(in) ::  indlmn(6,nlmn)
!!$  integer,intent(in) ::  nloalg(5)
!!$  real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn)
!!$  real(dp),intent(in) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
!!$  real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
!!$  real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3))
!!$  real(dp),intent(in) :: ph3d(2,npw,matblk)

  integer,  allocatable :: indlmn(:,:)
  integer,  allocatable :: nloalg(:)
  real(dp), allocatable :: ffnl(:,:,:)
  real(dp), allocatable :: gxfac(:,:,:,:)
  real(dp), allocatable :: gxfac_sij(:,:,:,:)
  real(dp), allocatable :: svect(:,:)
  real(dp), allocatable :: svect_perf(:,:)
  real(dp), allocatable :: svect_ref(:,:)
  real(dp), allocatable :: ph3d(:,:,:)

! integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,idir,matblk,ndgxdtfac,nincat
! integer,intent(in) :: nkpg,nlmn,npw,nspinor,paw_opt
! real(dp),intent(in) :: ucvol
! integer,intent(in) :: cplex_dgxdt(ndgxdtfac),indlmn(6,nlmn),nloalg(5)
! real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
! real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
! real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),gxfac(cplex_fac,nlmn,nincat,nspinor)
! real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
! real(dp),intent(in) :: kpg(npw,nkpg),ph3d(2,npw,matblk)
! real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3)),vect(2,npw*nspinor)

  real(dp), pointer, dimension(:,:,:,:,:) :: dgxdtfac => NULL()
  real(dp), pointer, dimension(:,:,:,:,:) :: dgxdtfac_sij => NULL()
  integer, pointer, dimension(:) :: CPLEX_DGXDT => NULL()
  real(dp), pointer, dimension(:,:) :: VECT => NULL()
  integer :: idir = 0
  real(dp), pointer, dimension(:,:) :: KPG => NULL()
  integer :: ndgxdtfac = 0
  integer :: nkpg = 0

  integer  :: i
  real(dp) :: t1,t2

  !1657
  !324971 
  !550000

  !open file 
  open( UNIT=11 , FILE ="mydata_paw3_550000.data", FORM ="unformatted", ACCESS ="sequential", ACTION ="read")

  !get scalars
  read(11) cplex
  read(11) cplex_fac
  read(11) dimffnl
  read(11) ia3
  read(11) nincat
  read(11) matblk
  read(11) nlmn,npw,nspinor
  read(11) paw_opt
  read(11) ucvol

  !allocate arrays
  allocate(indlmn(6,nlmn))
  allocate(nloalg(5))
  allocate(ffnl(npw,dimffnl,nlmn))
  allocate(gxfac(cplex_fac,nlmn,nincat,nspinor))
  allocate(gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3)))
  allocate(svect(2,npw*nspinor*(paw_opt/3)))
  allocate(svect_perf(2,npw*nspinor*(paw_opt/3)))
  allocate(svect_ref(2,npw*nspinor*(paw_opt/3)))
  allocate(ph3d(2,npw,matblk))

  !get arrays values
  read(11) indlmn
  read(11) nloalg
  read(11) ffnl
  read(11) gxfac
  read(11) gxfac_sij
  read(11) svect
  read(11) ph3d
  read(11) svect_ref

  svect_perf = svect

!subroutine opernlb_ylm(choice,cplex,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
!& ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
!& nspinor,paw_opt,ph3d,svect,ucvol,vect)

  t1 = wallclock()
  do i=0,5000
	call opernlb_ylm(1,2,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,ia3,idir,indlmn,&
                  kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,3,ph3d,svect_perf,ucvol,vect)
  end do
  t2 = wallclock()
  print *, "time (sec): ", t2-t1


  !check results
  call opernlb_ylm(1,2,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,ia3,idir,indlmn,&
                  kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,3,ph3d,svect_ref,ucvol,vect)
  print *, svect(1,15),   svect_ref(1,15)
  print *, svect(2,15),   svect_ref(2,15)
  print *, svect(1,150),  svect_ref(1,150)
  print *, svect(2,150),  svect_ref(2,150)
  print *, svect(1,1500), svect_ref(1,1500)
  print *, svect(2,1500), svect_ref(2,1500)


  !deallocate arrays
  deallocate(indlmn)
  deallocate(nloalg)
  deallocate(ffnl)
  deallocate(gxfac)
  deallocate(gxfac_sij)
  deallocate(svect)
  deallocate(svect_perf)
  deallocate(svect_ref)
  deallocate(ph3d)

  !close file
  close(11)


End Program testpaw3
