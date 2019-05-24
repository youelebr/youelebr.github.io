
!DIR$ MAQAO SPECIALIZE(choice=1, paw_opt=3, cplex=2)
!DIR$ MAQAO SPECIALIZE(choice=1, paw_opt<3, cplex=2)
!DIR$ MAQAO SPECIALIZE(choice=1, paw_opt>3, cplex=2)
subroutine opernlb_ylm(choice,cplex,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
& ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
& nspinor,paw_opt,ph3d,svect,ucvol,vect)

! use defs_basis
! use m_profiling_abi
! use m_errors




!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.


!End of the abilint section

 implicit none

  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp)) ! Complex should not be used presently
                                                   ! except for use of libraries
  complex(dpc), parameter :: czero=(0._dp,0._dp)

  real(dp), parameter :: zero=0._dp
  real(dp), parameter :: four=4._dp
  real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: two_pi=2._dp*pi
  real(dp), parameter :: four_pi=four*pi


!Arguments ------------------------------------
!scalars

 integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,idir,matblk,ndgxdtfac,nincat
 integer,intent(in) :: nkpg,nlmn,npw,nspinor,paw_opt
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: cplex_dgxdt(ndgxdtfac),indlmn(6,nlmn),nloalg(5)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3)),vect(2,npw*nspinor)
!Local variables-------------------------------
!Arrays
!scalars
 integer :: fdb,fdf,ia,iaph3d,ic,ii,il,ilmn,ipw,ipwshft,ispinor,jc,nthreads,k
 real(dp) :: scale,wt
 logical :: parity
!arrays
 integer,parameter :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 real(dp),allocatable :: dgxdtfac_(:,:,:),dgxdtfacs_(:,:,:),gxfac_(:,:),gxfacs_(:,:)
 complex(dpc),allocatable :: ztab(:)

! *************************************************************************



!Nothing to do when choice=4, 6 or 23
 if (choice==4.or.choice==6.or.choice==23) return

!DDK not compatible with istwkf > 1
 if(cplex==1.and.any(cplex_dgxdt(:)==2))then
  
 end if

!Inits
 wt=four_pi/sqrt(ucvol)
 nthreads=1




 if (paw_opt/=3) then
   allocate(gxfac_ (2,nlmn))
   gxfac_(:,:)=zero
   if (choice>1) then
     allocate(dgxdtfac_ (2,ndgxdtfac,nlmn))
     if(ndgxdtfac>0) dgxdtfac_(:,:,:)=zero
   end if
 end if
 if (paw_opt>=3) then
   allocate(gxfacs_ (2,nlmn))
   gxfacs_(:,:)=zero
   if (choice>1) then
     allocate(dgxdtfacs_ (2,ndgxdtfac,nlmn))
     if (ndgxdtfac>0) dgxdtfacs_(:,:,:)=zero
   end if
 end if

!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 if (nthreads==1) then

! Loop on spinorial components
   do ispinor=1,nspinor
     ipwshft=(ispinor-1)*npw

! Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

! Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxfac_(2,ilmn)=zero
           else
             gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
             else
               gxfac_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn)= scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfac_(jc,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfac_(jc,ii,ilmn)= scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
       end if

! Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       if (paw_opt>=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
             if (cplex==1) gxfacs_(2,ilmn)=zero
           else
             gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
             if (cplex==2) then
               gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
             else
               gxfacs_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn)= scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfacs_(jc,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfacs_(jc,ii,ilmn)= scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
       end if

       allocate(ztab (npw))

! Compute <g|Vnl|c> (or derivatives) for each plane wave:

       if (paw_opt/=3) then

         ztab(:)=czero

! ------
         if (choice==1) then ! <g|Vnl|c>
!D I R $  MAQAO TILE_INNER_IF_SPE_choicee1=8
           do ilmn=1,nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
do k=1,npw
ztab(k) = ztab(k)+ffnl(k,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
end do
           end do
         end if

! ------
         if (choice==2) then ! derivative w.r.t. atm. pos
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir)*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
& *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp)&
& -ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           else
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
& -ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end if
         end if

! ------
         if (choice==5) then ! full derivative w.r.t. k
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
& +ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
! -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) + &
& ffnl(:,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) - &
& ffnl(:,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if


! ------
         ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         vect(1,1+ipwshft:npw+ipwshft)=vect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
         vect(2,1+ipwshft:npw+ipwshft)=vect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
       end if

! Compute <g|S|c> (or derivatives) for each plane wave:

       if (paw_opt>=3) then

         ztab(:)=czero

! ------
         if (choice==1) then ! <g|S|c>
!DIR$ MAQAO TILE_INNER_IF_SPE_choicee1=8
           do ilmn=1,nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
do k=1,npw
ztab(k) = ztab(k)+ffnl(k,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
end do
           end do
         end if

! ------
         if (choice==2) then ! derivative w.r.t. atm. pos
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir)*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
& *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
& -ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           else
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
& -ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end if
         end if

! ------
         if (choice==5) then ! full derivative w.r.t. k
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
& +ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

! ------
         if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
! -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) + &
& ffnl(:,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) - &
& ffnl(:,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if


! ------
         ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
         svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
       end if

       deallocate(ztab)

! End loop on atoms
     end do
   end do ! End loop on spinors


! ==========================================================================
! ========== OPENMP VERSION ================================================
! ==========================================================================
 else

! Loop on spinorial components
   do ispinor=1,nspinor
     ipwshft=(ispinor-1)*npw

! Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

! Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxfac_(2,ilmn)=zero
           else
             gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
             else
               gxfac_(1,ilmn)=zero
             end if
           end if
         end do
!$OMP END DO
         if (choice>1) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn)= scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfac_(jc,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfac_(jc,ii,ilmn)= scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
!$OMP END DO
         end if
!$OMP END PARALLEL
       end if

! Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       if (paw_opt>=3) then
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
             if (cplex==1) gxfacs_(2,ilmn)=zero
           else
             gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
             if (cplex==2) then
               gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
             else
               gxfacs_(1,ilmn)=zero
             end if
           end if
         end do
!$OMP END DO
         if (choice>1) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn)= scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfacs_(jc,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfacs_(jc,ii,ilmn)= scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
!$OMP END DO
         end if
!$OMP END PARALLEL
       end if

       allocate(ztab (npw))

! Compute <g|Vnl|c> (or derivatives) for each plane wave:
       if (paw_opt/=3) then
!$OMP PARALLEL PRIVATE(ipw,ilmn,fdf,fdb)

! ------
         if (choice==1) then ! <g|Vnl|c>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==2) then ! derivative w.r.t. atm. pos
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir)*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn) &
& *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp) &
& -ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           else
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) &
& -ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           end if


! ------
         else if (choice==5) then ! full derivative w.r.t. k
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) &
& +ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
! -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
& +ffnl(ipw,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) &
& -ffnl(ipw,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO
         else
!$OMP WORKSHARE
           ztab(:)=czero
!$OMP END WORKSHARE
         end if

! ------
!$OMP DO
         do ipw=1,npw
           ztab(ipw)=ztab(ipw)*cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           vect(1,ipw+ipwshft)=vect(1,ipw+ipwshft)+real(ztab(ipw))
           vect(2,ipw+ipwshft)=vect(2,ipw+ipwshft)+aimag(ztab(ipw))
         end do
!$OMP END DO

!$OMP END PARALLEL
       end if

! Compute <g|S|c> (or derivatives) for each plane wave:
       if (paw_opt>=3) then
!$OMP PARALLEL PRIVATE(ilmn,ipw,fdf,fdb)

! ------
         if (choice==1) then ! <g|S|c>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==2) then ! derivative w.r.t. atm. pos
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir)*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn) &
& *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
& -ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           else
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) &
& -ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           end if

! ------
         else if (choice==5) then ! full derivative w.r.t. k
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) &
& +ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

! ------
         else if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi> -
! <G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
& +ffnl(ipw,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) &
& -ffnl(ipw,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO
         else
!$OMP WORKSHARE
           ztab(:)=czero
!$OMP END WORKSHARE
         end if


! ------
! The OMP WORKSHARE directive doesn't have a good performance with Intel Compiler
! !$OMP WORKSHARE
! ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
! !$OMP END WORKSHARE
! !$OMP WORKSHARE
! svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
! svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
! !$OMP END WORKSHARE
!$OMP DO
         do ipw=1,npw
           ztab(ipw)=ztab(ipw)*cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           svect(1,ipw+ipwshft)=svect(1,ipw+ipwshft)+real(ztab(ipw))
           svect(2,ipw+ipwshft)=svect(2,ipw+ipwshft)+aimag(ztab(ipw))
         end do
!$OMP END DO
!$OMP END PARALLEL
       end if

       deallocate(ztab)

! End loop on atoms
     end do
! End loop on spinors
   end do

! ==========================================================================
 end if

 if (paw_opt/=3) then
   deallocate(gxfac_)
   if (choice>1) then
     deallocate(dgxdtfac_)
   end if
 end if
 if (paw_opt>=3) then
   deallocate(gxfacs_)
   if (choice>1) then
     deallocate(dgxdtfacs_)
   end if
 end if




!Fake use of unused variable
! if (.false.) write(std_out,*) ipw


end subroutine opernlb_ylm
