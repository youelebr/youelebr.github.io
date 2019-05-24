SUBROUTINE opernlb_ylm(choice,cplex,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,paw_opt,ph3d,svect,ucvol,vect)
 IMPLICIT NONE

  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER, PARAMETER :: dpc = kind(((1.0_dp, 1.0_dp)))
  COMPLEX(kind=dpc), PARAMETER :: czero = (0._dp, 0._dp)
  REAL(kind=dp), PARAMETER :: zero = 0._dp
  REAL(kind=dp), PARAMETER :: four = 4._dp
  REAL(kind=dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER :: two_pi = 2._dp * pi
  REAL(kind=dp), PARAMETER :: four_pi = four * pi
  INTEGER, INTENT(IN) :: choice, cplex, cplex_fac, dimffnl, ia3, idir, matblk, ndgxdtfac, nincat
  INTEGER, INTENT(IN) :: nkpg, nlmn, npw, nspinor, paw_opt
  REAL(kind=dp), INTENT(IN) :: ucvol
  INTEGER, INTENT(IN) :: cplex_dgxdt(ndgxdtfac), indlmn(6,nlmn), nloalg(5)
  REAL(kind=dp), DIMENSION(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor), INTENT(IN) :: dgxdtfac
  REAL(kind=dp), DIMENSION(cplex,ndgxdtfac,nlmn,nincat * (paw_opt / 3),nspinor), INTENT(IN) :: dgxdtfac_sij
  REAL(kind=dp), INTENT(IN) :: ffnl(npw,dimffnl,nlmn), gxfac(cplex_fac,nlmn,nincat,nspinor)
  REAL(kind=dp), DIMENSION(cplex,nlmn,nincat,nspinor * (paw_opt / 3)), INTENT(IN) :: gxfac_sij
  REAL(kind=dp), INTENT(IN) :: kpg(npw,nkpg), ph3d(2,npw,matblk)
  REAL(kind=dp), INTENT(INOUT) :: svect(2,npw * nspinor * (paw_opt / 3)), vect(2,npw * nspinor)
!Local variables-------------------------------
!Arrays
!scalars
  INTEGER :: fdb, fdf, ia, iaph3d, ic, ii, il, ilmn, ipw, ipwshft, ispinor, jc, nthreads, k
  REAL(kind=dp) :: scale, wt
  LOGICAL :: parity
!arrays
  INTEGER, DIMENSION(6), PARAMETER :: ffnl_dir_dat = (/3,4,4,2,2,3/)
  REAL(kind=dp), ALLOCATABLE :: dgxdtfac_(:,:,:), dgxdtfacs_(:,:,:), gxfac_(:,:), gxfacs_(:,:)
  COMPLEX(kind=dpc), DIMENSION(:), ALLOCATABLE :: ztab
 IF (choice .EQ. 1 .AND. paw_opt > 3 .AND. cplex .EQ. 2)  THEN
CALL opernlb_ylm_che1_pae3_cpe2(cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn, kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,paw_opt,ph3d,svect, &
  & ucvol,vect)
 RETURN
 END IF
 IF (choice .EQ. 1 .AND. paw_opt < 3 .AND. cplex .EQ. 2)  THEN
CALL opernlb_ylm_ASSIST_choicee1_paw_opti3_cplexe2(cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,paw_opt,ph3d,svect, &
  & ucvol,vect)
 RETURN
 END IF
 IF (choice .EQ. 1 .AND. paw_opt .EQ. 3 .AND. cplex .EQ. 2)  THEN
CALL opernlb_ylm_ASSIST_choicee1_paw_opte3_cplexe2(cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,ph3d,svect,ucvol,vect)
 RETURN
 END IF
! *************************************************************************
!Nothing to do when choice=4, 6 or 23
 IF (choice .EQ. 4 .OR. (choice .EQ. 6 .OR. choice .EQ. 23)) RETURN
!DDK not compatible with istwkf > 1
 IF (cplex .EQ. 1 .AND. any(cplex_dgxdt(:) .EQ. 2))  THEN
 END IF
!Inits
  wt = four_pi / sqrt(ucvol)
  nthreads = 1
 IF (paw_opt .NE. 3)  THEN
   allocate( gxfac_(2,nlmn) )
   gxfac_(:,:) = zero
   IF (choice > 1)  THEN
     allocate( dgxdtfac_(2,ndgxdtfac,nlmn) )
     IF (ndgxdtfac > 0)        dgxdtfac_(:,:,:) = zero
   END IF
 END IF
 IF (paw_opt >= 3)  THEN
   allocate( gxfacs_(2,nlmn) )
   gxfacs_(:,:) = zero
   IF (choice > 1)  THEN
     allocate( dgxdtfacs_(2,ndgxdtfac,nlmn) )
     IF (ndgxdtfac > 0)        dgxdtfacs_(:,:,:) = zero
   END IF
 END IF
!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 IF (nthreads .EQ. 1)  THEN
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
! Scale gxfac with 4pi/sqr(omega).(-i)^l
       IF (paw_opt .NE. 3)  THEN
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfac_(1:cplex_fac,ilmn) = scale * gxfac(1:cplex_fac,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 1)                gxfac_(2,ilmn) = zero
           ELSE
             gxfac_(2,ilmn) = -scale * gxfac(1,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 2)  THEN
               gxfac_(1,ilmn) = scale * gxfac(2,ilmn,ia,ispinor)
             ELSE
               gxfac_(1,ilmn) = zero
             END IF
           END IF
         END DO

         IF (choice > 1)  THEN
           DO ilmn = 1, nlmn
             il = mod(indlmn(1,ilmn),4)
             parity = (mod(il,2) .EQ. 0)
             scale = wt
                      IF (il > 1)                         scale = -scale
             IF (parity)  THEN
               IF (cplex_fac .EQ. 2)  THEN
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn) = scale * dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfac_(ic,ii,ilmn) = scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn) = zero
                 END DO

               END IF
             ELSE
               IF (cplex_fac .EQ. 2)  THEN
                 DO ii = 1, ndgxdtfac
                   dgxdtfac_(2,ii,ilmn) = -scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn) = scale * dgxdtfac(2,ii,ilmn,ia,ispinor)
                 END DO

               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfac_(ic,ii,ilmn) = zero
                   IF (ic .EQ. 1)  THEN
                     dgxdtfac_(jc,ii,ilmn) = -scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   ELSE
                     dgxdtfac_(jc,ii,ilmn) = scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   END IF
                 END DO

               END IF
             END IF
           END DO

         END IF
       END IF
! Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       IF (paw_opt >= 3)  THEN
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfacs_(1:cplex,ilmn) = scale * gxfac_sij(1:cplex,ilmn,ia,ispinor)
             IF (cplex .EQ. 1)                gxfacs_(2,ilmn) = zero
           ELSE
             gxfacs_(2,ilmn) = -scale * gxfac_sij(1,ilmn,ia,ispinor)
             IF (cplex .EQ. 2)  THEN
               gxfacs_(1,ilmn) = scale * gxfac_sij(2,ilmn,ia,ispinor)
             ELSE
               gxfacs_(1,ilmn) = zero
             END IF
           END IF
         END DO

         IF (choice > 1)  THEN
           DO ilmn = 1, nlmn
             il = mod(indlmn(1,ilmn),4)
             parity = (mod(il,2) .EQ. 0)
             scale = wt
                      IF (il > 1)                         scale = -scale
             IF (parity)  THEN
               IF (cplex .EQ. 2)  THEN
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn) = scale * dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfacs_(ic,ii,ilmn) = scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn) = zero
                 END DO

               END IF
             ELSE
               IF (cplex .EQ. 2)  THEN
                 DO ii = 1, ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn) = -scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn) = scale * dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 END DO

               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfacs_(ic,ii,ilmn) = zero
                   IF (ic .EQ. 1)  THEN
                     dgxdtfacs_(jc,ii,ilmn) = -scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   ELSE
                     dgxdtfacs_(jc,ii,ilmn) = scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   END IF
                 END DO

               END IF
             END IF
           END DO

         END IF
       END IF
       allocate( ztab(npw) )
! Compute <g|Vnl|c> (or derivatives) for each plane wave:
       IF (paw_opt .NE. 3)  THEN
         ztab(:) = czero
! ------
! <g|Vnl|c>
         IF (choice .EQ. 1)  THEN
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 DO k = 1, npw
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 END DO

           END DO

         END IF
! ------
! derivative w.r.t. atm. pos
         IF (choice .EQ. 2)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
           END DO

           ztab(:) = two_pi * kpg(:,idir) * ztab(:)
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           END DO

         END IF
! ------
! derivative w.r.t. strain
         IF (choice .EQ. 3)  THEN
           IF (idir <= 3)  THEN
             DO ilmn = 1, nlmn
               ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn) - gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn) - &
               gxfac_(2,ilmn),kind=dp) - ffnl(:,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           ELSE
             DO ilmn = 1, nlmn
               ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) - &
               ffnl(:,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           END IF
         END IF
! ------
! full derivative w.r.t. k
         IF (choice .EQ. 5)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) + &
             ffnl(:,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           END DO

         END IF
! ------
! right derivative: <G|p>V<dp/dk|psi>
         IF (choice .EQ. 51)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           END DO

         END IF
! ------
! left derivative: <G|dp/dk>V<p|psi>
         IF (choice .EQ. 52)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           END DO

         END IF
! ------
! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
         IF (choice .EQ. 53)  THEN
! -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2 * idir - 1)
           fdb = ffnl_dir_dat(2 * idir)
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,fdf,ilmn) * cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) - &
             ffnl(:,fdb,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           END DO

         END IF
! ------
         ztab(:) = ztab(:) * cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         vect(1,(1 + ipwshft):(npw + ipwshft)) = vect(1,(1 + ipwshft):(npw + ipwshft)) + real(ztab(:))
         vect(2,(1 + ipwshft):(npw + ipwshft)) = vect(2,(1 + ipwshft):(npw + ipwshft)) + aimag(ztab(:))
       END IF
! Compute <g|S|c> (or derivatives) for each plane wave:
       IF (paw_opt >= 3)  THEN
         ztab(:) = czero
! ------
! <g|S|c>
         IF (choice .EQ. 1)  THEN
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 DO k = 1, npw
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 END DO

           END DO

         END IF
! ------
! derivative w.r.t. atm. pos
         IF (choice .EQ. 2)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
           END DO

           ztab(:) = two_pi * kpg(:,idir) * ztab(:)
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           END DO

         END IF
! ------
! derivative w.r.t. strain
         IF (choice .EQ. 3)  THEN
           IF (idir <= 3)  THEN
             DO ilmn = 1, nlmn
               ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn) - gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn) - &
               gxfacs_(2,ilmn),kind=dp) - ffnl(:,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           ELSE
             DO ilmn = 1, nlmn
               ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) - &
               ffnl(:,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           END IF
         END IF
! ------
! full derivative w.r.t. k
         IF (choice .EQ. 5)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) + &
             ffnl(:,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           END DO

         END IF
! ------
! right derivative: <G|p>V<dp/dk|psi>
         IF (choice .EQ. 51)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           END DO

         END IF
! ------
! left derivative: <G|dp/dk>V<p|psi>
         IF (choice .EQ. 52)  THEN
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           END DO

         END IF
! ------
! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
         IF (choice .EQ. 53)  THEN
! -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2 * idir - 1)
           fdb = ffnl_dir_dat(2 * idir)
           DO ilmn = 1, nlmn
             ztab(:) = ztab(:) + ffnl(:,fdf,ilmn) * cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) - &
             ffnl(:,fdb,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           END DO

         END IF
! ------
         ztab(:) = ztab(:) * cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         svect(1,(1 + ipwshft):(npw + ipwshft)) = svect(1,(1 + ipwshft):(npw + ipwshft)) + real(ztab(:))
         svect(2,(1 + ipwshft):(npw + ipwshft)) = svect(2,(1 + ipwshft):(npw + ipwshft)) + aimag(ztab(:))
       END IF
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 ELSE
! ==========================================================================
! ========== OPENMP VERSION ================================================
! ==========================================================================
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
! Scale gxfac with 4pi/sqr(omega).(-i)^l
       IF (paw_opt .NE. 3)  THEN
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfac_(1:cplex_fac,ilmn) = scale * gxfac(1:cplex_fac,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 1)                gxfac_(2,ilmn) = zero
           ELSE
             gxfac_(2,ilmn) = -scale * gxfac(1,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 2)  THEN
               gxfac_(1,ilmn) = scale * gxfac(2,ilmn,ia,ispinor)
             ELSE
               gxfac_(1,ilmn) = zero
             END IF
           END IF
         END DO

!$OMP END DO
         IF (choice > 1)  THEN
!$OMP DO
           DO ilmn = 1, nlmn
             il = mod(indlmn(1,ilmn),4)
             parity = (mod(il,2) .EQ. 0)
             scale = wt
                      IF (il > 1)                         scale = -scale
             IF (parity)  THEN
               IF (cplex_fac .EQ. 2)  THEN
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn) = scale * dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfac_(ic,ii,ilmn) = scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn) = zero
                 END DO

               END IF
             ELSE
               IF (cplex_fac .EQ. 2)  THEN
                 DO ii = 1, ndgxdtfac
                   dgxdtfac_(2,ii,ilmn) = -scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn) = scale * dgxdtfac(2,ii,ilmn,ia,ispinor)
                 END DO

               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfac_(ic,ii,ilmn) = zero
                   IF (ic .EQ. 1)  THEN
                     dgxdtfac_(jc,ii,ilmn) = -scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   ELSE
                     dgxdtfac_(jc,ii,ilmn) = scale * dgxdtfac(1,ii,ilmn,ia,ispinor)
                   END IF
                 END DO

               END IF
             END IF
           END DO

         END IF
!$OMP END DO
       END IF
!$OMP END PARALLEL
! Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       IF (paw_opt >= 3)  THEN
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfacs_(1:cplex,ilmn) = scale * gxfac_sij(1:cplex,ilmn,ia,ispinor)
             IF (cplex .EQ. 1)                gxfacs_(2,ilmn) = zero
           ELSE
             gxfacs_(2,ilmn) = -scale * gxfac_sij(1,ilmn,ia,ispinor)
             IF (cplex .EQ. 2)  THEN
               gxfacs_(1,ilmn) = scale * gxfac_sij(2,ilmn,ia,ispinor)
             ELSE
               gxfacs_(1,ilmn) = zero
             END IF
           END IF
         END DO

!$OMP END DO
         IF (choice > 1)  THEN
!$OMP DO
           DO ilmn = 1, nlmn
             il = mod(indlmn(1,ilmn),4)
             parity = (mod(il,2) .EQ. 0)
             scale = wt
                      IF (il > 1)                         scale = -scale
             IF (parity)  THEN
               IF (cplex .EQ. 2)  THEN
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn) = scale * dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfacs_(ic,ii,ilmn) = scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn) = zero
                 END DO

               END IF
             ELSE
               IF (cplex .EQ. 2)  THEN
                 DO ii = 1, ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn) = -scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn) = scale * dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 END DO

               ELSE
                 DO ii = 1, ndgxdtfac
                   ic = cplex_dgxdt(ii)
                   jc = 3 - ic
                   dgxdtfacs_(ic,ii,ilmn) = zero
                   IF (ic .EQ. 1)  THEN
                     dgxdtfacs_(jc,ii,ilmn) = -scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   ELSE
                     dgxdtfacs_(jc,ii,ilmn) = scale * dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   END IF
                 END DO

               END IF
             END IF
           END DO

         END IF
!$OMP END DO
       END IF
!$OMP END PARALLEL
       allocate( ztab(npw) )
! Compute <g|Vnl|c> (or derivatives) for each plane wave:
       IF (paw_opt .NE. 3)  THEN
!$OMP PARALLEL PRIVATE(ipw,ilmn,fdf,fdb)
! ------
! <g|Vnl|c>
         IF (choice .EQ. 1)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           END DO

         ELSE               IF (choice .EQ. 2)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
             END DO

             ztab(ipw) = two_pi * kpg(ipw,idir) * ztab(ipw)
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 3)  THEN
           IF (idir <= 3)  THEN
!$OMP DO
             DO ipw = 1, npw
               ztab(ipw) = czero
               DO ilmn = 1, nlmn
                 ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn) - gxfac_(1,ilmn),&
                 dgxdtfac_(2,1,ilmn) - gxfac_(2,ilmn),kind=dp) - ffnl(ipw,2,ilmn) * cmplx(gxfac_(1,ilmn),&
                 gxfac_(2,ilmn),kind=dp)
               END DO

             END DO

           ELSE
!$OMP END DO
!$OMP DO
             DO ipw = 1, npw
               ztab(ipw) = czero
               DO ilmn = 1, nlmn
                 ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) - &
                 ffnl(ipw,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               END DO

             END DO

           END IF
!$OMP END DO
              ELSE               IF (choice .EQ. 5)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) + &
               ffnl(ipw,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 51)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 52)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,2,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 53)  THEN
! -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
                fdf = ffnl_dir_dat(2 * idir - 1)
                fdb = ffnl_dir_dat(2 * idir)
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,fdf,ilmn) * cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) - &
               ffnl(ipw,fdb,ilmn) * cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             END DO

           END DO

              ELSE 
!$OMP END DO
!$OMP WORKSHARE
                ztab(:) = czero
         END IF
! ------
!$OMP DO
         DO ipw = 1, npw
           ztab(ipw) = ztab(ipw) * cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           vect(1,ipw + ipwshft) = vect(1,ipw + ipwshft) + real(ztab(ipw))
           vect(2,ipw + ipwshft) = vect(2,ipw + ipwshft) + aimag(ztab(ipw))
         END DO

       END IF
!$OMP END DO
!$OMP END PARALLEL
! Compute <g|S|c> (or derivatives) for each plane wave:
       IF (paw_opt >= 3)  THEN
!$OMP PARALLEL PRIVATE(ilmn,ipw,fdf,fdb)
! ------
! <g|S|c>
         IF (choice .EQ. 1)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           END DO

         ELSE               IF (choice .EQ. 2)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
             END DO

             ztab(ipw) = two_pi * kpg(ipw,idir) * ztab(ipw)
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 3)  THEN
           IF (idir <= 3)  THEN
!$OMP DO
             DO ipw = 1, npw
               ztab(ipw) = czero
               DO ilmn = 1, nlmn
                 ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn) - gxfacs_(1,ilmn),&
                 dgxdtfacs_(2,1,ilmn) - gxfacs_(2,ilmn),kind=dp) - ffnl(ipw,2,ilmn) * cmplx(gxfacs_(1,ilmn),&
                 gxfacs_(2,ilmn),kind=dp)
               END DO

             END DO

           ELSE
!$OMP END DO
!$OMP DO
             DO ipw = 1, npw
               ztab(ipw) = czero
               DO ilmn = 1, nlmn
                 ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) - &
                 ffnl(ipw,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               END DO

             END DO

           END IF
!$OMP END DO
              ELSE               IF (choice .EQ. 5)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) + &
               ffnl(ipw,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 51)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 52)  THEN
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,2,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           END DO

              ELSE               IF (choice .EQ. 53)  THEN
! <G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
                fdf = ffnl_dir_dat(2 * idir - 1)
                fdb = ffnl_dir_dat(2 * idir)
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,fdf,ilmn) * cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) - &
               ffnl(ipw,fdb,ilmn) * cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             END DO

           END DO

              ELSE 
!$OMP END DO
!$OMP WORKSHARE
                ztab(:) = czero
         END IF
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
         DO ipw = 1, npw
           ztab(ipw) = ztab(ipw) * cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           svect(1,ipw + ipwshft) = svect(1,ipw + ipwshft) + real(ztab(ipw))
           svect(2,ipw + ipwshft) = svect(2,ipw + ipwshft) + aimag(ztab(ipw))
         END DO

       END IF
!$OMP END DO
!$OMP END PARALLEL
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 END IF
! ==========================================================================
 IF (paw_opt .NE. 3)  THEN
   deallocate( gxfac_ )
   IF (choice > 1)  THEN
     deallocate( dgxdtfac_ )
   END IF
 END IF
 IF (paw_opt >= 3)  THEN
   deallocate( gxfacs_ )
   IF (choice > 1)  THEN
     deallocate( dgxdtfacs_ )
   END IF
 END IF
!Fake use of unused variable
! if (.false.) write(std_out,*) ipw
END SUBROUTINE opernlb_ylm

SUBROUTINE opernlb_ylm_ASSIST_choicee1_paw_opte3_cplexe2(cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,ph3d,svect,ucvol,vect)
 IMPLICIT NONE

  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER, PARAMETER :: dpc = kind(((1.0_dp, 1.0_dp)))
  COMPLEX(kind=dpc), PARAMETER :: czero = (0._dp, 0._dp)
  REAL(kind=dp), PARAMETER :: zero = 0._dp
  REAL(kind=dp), PARAMETER :: four = 4._dp
  REAL(kind=dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER :: two_pi = 2._dp * pi
  REAL(kind=dp), PARAMETER :: four_pi = four * pi
  INTEGER, INTENT(IN) :: cplex_fac
  INTEGER, INTENT(IN) :: dimffnl
  INTEGER, INTENT(IN) :: ia3
  INTEGER, INTENT(IN) :: idir
  INTEGER, INTENT(IN) :: matblk
  INTEGER, INTENT(IN) :: ndgxdtfac
  INTEGER, INTENT(IN) :: nincat
  INTEGER, INTENT(IN) :: nkpg
  INTEGER, INTENT(IN) :: nlmn
  INTEGER, INTENT(IN) :: npw
  INTEGER, INTENT(IN) :: nspinor
  REAL(kind=dp), INTENT(IN) :: ucvol
  INTEGER, DIMENSION(ndgxdtfac), INTENT(IN) :: cplex_dgxdt
  INTEGER, DIMENSION(6,nlmn), INTENT(IN) :: indlmn
  INTEGER, DIMENSION(5), INTENT(IN) :: nloalg
  REAL(kind=dp), DIMENSION(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor), INTENT(IN) :: dgxdtfac
  REAL(kind=dp), DIMENSION(2,ndgxdtfac,nlmn,nincat * (3 / 3),nspinor), INTENT(IN) :: dgxdtfac_sij
  REAL(kind=dp), DIMENSION(npw,dimffnl,nlmn), INTENT(IN) :: ffnl
  REAL(kind=dp), DIMENSION(cplex_fac,nlmn,nincat,nspinor), INTENT(IN) :: gxfac
  REAL(kind=dp), DIMENSION(2,nlmn,nincat,nspinor * (3 / 3)), INTENT(IN) :: gxfac_sij
  REAL(kind=dp), DIMENSION(npw,nkpg), INTENT(IN) :: kpg
  REAL(kind=dp), DIMENSION(2,npw,matblk), INTENT(IN) :: ph3d
  REAL(kind=dp), DIMENSION(2,npw * nspinor * (3 / 3)), INTENT(INOUT) :: svect
  REAL(kind=dp), DIMENSION(2,npw * nspinor), INTENT(INOUT) :: vect
  INTEGER :: fdb
  INTEGER :: fdf
  INTEGER :: ia
  INTEGER :: iaph3d
  INTEGER :: ic
  INTEGER :: ii
  INTEGER :: il
  INTEGER :: ilmn
  INTEGER :: ipw
  INTEGER :: ipwshft
  INTEGER :: ispinor
  INTEGER :: jc
  INTEGER :: nthreads
  INTEGER :: k
  REAL(kind=dp) :: scale
  REAL(kind=dp) :: wt
  LOGICAL :: parity
  INTEGER, DIMENSION(6), PARAMETER :: ffnl_dir_dat = (/3,4,4,2,2,3/)
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: dgxdtfac_
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: dgxdtfacs_
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxfac_
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxfacs_
  COMPLEX(kind=dpc), DIMENSION(:), ALLOCATABLE :: ztab
  INTEGER :: lt_var_k
  INTEGER :: lt_bound_npw
!Inits
  wt = four_pi / sqrt(ucvol)
  nthreads = 1
   allocate( gxfacs_(2,nlmn) )
  gxfacs_(:,:) = zero
!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 IF (nthreads .EQ. 1)  THEN
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfacs_(1:2,ilmn) = scale * gxfac_sij(1:2,ilmn,ia,ispinor)
           ELSE
             gxfacs_(2,ilmn) = -scale * gxfac_sij(1,ilmn,ia,ispinor)
  gxfacs_(1,ilmn) = scale * gxfac_sij(2,ilmn,ia,ispinor)
           END IF
         END DO

       allocate( ztab(npw) )
         ztab(:) = czero
  lt_bound_npw = (npw / 8) * 8
 DO lt_var_k = 1, lt_bound_npw, 8
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 DO k = lt_var_k, lt_var_k + (8 - 1)
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 END DO

           END DO

 END DO

 IF (lt_bound_npw < npw)  THEN
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 DO k = lt_bound_npw + 1, npw
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 END DO

           END DO

 END IF
! ------
         ztab(:) = ztab(:) * cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         svect(1,(1 + ipwshft):(npw + ipwshft)) = svect(1,(1 + ipwshft):(npw + ipwshft)) + real(ztab(:))
         svect(2,(1 + ipwshft):(npw + ipwshft)) = svect(2,(1 + ipwshft):(npw + ipwshft)) + aimag(ztab(:))
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 ELSE
! ==========================================================================
! ========== OPENMP VERSION ================================================
! ==========================================================================
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfacs_(1:2,ilmn) = scale * gxfac_sij(1:2,ilmn,ia,ispinor)
           ELSE
             gxfacs_(2,ilmn) = -scale * gxfac_sij(1,ilmn,ia,ispinor)
  gxfacs_(1,ilmn) = scale * gxfac_sij(2,ilmn,ia,ispinor)
           END IF
         END DO

!$OMP END DO
!$OMP END PARALLEL
       allocate( ztab(npw) )
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           END DO

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
         DO ipw = 1, npw
           ztab(ipw) = ztab(ipw) * cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           svect(1,ipw + ipwshft) = svect(1,ipw + ipwshft) + real(ztab(ipw))
           svect(2,ipw + ipwshft) = svect(2,ipw + ipwshft) + aimag(ztab(ipw))
         END DO

!$OMP END DO
!$OMP END PARALLEL
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 END IF
   deallocate( gxfacs_ )
END SUBROUTINE 

SUBROUTINE opernlb_ylm_ASSIST_choicee1_paw_opti3_cplexe2(cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,paw_opt,ph3d,svect,ucvol,vect)
 IMPLICIT NONE

  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER, PARAMETER :: dpc = kind(((1.0_dp, 1.0_dp)))
  COMPLEX(kind=dpc), PARAMETER :: czero = (0._dp, 0._dp)
  REAL(kind=dp), PARAMETER :: zero = 0._dp
  REAL(kind=dp), PARAMETER :: four = 4._dp
  REAL(kind=dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER :: two_pi = 2._dp * pi
  REAL(kind=dp), PARAMETER :: four_pi = four * pi
  INTEGER, INTENT(IN) :: cplex_fac
  INTEGER, INTENT(IN) :: dimffnl
  INTEGER, INTENT(IN) :: ia3
  INTEGER, INTENT(IN) :: idir
  INTEGER, INTENT(IN) :: matblk
  INTEGER, INTENT(IN) :: ndgxdtfac
  INTEGER, INTENT(IN) :: nincat
  INTEGER, INTENT(IN) :: nkpg
  INTEGER, INTENT(IN) :: nlmn
  INTEGER, INTENT(IN) :: npw
  INTEGER, INTENT(IN) :: nspinor
  INTEGER, INTENT(IN) :: paw_opt
  REAL(kind=dp), INTENT(IN) :: ucvol
  INTEGER, DIMENSION(ndgxdtfac), INTENT(IN) :: cplex_dgxdt
  INTEGER, DIMENSION(6,nlmn), INTENT(IN) :: indlmn
  INTEGER, DIMENSION(5), INTENT(IN) :: nloalg
  REAL(kind=dp), DIMENSION(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor), INTENT(IN) :: dgxdtfac
  REAL(kind=dp), DIMENSION(2,ndgxdtfac,nlmn,nincat * (paw_opt / 3),nspinor), INTENT(IN) :: dgxdtfac_sij
  REAL(kind=dp), DIMENSION(npw,dimffnl,nlmn), INTENT(IN) :: ffnl
  REAL(kind=dp), DIMENSION(cplex_fac,nlmn,nincat,nspinor), INTENT(IN) :: gxfac
  REAL(kind=dp), DIMENSION(2,nlmn,nincat,nspinor * (paw_opt / 3)), INTENT(IN) :: gxfac_sij
  REAL(kind=dp), DIMENSION(npw,nkpg), INTENT(IN) :: kpg
  REAL(kind=dp), DIMENSION(2,npw,matblk), INTENT(IN) :: ph3d
  REAL(kind=dp), DIMENSION(2,npw * nspinor * (paw_opt / 3)), INTENT(INOUT) :: svect
  REAL(kind=dp), DIMENSION(2,npw * nspinor), INTENT(INOUT) :: vect
  INTEGER :: fdb
  INTEGER :: fdf
  INTEGER :: ia
  INTEGER :: iaph3d
  INTEGER :: ic
  INTEGER :: ii
  INTEGER :: il
  INTEGER :: ilmn
  INTEGER :: ipw
  INTEGER :: ipwshft
  INTEGER :: ispinor
  INTEGER :: jc
  INTEGER :: nthreads
  INTEGER :: k
  REAL(kind=dp) :: scale
  REAL(kind=dp) :: wt
  LOGICAL :: parity
  INTEGER, DIMENSION(6), PARAMETER :: ffnl_dir_dat = (/3,4,4,2,2,3/)
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: dgxdtfac_
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: dgxdtfacs_
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxfac_
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxfacs_
  COMPLEX(kind=dpc), DIMENSION(:), ALLOCATABLE :: ztab
  INTEGER :: lt_var_k
  INTEGER :: lt_bound_npw
!Inits
  wt = four_pi / sqrt(ucvol)
  nthreads = 1
   allocate( gxfac_(2,nlmn) )
  gxfac_(:,:) = zero
!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 IF (nthreads .EQ. 1)  THEN
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfac_(1:cplex_fac,ilmn) = scale * gxfac(1:cplex_fac,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 1)                gxfac_(2,ilmn) = zero
           ELSE
             gxfac_(2,ilmn) = -scale * gxfac(1,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 2)  THEN
               gxfac_(1,ilmn) = scale * gxfac(2,ilmn,ia,ispinor)
             ELSE
               gxfac_(1,ilmn) = zero
             END IF
           END IF
         END DO

       allocate( ztab(npw) )
         ztab(:) = czero
  lt_bound_npw = (npw / 8) * 8
 DO lt_var_k = 1, lt_bound_npw, 8
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 DO k = lt_var_k, lt_var_k + (8 - 1)
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 END DO

           END DO

 END DO

 IF (lt_bound_npw < npw)  THEN
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 DO k = lt_bound_npw + 1, npw
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 END DO

           END DO

 END IF
! ------
         ztab(:) = ztab(:) * cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         vect(1,(1 + ipwshft):(npw + ipwshft)) = vect(1,(1 + ipwshft):(npw + ipwshft)) + real(ztab(:))
         vect(2,(1 + ipwshft):(npw + ipwshft)) = vect(2,(1 + ipwshft):(npw + ipwshft)) + aimag(ztab(:))
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 ELSE
! ==========================================================================
! ========== OPENMP VERSION ================================================
! ==========================================================================
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfac_(1:cplex_fac,ilmn) = scale * gxfac(1:cplex_fac,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 1)                gxfac_(2,ilmn) = zero
           ELSE
             gxfac_(2,ilmn) = -scale * gxfac(1,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 2)  THEN
               gxfac_(1,ilmn) = scale * gxfac(2,ilmn,ia,ispinor)
             ELSE
               gxfac_(1,ilmn) = zero
             END IF
           END IF
         END DO

!$OMP END DO
!$OMP END PARALLEL
       allocate( ztab(npw) )
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           END DO

! ------
!$OMP DO
         DO ipw = 1, npw
           ztab(ipw) = ztab(ipw) * cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           vect(1,ipw + ipwshft) = vect(1,ipw + ipwshft) + real(ztab(ipw))
           vect(2,ipw + ipwshft) = vect(2,ipw + ipwshft) + aimag(ztab(ipw))
         END DO

!$OMP END DO
!$OMP END PARALLEL
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 END IF
   deallocate( gxfac_ )
END SUBROUTINE 

SUBROUTINE opernlb_ylm_che1_pae3_cpe2(cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,&
  & gxfac_sij,ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,nspinor,paw_opt,ph3d,svect,ucvol,vect)
 IMPLICIT NONE

  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER, PARAMETER :: dpc = kind(((1.0_dp, 1.0_dp)))
  COMPLEX(kind=dpc), PARAMETER :: czero = (0._dp, 0._dp)
  REAL(kind=dp), PARAMETER :: zero = 0._dp
  REAL(kind=dp), PARAMETER :: four = 4._dp
  REAL(kind=dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER :: two_pi = 2._dp * pi
  REAL(kind=dp), PARAMETER :: four_pi = four * pi
  INTEGER, INTENT(IN) :: cplex_fac
  INTEGER, INTENT(IN) :: dimffnl
  INTEGER, INTENT(IN) :: ia3
  INTEGER, INTENT(IN) :: idir
  INTEGER, INTENT(IN) :: matblk
  INTEGER, INTENT(IN) :: ndgxdtfac
  INTEGER, INTENT(IN) :: nincat
  INTEGER, INTENT(IN) :: nkpg
  INTEGER, INTENT(IN) :: nlmn
  INTEGER, INTENT(IN) :: npw
  INTEGER, INTENT(IN) :: nspinor
  INTEGER, INTENT(IN) :: paw_opt
  REAL(kind=dp), INTENT(IN) :: ucvol
  INTEGER, DIMENSION(ndgxdtfac), INTENT(IN) :: cplex_dgxdt
  INTEGER, DIMENSION(6,nlmn), INTENT(IN) :: indlmn
  INTEGER, DIMENSION(5), INTENT(IN) :: nloalg
  REAL(kind=dp), DIMENSION(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor), INTENT(IN) :: dgxdtfac
  REAL(kind=dp), DIMENSION(2,ndgxdtfac,nlmn,nincat * (paw_opt / 3),nspinor), INTENT(IN) :: dgxdtfac_sij
  REAL(kind=dp), DIMENSION(npw,dimffnl,nlmn), INTENT(IN) :: ffnl
  REAL(kind=dp), DIMENSION(cplex_fac,nlmn,nincat,nspinor), INTENT(IN) :: gxfac
  REAL(kind=dp), DIMENSION(2,nlmn,nincat,nspinor * (paw_opt / 3)), INTENT(IN) :: gxfac_sij
  REAL(kind=dp), DIMENSION(npw,nkpg), INTENT(IN) :: kpg
  REAL(kind=dp), DIMENSION(2,npw,matblk), INTENT(IN) :: ph3d
  REAL(kind=dp), DIMENSION(2,npw * nspinor * (paw_opt / 3)), INTENT(INOUT) :: svect
  REAL(kind=dp), DIMENSION(2,npw * nspinor), INTENT(INOUT) :: vect
  INTEGER :: fdb
  INTEGER :: fdf
  INTEGER :: ia
  INTEGER :: iaph3d
  INTEGER :: ic
  INTEGER :: ii
  INTEGER :: il
  INTEGER :: ilmn
  INTEGER :: ipw
  INTEGER :: ipwshft
  INTEGER :: ispinor
  INTEGER :: jc
  INTEGER :: nthreads
  INTEGER :: k
  REAL(kind=dp) :: scale
  REAL(kind=dp) :: wt
  LOGICAL :: parity
  INTEGER, DIMENSION(6), PARAMETER :: ffnl_dir_dat = (/3,4,4,2,2,3/)
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: dgxdtfac_
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: dgxdtfacs_
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxfac_
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxfacs_
  COMPLEX(kind=dpc), DIMENSION(:), ALLOCATABLE :: ztab
  INTEGER :: lt_var_k
  INTEGER :: lt_bound_npw
!Inits
  wt = four_pi / sqrt(ucvol)
  nthreads = 1
   allocate( gxfac_(2,nlmn) )
  gxfac_(:,:) = zero
   allocate( gxfacs_(2,nlmn) )
  gxfacs_(:,:) = zero
!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 IF (nthreads .EQ. 1)  THEN
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfac_(1:cplex_fac,ilmn) = scale * gxfac(1:cplex_fac,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 1)                gxfac_(2,ilmn) = zero
           ELSE
             gxfac_(2,ilmn) = -scale * gxfac(1,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 2)  THEN
               gxfac_(1,ilmn) = scale * gxfac(2,ilmn,ia,ispinor)
             ELSE
               gxfac_(1,ilmn) = zero
             END IF
           END IF
         END DO

         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfacs_(1:2,ilmn) = scale * gxfac_sij(1:2,ilmn,ia,ispinor)
           ELSE
             gxfacs_(2,ilmn) = -scale * gxfac_sij(1,ilmn,ia,ispinor)
  gxfacs_(1,ilmn) = scale * gxfac_sij(2,ilmn,ia,ispinor)
           END IF
         END DO

       allocate( ztab(npw) )
         ztab(:) = czero
  lt_bound_npw = (npw / 8) * 8
 DO lt_var_k = 1, lt_bound_npw, 8
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 DO k = lt_var_k, lt_var_k + (8 - 1)
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 END DO

           END DO

 END DO

 IF (lt_bound_npw < npw)  THEN
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 DO k = lt_bound_npw + 1, npw
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
 END DO

           END DO

 END IF
! ------
         ztab(:) = ztab(:) * cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         vect(1,(1 + ipwshft):(npw + ipwshft)) = vect(1,(1 + ipwshft):(npw + ipwshft)) + real(ztab(:))
         vect(2,(1 + ipwshft):(npw + ipwshft)) = vect(2,(1 + ipwshft):(npw + ipwshft)) + aimag(ztab(:))
         ztab(:) = czero
  lt_bound_npw = (npw / 8) * 8
 DO lt_var_k = 1, lt_bound_npw, 8
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 DO k = lt_var_k, lt_var_k + (8 - 1)
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 END DO

           END DO

 END DO

 IF (lt_bound_npw < npw)  THEN
           DO ilmn = 1, nlmn
!             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 DO k = lt_bound_npw + 1, npw
  ztab(k) = ztab(k) + ffnl(k,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
 END DO

           END DO

 END IF
! ------
         ztab(:) = ztab(:) * cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         svect(1,(1 + ipwshft):(npw + ipwshft)) = svect(1,(1 + ipwshft):(npw + ipwshft)) + real(ztab(:))
         svect(2,(1 + ipwshft):(npw + ipwshft)) = svect(2,(1 + ipwshft):(npw + ipwshft)) + aimag(ztab(:))
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 ELSE
! ==========================================================================
! ========== OPENMP VERSION ================================================
! ==========================================================================
! Loop on spinorial components
   DO ispinor = 1, nspinor
     ipwshft = (ispinor - 1) * npw
! Loop on atoms (blocking)
     DO ia = 1, nincat
       iaph3d = ia
 IF (nloalg(1) > 0)                    iaph3d = ia + ia3 - 1
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfac_(1:cplex_fac,ilmn) = scale * gxfac(1:cplex_fac,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 1)                gxfac_(2,ilmn) = zero
           ELSE
             gxfac_(2,ilmn) = -scale * gxfac(1,ilmn,ia,ispinor)
             IF (cplex_fac .EQ. 2)  THEN
               gxfac_(1,ilmn) = scale * gxfac(2,ilmn,ia,ispinor)
             ELSE
               gxfac_(1,ilmn) = zero
             END IF
           END IF
         END DO

!$OMP END DO
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         DO ilmn = 1, nlmn
           il = mod(indlmn(1,ilmn),4)
           parity = (mod(il,2) .EQ. 0)
           scale = wt
                    IF (il > 1)                       scale = -scale
           IF (parity)  THEN
             gxfacs_(1:2,ilmn) = scale * gxfac_sij(1:2,ilmn,ia,ispinor)
           ELSE
             gxfacs_(2,ilmn) = -scale * gxfac_sij(1,ilmn,ia,ispinor)
  gxfacs_(1,ilmn) = scale * gxfac_sij(2,ilmn,ia,ispinor)
           END IF
         END DO

!$OMP END DO
!$OMP END PARALLEL
       allocate( ztab(npw) )
!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             END DO

           END DO

! ------
!$OMP DO
         DO ipw = 1, npw
           ztab(ipw) = ztab(ipw) * cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           vect(1,ipw + ipwshft) = vect(1,ipw + ipwshft) + real(ztab(ipw))
           vect(2,ipw + ipwshft) = vect(2,ipw + ipwshft) + aimag(ztab(ipw))
         END DO

!$OMP DO
           DO ipw = 1, npw
             ztab(ipw) = czero
             DO ilmn = 1, nlmn
               ztab(ipw) = ztab(ipw) + ffnl(ipw,1,ilmn) * cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             END DO

           END DO

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
         DO ipw = 1, npw
           ztab(ipw) = ztab(ipw) * cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           svect(1,ipw + ipwshft) = svect(1,ipw + ipwshft) + real(ztab(ipw))
           svect(2,ipw + ipwshft) = svect(2,ipw + ipwshft) + aimag(ztab(ipw))
         END DO

!$OMP END DO
!$OMP END PARALLEL
       deallocate( ztab )
! End loop on atoms
     END DO

! End loop on spinors
   END DO

 END IF
   deallocate( gxfac_ )
   deallocate( gxfacs_ )
END SUBROUTINE 

