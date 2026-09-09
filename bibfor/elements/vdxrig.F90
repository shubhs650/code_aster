! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine vdxrig(plateCara, plateOrie, &
                  nomte, nodeCoor, matrRigi, nb1, &
                  indm, indf)
!
    use plate_type
    implicit none
!
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btkb.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrc.h"
#include "asterfort/matrkb.h"
#include "asterfort/moytem.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "blas/dscal.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: nodeCoor(3, 9)
    real(kind=8), intent(out) :: matrRigi(51, 51)
    integer(kind=8), intent(out) :: nb1
    integer(kind=8), intent(in) :: indm, indf
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: npge = 2
    integer(kind=8) :: nb2, nddle, npgsr, npgsn, itab(8), npg, nbLayer
    real(kind=8) :: tempMoye
    real(kind=8) :: vectBaseKpg(3, 3), vectTangKpg(2, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42), wmatcb(5, 42)
    real(kind=8) :: matrElas(5, 5), ktild(42, 42)
    real(kind=8) :: ctor, epais, kappa
    real(kind=8) :: ktildi(42, 42)
    integer(kind=8) :: i, inte, kpgsn, kpgsr, iret
    integer(kind=8) :: j, jcrf, kwgt, lzi, lzr, nddlet
    real(kind=8) :: coef

    real(kind=8) :: epsval(npge), ksi3s2
    blas_int :: b_incx, b_n
    data epsval/-0.577350269189626d0, 0.577350269189626d0/
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg)

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get plate parameters
    nbLayer = plateCara%nbLayer
    epais = plateCara%thick
    kappa = plateCara%shearCoef
    ctor = plateCara%coefRigiDRZ
!
    nddle = 5*nb1+2

! - Compute mean temperature (on all point and "sous-point" gauss)
    call moytem(fami, npg, 3*nbLayer, '+', tempMoye, iret)

    ktild = 0.d0
    kwgt = 0
    do inte = 1, npge
        ksi3s2 = epsval(inte)/2.d0

! ----- MEMBRANE ET CISAILLEMENT
        do kpgsr = 1, npgsr
            call mahsms(plateOrie, &
                        0, nb1, &
                        nodeCoor, ksi3s2, kpgsr, &
                        zr(lzr), epais, &
                        vectBaseKpg, vectTangKpg, &
                        hsfm, hss)
            call hsj1ms(epais, vectTangKpg, vectBaseKpg, hsfm, hss, &
                        hsj1m, hsj1s)
            call btdmsr(nb1, nb2, ksi3s2, kpgsr, zr(lzr), &
                        epais, plateOrie%vectTang, hsj1m, hsj1s, btdm, &
                        btds)
        end do

! ----- POUR L ENERGIE DE DEFORMATION DE MEMBRANE PAS DE CISAILLEMENT
        if (indm .eq. 1) then
            call r8inir(4*2*42, 0.d0, btds, 1)
        end if

! ----- POUR L ENERGIE DE DEFORMATION DE FLEXION
        if (indf .eq. 1) then
            call r8inir(4*3*42, 0.d0, btdm, 1)
            call r8inir(4*2*42, 0.d0, btds, 1)
        end if

! ----- MEMBRANE ET FLEXION
        do kpgsn = 1, npgsn
!
!     CALCUL DE BTDFN : F=FLEXION , N=NORMAL
!     ET DEFINITION DE WGT=PRODUIT DES POIDS ASSOCIES AUX PTS DE GAUSS
!                          (NORMAL) ET DU DETERMINANT DU JACOBIEN
            call mahsf(plateOrie, &
                       1, nb1, &
                       nodeCoor, ksi3s2, kpgsn, &
                       zr(lzr), epais, &
                       vectBaseKpg, vectTangKpg, &
                       hsf)
            call hsj1f(kpgsn, zr(lzr), epais, vectTangKpg, vectBaseKpg, &
                       hsf, kwgt, hsj1fx, wgt)
            call btdfn(1, nb1, nb2, ksi3s2, kpgsn, &
                       zr(lzr), epais, plateOrie%vectTang, hsj1fx, btdf)
!
!     CALCUL DE BTDMN, BTDSN
!     ET
!     FORMATION DE BTILD
!
!
!---- POUR L ENERGIE DE DEFORMATION DE MEMBRANE
!
            if (indm .eq. 1) then
                call r8inir(3*42, 0.d0, btdf, 1)
            end if
!
!
            call btdmsn(1, nb1, kpgsn, npgsr, zr(lzr), &
                        btdm, btdf, btds, btild)

! --------- Compute elastic matrix
            call matrc(plateOrie, vectBaseKpg, tempMoye, kappa, matrElas)
!
            b_n = to_blas_int(25)
            b_incx = to_blas_int(1)
            call dscal(b_n, wgt, matrElas, b_incx)
!
            call btkb(5, 42, nddle, matrElas, btild, &
                      wmatcb, ktildi)
!
            do i = 1, nddle
                do j = 1, nddle
                    ktild(i, j) = ktild(i, j)+ktildi(i, j)
                end do
            end do
        end do
    end do
!
!     EXPANSION DE LA MATRICE : AJOUTER DE LA ROTATION FICTIVE
!
!
    nddlet = 6*nb1+3
    call matrkb(nb1, 42, 51, nddlet, ktild, &
                ctor, matrRigi, coef)
!
!     AFFECTATION DU COEF POUR LA CONTRIBUTION DES ROTATIONS FICTIVES
!     POUR LE CALCUL NON LINEAIRE
!     (CETTE AFFECTATION N'A LIEU QUE DANS LE CAS OU ON PREND LA
!     MATRICE ELASTIQUE AU LIEU DE LA MATRICE TANGENTE)
!
    call tecach('NNO', 'PCACO3D', 'E', iret, nval=8, itab=itab)
    jcrf = itab(1)
    if (jcrf .ne. 0) zr(jcrf) = coef
!
end subroutine
