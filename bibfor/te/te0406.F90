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
subroutine te0406(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/btkb.h"
#include "asterfort/dxroep.h"
#include "asterfort/jacbm1.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matrn.h"
#include "asterfort/pmavec.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "asterfort/transp.h"
#include "asterfort/utmess.h"
#include "asterfort/vectgt.h"
#include "blas/ddot.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DES OBJETS ELEMENTS FINIS EN DYNAMIQUE LINEAIRE COQUE_3D
!
!     OPTIONS   :  MASS_MECA      MATRICE DE MASSE COHERENTE
!                  M_GAMMA        FORCE NODALE D INERTIE
!                  ECIN_ELEM ENERGIE CINETIQUE D UN MODE PROPRE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvGeom
    integer(kind=8) :: lzi, lzr
    integer(kind=8) :: nb1, nb2
    integer(kind=8) :: intsn, npgsn
    integer(kind=8) :: inte
    real(kind=8) :: rho, epais, ctor
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: jm1(3, 3), detj
    integer(kind=8) :: i, j, iret
    integer(kind=8) :: jd
    integer(kind=8) :: kompt
    integer(kind=8) :: imatuu, iacce, ivect
    integer(kind=8) :: jener, jfreq, iu, iv
    real(kind=8) :: mas(2601), masu(51), masv(51)
    real(kind=8) :: mantn(2601)
    real(kind=8) :: bid33(3, 3)
    real(kind=8) :: matn(3, 51), matnt(51, 3)
    real(kind=8) :: ksi3s2
    integer(kind=8), parameter :: npge = 2
    real(kind=8), parameter :: epsval(npge) = (/-1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0)/)
    real(kind=8) :: xmin
    integer(kind=8) :: in, icompo
    integer(kind=8) :: ii, jj
    integer(kind=8) :: imas
    real(kind=8) :: lam0(3, 3)
    real(kind=8) :: masrg(3, 3)
    real(kind=8) :: masrl(3, 3)
    real(kind=8) :: mnn
    blas_int :: b_incx, b_incy, b_n
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
!
!---- TEST D'EXISTENCE "COMPOR"
!
    call tecach('NNO', 'PCOMPOR', 'L', iret, iad=icompo)
    if (icompo .ne. 0) then
        if (zk16(icompo+2) .eq. 'GROT_GDEP') then
            call utmess('F', 'ELEMENTS3_91')
        end if
    end if

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get plate parameters
    call dxroep(plateCara, rho, epais)
    ctor = plateCara%coefRigiDRZ

    call r8inir(51*51, 0.d0, mas, 1)
    call r8inir(51*51, 0.d0, mantn, 1)
    do inte = 1, npge
! ----- COORDONNEE ISOPARAMETRIQUE SUR L EPAISSEUR  DIVISEE PAR DEUX
        ksi3s2 = epsval(inte)/2.d0

        do intsn = 1, npgsn
! --------- Compute local basis at current integration point
            call vectgt(plateOrie, 1, nb1, &
                        zr(jvGeom), ksi3s2, intsn, &
                        epais, zr(lzr), &
                        vectBaseKpg, vectTangKpg)

!---------- CALCUL DE DETJ
            call jacbm1(epais, vectTangKpg, vectBaseKpg, bid33, jm1, &
                        detj)

!---------  MATRICE N
            call matrn(nb1, nb2, zr(lzr), ksi3s2, epais, &
                       intsn, plateOrie%vectNorm, matn)

!---------- TRANSPOSE DE MATN
            call transp(matn, 3, 3, 6*nb1+3, matnt, &
                        6*nb1+3)

!---------- PRODUIT MANT * MATN
            call promat(matnt, 6*nb1+3, 6*nb1+3, 3, matn, &
                        3, 3, 6*nb1+3, mantn)

!---------- INTEGRATION NUMERIQUE
            do j = 1, 6*nb1+3
                do i = 1, 6*nb1+3
                    jd = (6*nb1+3)*(j-1)+i
                    mas(jd) = mas(jd)+(rho*mantn(jd)*zr(lzr-1+127+intsn-1)*detj*1.d0)
                end do
            end do
        end do
    end do
!
!
    xmin = 1.d0/r8prem()
!
!---- EN CHAQUE NOEUD
!
    do in = 1, nb2
!
!------- ON CONSTRUIT LAMBDA0
!
        do ii = 1, 3
            lam0(ii, 1) = plateOrie%vectTang(in, 1, ii)
            lam0(ii, 2) = plateOrie%vectTang(in, 2, ii)
            lam0(ii, 3) = plateOrie%vectNorm(in, ii)
        end do
!
!------- ON CONSTRUIT MASRG
!
        if (in .le. nb1) then
!
!-------------- NOEUDS DE SERENDIP
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*(in-1)+jj+3
                    i = 6*(in-1)+ii+3
                    imas = (6*nb1+3)*(j-1)+i
                    masrg(ii, jj) = mas(imas)
                end do
            end do
!
        else
!
!-------------- SUPERNOEUD
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*nb1+jj
                    i = 6*nb1+ii
                    imas = (6*nb1+3)*(j-1)+i
                    masrg(ii, jj) = mas(imas)
                end do
            end do
!
        end if
!
!------- ROTATION DE MASRG : LOCALES --> GLOBALES
!
!        MASRL =  ( LAMBDA0 )   * SIGMT * ( LAMBDA0 ) T
!
        call btkb(3, 3, 3, masrg, lam0, &
                  bid33, masrl)
!
!------- ON COMPARE LES DEUX PREMIERS TERMES DIAGONAUX DE MASRL
!
        if (masrl(1, 1) .lt. xmin) xmin = masrl(1, 1)
        if (masrl(2, 2) .lt. xmin) xmin = masrl(2, 2)
!
    end do
!
    mnn = ctor*xmin
    do in = 1, nb2
        if (in .le. nb1) then
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*(in-1)+jj+3
                    i = 6*(in-1)+ii+3
                    mas((6*nb1+3)*(j-1)+i) = mas((6*nb1+3)*(j-1)+i)+ &
                                             mnn* &
                                             plateOrie%vectNorm(in, ii)*plateOrie%vectNorm(in, jj)
                end do
            end do
        else
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*nb1+jj
                    i = 6*nb1+ii
                    mas((6*nb1+3)*(j-1)+i) = mas((6*nb1+3)*(j-1)+i)+ &
                                             mnn* &
                                             plateOrie%vectNorm(in, ii)*plateOrie%vectNorm(in, jj)
                end do
            end do
        end if
    end do

    if (option .eq. 'MASS_MECA') then
        call jevech('PMATUUR', 'E', imatuu)
        kompt = 0
        do j = 1, 6*nb1+3
            do i = 1, j
                kompt = kompt+1
                zr(imatuu-1+kompt) = mas((6*nb1+3)*(j-1)+i)
            end do
        end do

    else if (option .eq. 'M_GAMMA') then
        call jevech('PACCELR', 'L', iacce)
        call jevech('PVECTUR', 'E', ivect)
        call pmavec('ZERO', 6*nb1+3, mas, zr(iacce), zr(ivect))

    else if (option .eq. 'ECIN_ELEM') then
        call jevech('PENERCR', 'E', jener)
        call tecach('ONO', 'PVITESR', 'L', iret, iad=iv)
        if (iret .eq. 0) then
            call r8inir(51, 0.d0, masv, 1)
            call pmavec('ZERO', 6*nb1+3, mas, zr(iv), masv)
            b_n = to_blas_int(6*nb1+3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(jener) = 5.d-1*ddot(b_n, zr(iv), b_incx, masv, b_incy)
        else
            call tecach('ONO', 'PDEPLAR', 'L', iret, iad=iu)
            if (iret .eq. 0) then
                call jevech('POMEGA2', 'L', jfreq)
                call r8inir(51, 0.d0, masu, 1)
                call pmavec('ZERO', 6*nb1+3, mas, zr(iu), masu)
                b_n = to_blas_int(6*nb1+3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(jener) = ddot(b_n, zr(iu), b_incx, masu, b_incy)
!
!--------- VITESSE = OMEGA * MODE
!
                zr(jener) = 0.5d0*zr(jfreq)*zr(jener)
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
!------- ENERGIE DE MEMBRANE = ENERGIE TOTALE
!        ENERGIE DE FLEXION  = ENERGIE TOTALE
        call r8inir(2, zr(jener), zr(jener+1), 1)

    else
        ASSERT(.false.)
    end if
!
end subroutine
