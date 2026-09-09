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
subroutine pk2cau(plateOrie, &
                  nbLayer, eptot, &
                  nomte, ncmp, &
                  pk2, sigma)
!
    use plate_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/btkb.h"
#include "asterfort/jacbm1.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/jm1dn1.h"
#include "asterfort/promat.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "asterfort/vectgt.h"
#include "asterfort/vectpe.h"
#include "asterfort/vectrn.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    integer(kind=8), intent(in) :: nbLayer
    real(kind=8), intent(in) :: eptot
    character(len=16), intent(in) :: nomte
    integer(kind=8), intent(in) :: ncmp
    real(kind=8), intent(in) :: pk2(ncmp, *)
    real(kind=8), intent(out) :: sigma(ncmp, *)
!
! --------------------------------------------------------------------------------------------------
!
!      PK2CAU  -- CALCUL DES CONTAINTES DE CAUCHY A PARTIR DES
!                 CONTRAINTES DE PIOLA-KIRCHHOFF DE SECONDE ESPECE
!                 A PARTIR DE LA FORMULE :
!
!            SIGMA = (1/DET[F])*([F]*[PK2]*[F]T)
!             OU [F] EST LA MATRICE DU GRADIENT DES DEFORMATIONS
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMTE          IN     K16      NOM DU TYPE D'ELEMENT
!    NCMP           IN     I        NOMBRE DE COMPOSANTES DU TENSEUR
!                                   DES CONTRAINTES
!    PK2(NCMP,1)    IN     R        TENSEUR DES CONTRAINTES
!                                   DE PIOLA-KIRCHHOFF DE SECONDE ESPECE
!    SIGMA(NCMP,1)  VAR    R        TENSEUR DES CONTRAINTES DE CAUCHY
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3, nbinco = 51
    real(kind=8), parameter :: un = 1.d0, deux = 2.d0
    integer(kind=8) :: i, iLayer, jvDisp, jvGeom, ii
    integer(kind=8) :: in, inte, intsn, iret, j, kpgs, lzi
    integer(kind=8) :: lzr, nb1, nb2, npgsn
    real(kind=8) :: cof11, cof21, cof31, detf, detfm1, detj
    real(kind=8) :: hLayer, zic, zmin
    real(kind=8) :: vectDisp(8, 3), vectRota(9, 3), vecnph(9, 3)
    real(kind=8) :: vectBaseInvKpg(3, 3), jm1(3, 3)
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: vecpe(nbinco), blam(9, 3, 3), bid33(3, 3)
    real(kind=8) :: xab(3, 3), dudx(3), dudy(3), dudz(3)
    real(kind=8) :: jdn1nc(9, nbinco), dudxnc(9), sigmKpgGlob(3, 3)
    real(kind=8) :: ft(3, 3), sigmKpgLoca(3, 3), pk2KpgLoca(3, 3), pk2KpgGlob(3, 3)
    real(kind=8) :: ksi3s2

!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbLayer .ge. 1)
    zmin = -eptot/deux
    hLayer = eptot/nbLayer

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi+1-1)
    nb2 = zi(lzi+2-1)
    npgsn = zi(lzi+4-1)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get displacements
    call tecach('NNO', 'PDEPLAR', 'L', iret, iad=jvDisp)
    if (iret .ne. 0) then
        call tecach('NNO', 'PDEPPLU', 'L', iret, iad=jvDisp)
        ASSERT(iret .eq. 0)
    end if

! - Separation of displacements and rotations
    vectDisp = 0.d0
    vectRota = 0.d0
    do in = 1, nb1
        do ii = 1, 3
            vectDisp(in, ii) = zr(jvDisp+6*(in-1)+ii-1)
            vectRota(in, ii) = zr(jvDisp+6*(in-1)+ii+3-1)
        end do
    end do
    do ii = 1, 3
        vectRota(nb2, ii) = zr(jvDisp+6*nb1+ii-1)
    end do

! - DETERMINATION AUX NOEUDS DES VECTEURS VECNPH QUI SONT LA
! - TRANSFORMEE APRES DEFORMATION DES VECTEURS VECTN NORMAUX
! - AU PLAN MOYEN INITIAL ET DES MATRICES DE ROTATION BLAM FAISANT
! - PASSER DES VECTEURS VECTN AUX VECTEURS VECNPH :
    call vectrn(nb2, plateOrie%vectTang, plateOrie%vectNorm, vectRota, vecnph, &
                blam)

! - DETERMINATION DU VECTEUR DE DEPLACEMENT AUX NOEUDS VECPE
! - DEFINI PAR VECPE = <U V W (NPHI-N)_X (NPHI-N)_Y (NPHI-N)_Z>
! - OU U, V, W SONT LES 3 DDLS DE TRANSLATION
! - NPHI EST LE VECTEUR VECNPH ET N LE VECTEUR VECTN :
    call vectpe(nb1, nb2, vectDisp, plateOrie%vectNorm, vecnph, &
                vecpe)

    kpgs = 0
    do iLayer = 1, nbLayer
        do inte = 1, npge
            if (inte .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
            else if (inte .eq. 2) then
                zic = zmin+hLayer/deux+(iLayer-1)*hLayer
            else if (inte .eq. 3) then
                zic = zmin+hLayer+(iLayer-1)*hLayer
            end if
            ksi3s2 = zic/hLayer
            do intsn = 1, npgsn
                kpgs = kpgs+1
! ------------- Compute local base at integration point
                call vectgt(plateOrie, 1, nb1, &
                            zr(jvGeom), ksi3s2, intsn, &
                            hLayer, zr(lzr), &
                            vectBaseKpg, vectTangKpg)
                call jacbm1(hLayer, vectTangKpg, vectBaseKpg, bid33, jm1, &
                            detj)

! ------------- CALCUL DU VECTEUR JDN1NC QUI EST < DU/DQSI> (I.E.
! ------------- <DU/DQSI1,DU/DQSI2,DU/DQSI3,DV/DQSI1,DV/DQSI2,DV/DQSI3,
! ------------- DW/DQSI1,DW/DQSI2,DW/DQSI3> )
                call jm1dn1(1, 1, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsn, jm1, jdn1nc)

! ------------- CALCUL DU VECTEUR DUDXNC QUI EST < DU/DX> (I.E.
! ------------- <DU/DX,DU/DY,DU/DZ,DV/DX,DV/DY,DV/DZ,DW/DX,DW/DY,DW/DZ> )
                call promat(jdn1nc, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxnc)
!
                do i = 1, 3
                    dudx(i) = dudxnc(1+3*(i-1))
                    dudy(i) = dudxnc(2+3*(i-1))
                    dudz(i) = dudxnc(3+3*(i-1))
                end do

! ------------- CONSTRUCTION DE LA MATRICE [F] DU GRADIENT DES
! ------------- DEFORMATIONS AU POINT D'INTEGRATION COURANT.
!                          | 1 0 0 |   | DU/DX DU/DY DU/DZ |
!                    [F] = | 0 1 0 | + | DV/DX DV/DY DV/DZ |
!                          | 0 0 1 |   | DW/DX DW/DY DW/DZ |
! ------------- PAR COMMODITE, ON UTILISE PLUTOT [FT] , LA MATRICE
! ------------- TRANSPOSEE DE [F]
                do i = 1, 3
                    ft(1, i) = dudx(i)
                    ft(2, i) = dudy(i)
                    ft(3, i) = dudz(i)
                end do
                ft(1, 1) = ft(1, 1)+un
                ft(2, 2) = ft(2, 2)+un
                ft(3, 3) = ft(3, 3)+un

! ------------- CALCUL DU DETERMINANT DE [F] ( = DET [FT] )
                cof11 = ft(2, 2)*ft(3, 3)-ft(2, 3)*ft(3, 2)
                cof21 = ft(3, 1)*ft(2, 3)-ft(2, 1)*ft(3, 3)
                cof31 = ft(2, 1)*ft(3, 2)-ft(3, 1)*ft(2, 2)
                detf = ft(1, 1)*cof11+ft(1, 2)*cof21+ft(1, 3)*cof31
                detfm1 = un/(detf+r8prem())

! ------------- CONSTRUCTION DU TENSEUR DES CONTRAINTES PK2
                pk2KpgLoca(1, 1) = pk2(1, kpgs)
                pk2KpgLoca(2, 2) = pk2(2, kpgs)
                pk2KpgLoca(3, 3) = pk2(3, kpgs)
                pk2KpgLoca(1, 2) = pk2(4, kpgs)
                pk2KpgLoca(1, 3) = pk2(5, kpgs)
                pk2KpgLoca(2, 3) = pk2(6, kpgs)
                pk2KpgLoca(2, 1) = pk2KpgLoca(1, 2)
                pk2KpgLoca(3, 1) = pk2KpgLoca(1, 3)
                pk2KpgLoca(3, 2) = pk2KpgLoca(2, 3)

! ------------- PK2 local => global
                call btkb(3, 3, 3, pk2KpgLoca, vectBaseKpg, &
                          bid33, pk2KpgGlob)

! ------------- Cauchy: [SIGMAG] = [F]*[PK2]*[FT]
                call utbtab('ZERO', 3, 3, pk2KpgGlob, ft, xab, sigmKpgGlob)

! ------------- Compute inverse matrix (global => local)
                do i = 1, 3
                    do j = 1, 3
                        vectBaseInvKpg(i, j) = vectBaseKpg(j, i)
                    end do
                end do

! ------------- Cauchy global => local
                call btkb(3, 3, 3, sigmKpgGlob, vectBaseInvKpg, &
                          bid33, sigmKpgLoca)

                sigma(1, kpgs) = sigmKpgLoca(1, 1)*detfm1
                sigma(2, kpgs) = sigmKpgLoca(2, 2)*detfm1
                sigma(3, kpgs) = sigmKpgLoca(3, 3)*detfm1
                sigma(4, kpgs) = sigmKpgLoca(1, 2)*detfm1
                sigma(5, kpgs) = sigmKpgLoca(1, 3)*detfm1
                sigma(6, kpgs) = sigmKpgLoca(2, 3)*detfm1
            end do
        end do
    end do
!
end subroutine
