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
! aslint: disable=W0413
! => real zero (init by calcul.F90)
!
subroutine forngr(plateCara, plateOrie, &
                  option, nomte)
!
    use plate_type
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/btsig.h"
#include "asterfort/cosiro.h"
#include "asterfort/jacbm1.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/jm1dn1.h"
#include "asterfort/jm1dn2.h"
#include "asterfort/matbmn.h"
#include "asterfort/matbmr.h"
#include "asterfort/matbsr.h"
#include "asterfort/matbsu.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vectgt.h"
#include "asterfort/vectpe.h"
#include "asterfort/vectrn.h"
#include "blas/daxpy.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     FONCTION  :  FORC_NODA DES COQUE_3D
!                  GEOMETRIQUE AVEC GRANDES ROTATIONS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3
    real(kind=8) :: bid33(3, 3)
    integer(kind=8) :: i, j, in, ii, nval, k1, iret, itab(7)
    real(kind=8) :: stild(5)
    integer(kind=8) :: jvGeom, icontm, ivectu
    integer(kind=8) :: lzi, lzr
    integer(kind=8) :: nb1, nb2
    integer(kind=8) :: inte, intsr, intsn
    real(kind=8) :: eptot
    integer(kind=8) :: npgsr, npgsn
    real(kind=8) :: vecnph(9, 3)
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: jm1(3, 3), detj
    real(kind=8) :: jdn1ri(9, 51), jdn1rc(9, 51)
    real(kind=8) :: jdn1ni(9, 51), jdn1nc(9, 51)
    real(kind=8) :: jdn2rc(9, 51)
    real(kind=8) :: jdn2nc(9, 51)
    real(kind=8) :: ksi3s2
    integer(kind=8) :: nbLayer, nbsp
    integer(kind=8) :: iLayer
    real(kind=8) :: zic, zmin, hLayer, coef
    real(kind=8) :: sigref
    integer(kind=8) :: jvDisp
    real(kind=8) :: b1su(5, 51), b2su(5, 51)
    real(kind=8) :: b1src(2, 51, 4)
    real(kind=8) :: b2src(2, 51, 4)
    real(kind=8) :: b1mnc(3, 51), b1mni(3, 51)
    real(kind=8) :: b2mnc(3, 51), b2mni(3, 51)
    real(kind=8) :: b1mri(3, 51, 4)
    real(kind=8) :: b2mri(3, 51, 4)
    real(kind=8) :: dudxri(9), dudxni(9)
    real(kind=8) :: dudxrc(9), dudxnc(9)
    real(kind=8) :: vectDisp(8, 3), vectRota(9, 3)
    real(kind=8) :: vecpe(51)
    real(kind=8) :: sigtmp(5), ftemp(51), effint(51)
    character(len=16) :: kmess(2)
    real(kind=8) :: blam(9, 3, 3)
    blas_int :: b_incx, b_incy, b_n
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!

! - Get plate parameters
    nbLayer = plateCara%nbLayer
    ASSERT(nbLayer .ge. 1)
    eptot = plateCara%thick
    zmin = -eptot/2.d0
    hLayer = eptot/nbLayer

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get stress (in good frame !)
    if (option .eq. 'FORC_NODA') then
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, itab=itab)
        nbsp = itab(7)
        if (nbsp .ne. npge*nbLayer) then
            call utmess('F', 'ELEMENTS_4')
        end if
        call cosiro(plateCara, plateOrie, &
                    'PSIEFR', 'L', 'UI', 'G', &
                    icontm)
    else if (option .eq. 'REFE_FORC_NODA') then
        call refe%Init(nomte)
        sigref = refe%GetRef('SIGM')
        call refe%Check()
    end if

! - Get displacements
    if (option .eq. "FORC_NODA") then
        call jevech('PDEPLAR', 'L', jvDisp)
    else
        call jevech('PDEPLMR', 'L', jvDisp)
    end if

! - DEPLACEMENT TOTAL AUX NOEUDS DE SERENDIP
    vectDisp = 0.d0
    do in = 1, nb1
        do ii = 1, 3
            vectDisp(in, ii) = zr(jvDisp-1+6*(in-1)+ii)
        end do
    end do

! - ROTATION TOTALE AUX NOEUDS
    vectRota = 0.d0
    do in = 1, nb1
        do ii = 1, 3
            vectRota(in, ii) = zr(jvDisp-1+6*(in-1)+ii+3)
        end do
    end do
    do ii = 1, 3
        vectRota(nb2, ii) = zr(jvDisp-1+6*nb1+ii)
    end do

! - TRANSFORMEES NORMALES ET MATRICES DE ROTATION AUX NOEUDS
    call vectrn(nb2, plateOrie%vectTang, plateOrie%vectNorm, vectRota, vecnph, &
                blam)

! - VECTEUR PE DES VARIABLES NODALES TOTALES GENERALISEES
    call vectpe(nb1, nb2, vectDisp, plateOrie%vectNorm, vecnph, &
                vecpe)
!
!______________________________________________________________________
!
!---- INITIALISATION DES OPERATEURS DE DEFORMATION A EXTRAPOLER
!
!---- MEMBRANE REDUIT INCOMPLET
!
    call r8inir(3*51*4, 0.d0, b1mri, 1)
!
    call r8inir(3*51*4, 0.d0, b2mri, 1)
!
!---- SHEAR    REDUIT   COMPLET
!
    call r8inir(2*51*4, 0.d0, b1src, 1)
!
    call r8inir(2*51*4, 0.d0, b2src, 1)

    call jevech('PVECTUR', 'E', ivectu)
    ftemp = 0.d0
    do iLayer = 1, nbLayer
        do inte = 1, npge
! --------- POSITION SUR L EPAISSEUR ET POIDS D INTEGRATION
            if (inte .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
                coef = 1.d0/3.d0
            else if (inte .eq. 2) then
                zic = zmin+hLayer/2.d0+(iLayer-1)*hLayer
                coef = 4.d0/3.d0
            else
                zic = zmin+hLayer+(iLayer-1)*hLayer
                coef = 1.d0/3.d0
            end if

!---------- COORDONNEE ISOP.  SUR L EPAISSEUR  DIVISEE PAR DEUX
            ksi3s2 = zic/hLayer

            do intsr = 1, npgsr
! ------------- Compute local base at integration point
                call vectgt(plateOrie, 0, nb1, &
                            zr(jvGeom), ksi3s2, intsr, &
                            hLayer, zr(lzr), &
                            vectBaseKpg, vectTangKpg)
!
                call jacbm1(hLayer, vectTangKpg, vectBaseKpg, bid33, jm1, &
                            detj)
!
!------------- J1DN1RI ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 0 INCOMPLET
                call jm1dn1(0, 0, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsr, jm1, jdn1ri)
!
!------------- CALCUL DE    DUDXRI ( 9 ) REDUIT INCOMPLET
!
                call promat(jdn1ri, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxri)
!
!+++++++++++++ B1MRI ( 3 , 51 , 4 ) MEMBRANE REDUIT INCOMPLET
!              B2MRI ( 3 , 51 , 4 )
!
                call matbmr(nb1, vectBaseKpg, dudxri, intsr, jdn1ri, &
                            b1mri, b2mri)
!
!------------- J1DN1RC ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 1 COMPLET
!
                call jm1dn1(0, 1, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsr, jm1, jdn1rc)
!
!------------- CALCUL DE    DUDXRC ( 9 ) REDUIT COMPLET
!
                call promat(jdn1rc, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxrc)
!
!------------- J1DN2RC ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 1 COMPLET
!
                call jm1dn2(0, 1, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsr, vecnph, jm1, &
                            jdn2rc)
!
!+++++++++++++ B1SRC ( 2 , 51 , 4 ) SHEAR REDUIT COMPLET
!              B2SRC ( 2 , 51 , 4 )
!
                call matbsr(nb1, vectBaseKpg, dudxrc, intsr, jdn1rc, &
                            jdn2rc, b1src, b2src)
            end do

            do intsn = 1, npgsn
! ------------- Compute local base at integration point
                call vectgt(plateOrie, 1, nb1, &
                            zr(jvGeom), ksi3s2, intsn, &
                            hLayer, zr(lzr), &
                            vectBaseKpg, vectTangKpg)
!
                call jacbm1(hLayer, vectTangKpg, vectBaseKpg, bid33, jm1, &
                            detj)
!
!------------- J1DN1NC ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 1 COMPLET
!
                call jm1dn1(1, 1, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsn, jm1, jdn1nc)
!
!------------- CALCUL DE     DUDXNC ( 9 ) NORMAL COMPLET
!
                call promat(jdn1nc, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxnc)
!
!------------- J1DN2NC ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 1 COMPLET
!
                call jm1dn2(1, 1, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsn, vecnph, jm1, &
                            jdn2nc)
!
!+++++++++++++ B1MNC ( 3 , 51 ) MEMBRANE NORMAL COMPLET
!              B2MNC ( 3 , 51 )
!
                call matbmn(nb1, vectBaseKpg, dudxnc, jdn1nc, jdn2nc, &
                            b1mnc, b2mnc)
!
!------------- J1DN1NI ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 0 INCOMPLET
!
                call jm1dn1(1, 0, nb1, nb2, zr(lzr), &
                            hLayer, ksi3s2, intsn, jm1, jdn1ni)
!
!------------- CALCUL DE     DUDXNI ( 9 ) NORMAL INCOMPLET
!
                call promat(jdn1ni, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxni)
!
!+++++++++++++ B1MNI ( 3 , 51 ) MEMBRANE NORMAL INCOMPLET
!              B2MNI ( 3 , 51 )
!
                call matbmn(nb1, vectBaseKpg, dudxni, jdn1ni, jdn1ni, &
                            b1mni, b2mni)
!
!============= B1SU ( 5 , 51 ) SUBSTITUTION TOTAL
!              B2SU ( 5 , 51 ) SUBSTITUTION DIFFERENTIEL
!
                call matbsu(nb1, zr(lzr), npgsr, intsn, b1mnc, &
                            b2mnc, b1mni, b2mni, b1mri, b2mri, &
                            b1src, b2src, b1su, b2su)
!
                if (option .eq. 'FORC_NODA') then
!
!------- CONTRAINTES DE CAUCHY = PK2 AUX POINTS DE GAUSS
!
                    k1 = 6*((intsn-1)*npge*nbLayer+(iLayer-1)*npge+inte-1)
                    stild(1) = zr(icontm-1+k1+1)
                    stild(2) = zr(icontm-1+k1+2)
                    stild(3) = zr(icontm-1+k1+4)
                    stild(4) = zr(icontm-1+k1+5)
                    stild(5) = zr(icontm-1+k1+6)
!
!------------- FINT ( 6 * NB1 + 3 )  =     INTEGRALE  DE
!              ( B2SU ( 5 , 6 * NB1 + 3 ) ) T * STILD ( 5 ) *
!              POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btsig(6*nb1+3, 5, zr(lzr-1+127+intsn-1)*detj*coef, b2su, stild, &
                               zr(ivectu))
!
!------------- VARIABLES INTERNES INACTIVES COMPORTEMENT NON PLASTIQUE
!
                else if (option .eq. 'REFE_FORC_NODA') then
                    call r8inir(5, 0.d0, sigtmp, 1)
                    call r8inir(51, 0.d0, effint, 1)
!
                    do i = 1, 5
                        sigtmp(i) = sigref
                        call btsig(6*nb1+3, 5, zr(lzr-1+127+intsn-1)*detj*coef, b2su, sigtmp, &
                                   effint)
                        sigtmp(i) = 0.d0
                        do j = 1, 51
                            ftemp(j) = ftemp(j)+abs(effint(j))
                        end do
                    end do
                end if
            end do
        end do
    end do
!
!      ON PREND LA VALEUR MOYENNE DES FORCES NODALES DE REFERENCE
!
    if (option .eq. 'REFE_FORC_NODA') then
        nval = nbLayer*npge*npgsn*5
        b_n = to_blas_int(51)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/nval, ftemp, b_incx, zr(ivectu), b_incy)
        do j = 1, 51
            if (zr(ivectu+j-1) .eq. 0.) then
                kmess(1) = 'COQUE3D'
                kmess(2) = 'SIGM_REFE'
                call utmess('F', 'MECANONLINE5_59', nk=2, valk=kmess)
            end if
        end do
    end if
!
end subroutine
