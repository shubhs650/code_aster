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
subroutine fornpd(plateCara, plateOrie, &
                  option, nomte)
!
    use plate_type
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/cosiro.h"
#include "asterfort/epseff.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/tecach.h"
#include "asterfort/trndgl.h"
#include "asterfort/trnflg.h"
#include "asterfort/utmess.h"
#include "asterfort/vexpan.h"
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
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3
    character(len=10) :: elasKeyword
    integer(kind=8) :: i, ib, iLayer, inte, intsn, intsr, j, k1, kpgs, kwgt, itab(7), iret
    integer(kind=8) :: icontm, jvDisp, jvMaterc, ivectu, jvGeom, lzi, lzr
    integer(kind=8) :: nb1, nb2, npgsr, npgsn, nbLayer, nval, nbsp
    real(kind=8) :: vectBase(9, 3, 3)
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: eptot
    real(kind=8) :: rotfm(9)
    real(kind=8) :: deplm(42), effint(42), vecl(48), vecll(51)
    real(kind=8) :: sgmtd(5)
    real(kind=8) :: ksi3s2
    real(kind=8) :: sigtmp(5), ftemp(40), sigref
    real(kind=8) :: zero, zic, zmin, coef, hLayer
!
    character(len=16) :: kmess(2)
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

! - Access to material parameters
    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
    if (elasKeyword .ne. 'ELAS' .and. elasKeyword .ne. 'ELAS_ORTH') then
        call utmess('F', 'ELEMENTS_44', sk=elasKeyword)
    end if

! - Change coordinates of displacements/rotations
    call trndgl(nb2, plateOrie%vectNorm, plateOrie%vectTang, zr(jvDisp), &
                deplm, rotfm)

    effint = 0.d0
    ftemp = 0.d0
    kwgt = 0
    kpgs = 0
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
                call mahsms(plateOrie, &
                            0, nb1, &
                            zr(jvGeom), ksi3s2, intsr, &
                            zr(lzr), hLayer, &
                            vectBaseKpg, vectTangKpg, &
                            hsfm, hss)
                call hsj1ms(hLayer, vectTangKpg, vectBaseKpg, hsfm, hss, &
                            hsj1m, hsj1s)
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            hLayer, plateOrie%vectTang, hsj1m, hsj1s, btdm, &
                            btds)
            end do

            do intsn = 1, npgsn
                call mahsf(plateOrie, &
                           1, nb1, &
                           zr(jvGeom), ksi3s2, intsn, &
                           zr(lzr), hLayer, &
                           vectBaseKpg, vectTangKpg, &
                           hsf)
                call hsj1f(intsn, zr(lzr), hLayer, vectTangKpg, vectBaseKpg, &
                           hsf, kwgt, hsj1fx, wgt)
                wgt = coef*wgt
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), hLayer, plateOrie%vectTang, hsj1fx, btdf)
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbLayer+(iLayer-1)*npge+inte-1)
                if (option .eq. 'FORC_NODA') then
                    sgmtd(1) = zr(icontm-1+k1+1)
                    sgmtd(2) = zr(icontm-1+k1+2)
                    sgmtd(3) = zr(icontm-1+k1+4)
                    sgmtd(4) = zr(icontm-1+k1+5)
                    sgmtd(5) = zr(icontm-1+k1+6)
                    call epseff('EFFORI', nb1, [0.d0], btild, sgmtd, &
                                [0.d0], wgt, effint)

                else if (option .eq. 'REFE_FORC_NODA') then
!
!      CALCUL DES FORCES NODALES DE REFERENCE
!      EN AFFECTANT LA VALEUR SIGM_REFE A CHAQUE CMP SUCCESSIVEMENT
!      POUR CHAQUE POINT D'INTEGRATION
!
                    call r8inir(5, 0.d0, sigtmp, 1)
!
                    do i = 1, 5
                        sigtmp(i) = sigref
                        call epseff('EFFORI', nb1, [0.d0], btild, sigtmp, &
                                    [0.d0], wgt, effint)
                        sigtmp(i) = 0.d0
                        do j = 1, nb1*5
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
        b_n = to_blas_int(nb1*5)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/nval, ftemp, b_incx, effint, b_incy)
    end if

! - EXPANSION DU CHAMP
    call vexpan(nb1, effint, vecl)
!
    do i = 1, 6*nb1
        vecll(i) = vecl(i)
    end do
    vecll(6*nb1+1) = effint(5*nb1+1)
    vecll(6*nb1+2) = effint(5*nb1+2)
!        VECLL(6*NB1+3)=0.D0
!
!     ICI PAS DE CONTRIBUTION DES DDL DE LA ROTATION FICTIVE DANS EFFINT
!
    zero = 0.d0
    do i = 1, nb1
        vecll(6*i) = zero*rotfm(i)
    end do
    i = nb2
    vecll(6*nb1+3) = zero*rotfm(nb2)

! - Fuse tangents and normal in same object
    do ib = 1, nb2
        vectBase(ib, 1:2, 1:3) = plateOrie%vectTang(ib, 1:2, 1:3)
        vectBase(ib, 3, 1:3) = plateOrie%vectNorm(ib, 1:3)
    end do
!
    call jevech('PVECTUR', 'E', ivectu)
!
    call trnflg(nb2, vectBase, vecll, zr(ivectu))
!
    if (option .eq. 'REFE_FORC_NODA') then
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
