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
subroutine te0393(option, nomte)
!
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "jeveux.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gdfint.h"
#include "asterfort/gdjrg0.h"
#include "asterfort/jevech.h"
#include "asterfort/marota.h"
#include "asterfort/promat.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES FORCES NODALES DE MECA_POU_D_T_GD
!                          OPTION : 'FORC_NODA' OU 'REFE_FORC_NODA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    character(len=8) :: elrefe
    real(kind=8) :: en(3, 2), enprim(3, 2), fint(6, 3), y0(3), x00(3, 3)
    real(kind=8) :: x0k(3, 3), qik(3, 3), x0pg(3), qig(3), rot(3, 3), rot0(3, 3)
    real(kind=8) :: rotabs(3, 3), gn(3), gm(3), pn(3), pm(3), zero, un, unsurj
    real(kind=8) :: pjacob, ajacob
!
    real(kind=8) :: forref, momref
    integer(kind=8) :: nno, nc, ino, i, ndim, nnos, npg, ipoids, ivf, idfdk, jgano, kp
    integer(kind=8) :: ne, ic, kc, k0, k1, ico
    integer(kind=8) :: ivectu, jvDisp, igeom, jvSief, lorien, jefint
    integer(kind=8) :: ifint
!
    parameter(zero=0.0d0, un=1.0d0)
    type(RESI_REFE):: refe
! ----------------------------------------------------------------------
!
!
    if (option .eq. 'REFE_FORC_NODA') then
        nno = 2
        nc = 6
        call refe%Init(nomte)
        forref = refe%GetRef('EFFORT')
        momref = refe%GetRef('MOMENT')
        call refe%Check()

        call jevech('PVECTUR', 'E', ivectu)
        do ino = 1, nno
            do i = 1, 3
                zr(ivectu+(ino-1)*nc+i-1) = forref
            end do
            do i = 4, nc
                zr(ivectu+(ino-1)*nc+i-1) = momref
            end do
        end do
!
    else if (option .eq. 'FORC_NODA') then
        call elref1(elrefe)
!
!        PARAMETRES EN ENTREE
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
        ico = 0
        do kp = 1, npg
            do ne = 1, nno
                ico = ico+1
                en(ne, kp) = zr(ivf-1+ico)
                enprim(ne, kp) = zr(idfdk-1+ico)
            end do
        end do
!
!        CALL JEVECH('PTEMPPR','L',ITEMPR)
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PDEPLAR', 'L', jvDisp)
        call jevech('PSIEFR', 'L', jvSief)
!
!        --- RECUPERATION DES ORIENTATIONS INITIALES Y0(1), Y0(2), Y0(3)
        call jevech('PCAORIE', 'L', lorien)
        y0(1) = zr(lorien)
        y0(2) = zr(lorien+1)
        y0(3) = zr(lorien+2)
!
!        PARAMETRES EN SORTIE
        call jevech('PVECTUR', 'E', jefint)
        do ne = 1, nno
            do kc = 1, 6
                fint(kc, ne) = zero
            end do
        end do
!
!       DO 21 NE=1,NNO
!         TEMPN(NE) = ZR(ITEMPR-1+NE)
!21      CONTINUE
!
        k0 = igeom-1
        k1 = jvDisp-1
!
        do ne = 1, nno
            do kc = 1, 3
                k0 = k0+1
                k1 = k1+1
                x00(kc, ne) = zr(k0)
                x0k(kc, ne) = zr(k0)+zr(k1)
            end do
            do kc = 1, 3
                k1 = k1+1
                qik(kc, ne) = zr(k1)
            end do
        end do
!
!        BOUCLE SUR LES POINTS DE GAUSS
        do kp = 1, npg
            call gdjrg0(kp, nno, enprim, x00, y0, &
                        ajacob, rot0)
            pjacob = zr(ipoids-1+kp)*ajacob
!
            do ic = 1, 3
                x0pg(ic) = zero
                qig(ic) = zero
            end do
!           TEMPG = ZERO
            unsurj = un/ajacob
            do ic = 1, 3
                do ne = 1, nno
                    x0pg(ic) = x0pg(ic)+unsurj*enprim(ne, kp)*x0k(ic, ne)
                    qig(ic) = qig(ic)+en(ne, kp)*qik(ic, ne)
                end do
            end do
!         DO 45 NE=1,NNO
!          TEMPG = TEMPG + EN(NE,KP)*TEMPN(NE)
!45       CONTINUE
!
            call marota(qig, rot)
            call promat(rot, 3, 3, 3, rot0, &
                        3, 3, 3, rotabs)
            do ic = 1, 3
                gn(ic) = zr(jvSief-1+6*(kp-1)+ic)
                gm(ic) = zr(jvSief+2+6*(kp-1)+ic)
            end do
            call promat(rotabs, 3, 3, 3, gn, &
                        3, 3, 1, pn)
            call promat(rotabs, 3, 3, 3, gm, &
                        3, 3, 1, pm)
!
            call gdfint(kp, nno, ajacob, pjacob, en, &
                        enprim, x0pg, pn, pm, fint)
!
        end do
!        FIN DE BOUCLE SUR LES POINTS DE GAUSS
        ifint = jefint-1
        do ne = 1, nno
            do kc = 1, 6
                ifint = ifint+1
                zr(ifint) = fint(kc, ne)
            end do
        end do
    end if
!
end subroutine
