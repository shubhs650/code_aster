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
!
subroutine te0247(option, nomte)
!
    use Behaviour_module, only: behaviourOption
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matela.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/nmpoel.h"
#include "asterfort/porea1.h"
#include "asterfort/porigi.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ptkg00.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/verifm.h"
#include "asterfort/get_value_mode_local.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES OPTIONS FULL_MECA OU RAPH_MECA POUR
!     LES ELEMENTS DE POUTRE 'MECA_POU_D_E/D_T'
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        'MECA_POU_D_E' : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!        'MECA_POU_D_T' : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: igeom, jvMaterc, iorien, nd, nk
    integer(kind=8) :: iinstm, iinstp, icontm, ideplm, ideplp, imatuu
    integer(kind=8) :: ivectu, icontp, itype, nno, nc, ivarim, ivarip, itemp, i
    integer(kind=8) :: jcret, iretm, iretp
    integer(kind=8) :: npg, ndimel, nnoel, nnosel
    integer(kind=8) :: istrxm, istrxp, ldep, codret
    parameter(nno=2, nc=6, nd=nc*nno, nk=nd*(nd+1)/2)
    real(kind=8) :: e, nu, em, num
    real(kind=8) :: a, xiy, xiz, alfay, alfaz, xjx, ez, ey, xfly, xflz
    real(kind=8) :: a2, xiy2, xiz2, alfay2, alfaz2, xjx2, xl
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relaComp, defoComp
    real(kind=8) :: pgl(3, 3), fl(nd), klv(nk)
    real(kind=8) :: tempm, tempp
    real(kind=8) :: epsthe
    real(kind=8) :: sigma(nd), rgeom(nk), gamma, angp(3)
    aster_logical :: reactu
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: retp(2), iret
    real(kind=8) :: valr(2)
    character(len=8) :: valp(2)
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_cara = 17
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8), parameter :: noms_cara(nb_cara) = (/'A1  ', 'IY1 ', 'IZ1 ', &
                                                          'AY1 ', 'AZ1 ', 'EY1 ', &
                                                          'EZ1 ', 'JX1 ', 'A2  ', &
                                                          'IY2 ', 'IZ2 ', 'AY2 ', &
                                                          'AZ2 ', 'EY2 ', 'EZ2 ', &
                                                          'JX2 ', 'TVAR'/)
! --------------------------------------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PCAORIE', 'L', iorien)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PDEPLMR', 'L', ideplm)
!
!   la presence du champ de deplacement a l instant t+ devrait etre conditionne par l'option
!   (avec rigi_meca_tang ca n a pas de sens).
!   cependant ce champ est initialise a 0 par la routine nmmatr.
    call jevech('PDEPLPR', 'L', ideplp)
!   POINT DE GAUSS DE L'ELEMENT
    call elrefe_info(fami='RIGI', ndim=ndimel, nno=nnoel, nnos=nnosel, npg=npg)
    ASSERT((npg .eq. 2) .or. (npg .eq. 3))

! - Get fields for non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)

! - Properties of behaviour
    relaComp = compor(RELA_NAME)
    defoComp = compor(DEFO)
    reactu = defoComp .eq. 'GROT_GDEP'

! - Some checks
    if (defoComp .ne. 'PETIT' .and. defoComp .ne. 'GROT_GDEP') then
        call utmess('F', 'POUTRE0_40', sk=defoComp)
    end if
    if (relaComp .ne. 'ELAS') then
        call utmess('F', 'POUTRE0_61', sk=relaComp)
    end if

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

    if (option .eq. 'RIGI_MECA_TANG') then
        lVect = ASTER_FALSE
        lSigm = ASTER_FALSE
        lVari = ASTER_FALSE
    end if

!   parametres en sortie
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
!
    call matrot(zr(iorien), pgl)
    xl = lonele()
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
!
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    alfay = vale_cara(4)
    alfaz = vale_cara(5)
    xjx = vale_cara(8)
    a2 = vale_cara(9)
    xiy2 = vale_cara(10)
    xiz2 = vale_cara(11)
    alfay2 = vale_cara(12)
    alfaz2 = vale_cara(13)
    xjx2 = vale_cara(16)
    ey = (vale_cara(6)+vale_cara(14))/2.d0
    ez = (vale_cara(7)+vale_cara(15))/2.d0
    itype = nint(vale_cara(17))
!

    if (reactu) then
!       recuperation du 3eme angle nautique au temps t-
        call jevech('PSTRXMR', 'L', istrxm)
        gamma = zr(istrxm+3-1)
!       calcul de pgl,xl et angp
        call porea1(nno, nc, zr(ideplm), zr(ideplp), zr(igeom), gamma, lVect, pgl, xl, angp)
!       sauvegarde des angles nautiques
        if (lVect) then
            call jevech('PSTRXPR', 'E', istrxp)
            zr(istrxp+1-1) = angp(1)
            zr(istrxp+2-1) = angp(2)
            zr(istrxp+3-1) = angp(3)
        end if
!
    end if

! - Get average temperatures
    call moytem('RIGI', npg, 1, '+', tempp, iretp)
    call moytem('RIGI', npg, 1, '-', tempm, iretm)
    itemp = 0
    if ((iretp+iretm) .eq. 0) itemp = 1

! - Recuperation des caracteristiques elastiques
    call matela(zi(jvMaterc), ' ', itemp, tempp, e, nu)
    call matela(zi(jvMaterc), ' ', itemp, tempm, em, num)
    call verifm('RIGI', npg, 1, 'T', zi(jvMaterc), epsthe, iret)

!  - Calcul des matrices elementaires
    call porigi(nomte, e, nu, xl, klv)
!
    if (option .eq. 'RAPH_MECA' .or. option .eq. 'FULL_MECA') then
        if ((nomte .ne. 'MECA_POU_D_E') .and. (itemp .ne. 0) .and. (nu .ne. num)) then
            call utmess('A', 'POUTRE0_59')
        end if
        call nmpoel(npg, klv, xl, nno, nc, pgl, zr(ideplp), &
                    epsthe, e, em, zr(icontm), fl, zr(icontp))
    end if
!
    if (lMatr) then
!       calcul de la matrice de rigidite geometrique
        if (reactu) then
            ! Avec GROT_GDEP pas possible
            valp(1:2) = ['C_FLEX_Y', 'C_FLEX_Z']
            call get_value_mode_local('PCAARPO', valp, valr, iret, retpara_=retp)
            xfly = 1.0; xflz = 1.0
            if (retp(1) .eq. 0) xfly = valr(1)
            if (retp(2) .eq. 0) xflz = valr(2)
            if ((abs(xfly-1.0) .gt. r8miem()) .or. &
                (abs(xflz-1.0) .gt. r8miem())) then
                call utmess('F', 'POUTRE0_64', nr=2, valr=[xfly, xflz])
            end if
            if (option .eq. 'FULL_MECA') then
                ldep = icontp
            else
                ldep = icontm
            end if
            if (npg .eq. 2) then
                do i = 1, nc
                    sigma(i) = zr(ldep+i-1)
                    sigma(i+nc) = zr(ldep+nc+i-1)
                end do
            else
                do i = 1, nc
                    sigma(i) = zr(ldep+i-1)
                    sigma(i+nc) = zr(ldep+nc+nc+i-1)
                end do
            end if
            call r8inir(nk, 0.0d0, rgeom, 1)
            call ptkg00(sigma, a, a2, xiz, xiz2, xiy, xiy2, xl, ey, ez, rgeom)
            klv = klv+rgeom
        end if
    end if

!   passage du repere local au repere global
    if (lMatr) then
        call utpslg(nno, nc, pgl, klv, zr(imatuu))
    end if
    if (lVect) then
        call utpvlg(nno, nc, pgl, fl, zr(ivectu))
    end if
    if (lSigm) then
        zi(jcret) = codret
    end if
!
end subroutine
