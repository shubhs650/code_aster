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
subroutine te0206(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gedisc.h"
#include "asterfort/jevech.h"
#include "asterfort/nmfi3d.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_JOINT
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: nno, npg, nddl
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: jvGeom, jvMaterc, jvCarcri, jvDeplmr, jvDeplpr, icoret
    integer(kind=8) :: icontm, icontp, ivect, imatr
    integer(kind=8) :: ivarim, ivarip, jtab(7), iret, iinstm, iinstp
    integer(kind=8) :: lgpg, codret
!     COORDONNEES POINT DE GAUSS + POIDS : X,Y,Z,W => 1ER INDICE
    real(kind=8) :: coopg(4, 4)
    character(len=8), parameter :: typmod(2) = (/'3D      ', 'ELEMJOIN'/)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relaComp
    aster_logical :: matsym
    aster_logical :: lVect, lMatr, lVari, lSigm
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    ivarip = 1
    icoret = 1
    icontp = 1
    ivect = 1
    icoret = 1
    imatr = 1
!
! -  FONCTIONS DE FORMES ET POINTS DE GAUSS : ATTENTION CELA CORRESPOND
!    ICI AUX FONCTIONS DE FORMES 2D DES FACES DES MAILLES JOINT 3D
!    PAR EXEMPLE FONCTION DE FORME DU QUAD4 POUR LES HEXA8.
!
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
    nddl = 6*nno
!
    ASSERT(nno .le. 4)
    ASSERT(npg .le. 4)

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', jvDeplmr)
    call jevech('PDEPLPR', 'L', jvDeplpr)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Get behaviour parameters
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)
    relaComp = compor(RELA_NAME)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              zr(iinstm), zr(iinstp), &
                              materPara, BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - CALCUL DES COORDONNEES DES POINTS DE GAUSS, POIDS=0
    call gedisc(3, nno, npg, zr(ivf), zr(jvGeom), coopg)

! - Get output fields
    if (lMatr) then
        matsym = ASTER_TRUE
        if (relaComp .eq. 'JOINT_MECA_RUPT') matsym = ASTER_FALSE
        if (relaComp .eq. 'JOINT_MECA_FROT') matsym = ASTER_FALSE
        if (matsym) then
            call jevech('PMATUUR', 'E', imatr)
        else
            call jevech('PMATUNS', 'E', imatr)
        end if
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivect)
    end if
!
    call nmfi3d(BEHInteg, typmod, &
                nno, nddl, npg, lgpg, zr(ipoids), &
                zr(ivf), zr(idfde), option, zr(jvGeom), &
                zr(jvDeplmr), zr(jvDeplpr), zr(icontm), zr(icontp), zr(ivect), &
                zr(imatr), zr(ivarim), zr(ivarip), zr(jvCarcri), compor, &
                matsym, coopg, zr(iinstm), zr(iinstp), lMatr, lVect, lSigm, &
                codret)

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', icoret)
        zi(icoret) = codret
    end if
!
end subroutine
