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
subroutine te0322(option, nomte)
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
#include "asterfort/ejinit.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmfihm.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_JOINT_HYME
!           PLAN_JOINT_HYME
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
    integer(kind=8) :: ndim, nno1, nno2, npg, nddl, ntrou
    integer(kind=8) :: iw, ivf1, ivf2, idf2
    integer(kind=8) :: jvGeom, jvMaterc, jvCarcri, iddlm, iddld
    integer(kind=8) :: icontm, icontp, ivect, imatr, iu(3, 16), ip(8)
    integer(kind=8) :: ivarim, ivarip, jtab(7), iret, iinstm, iinstp
    integer(kind=8) :: lgpg
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: typmod(2), lielrf(10)
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: codret
    integer(kind=8) :: jv_codret
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    ivarip = 1
    icontp = 1
    ivect = 1

! - Get element parameters
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami=fami, ndim=ndim, nno=nno1, &
                     jvf=ivf1)
    call elrefe_info(elrefe=lielrf(1), fami=fami, ndim=ndim, nno=nno2, &
                     npg=npg, jpoids=iw, jvf=ivf2, jdfde=idf2)

! LA DIMENSION DE L'ESPACE EST CELLE DE L'ELEM DE REF SURFACIQUE PLUS 1
    ndim = ndim+1
!
! NB DE DDL = NDIM PAR NOEUD DE DEPL + UN PAR NOEUD DE PRES
    nddl = ndim*2*nno1+nno2
!
! DECALAGE D'INDICE POUR LES ELEMENTS DE JOINT
    call ejinit(nomte, iu, ip)

! - Type of finite element
    if (ndim == 3) then
        typmod(1) = '3D'
    elseif (ndim == 2) then
        typmod(1) = 'PLAN'
    else
        ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
    end if
    if (lteatt('TYPMOD2', 'EJ_HYME')) then
        typmod(2) = 'EJ_HYME'
    elseif (lteatt('TYPMOD2', 'ELEMJOIN')) then
        typmod(2) = 'ELEMJOIN'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PDEPLMR', 'L', iddlm)
    call jevech('PDEPLPR', 'L', iddld)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCONTMR', 'L', icontm)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Get fields for non-linear behaviour
    call jevech('PCARCRI', 'L', jvCarcri)
    call jevech('PCOMPOR', 'L', vk16=compor)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              zr(iinstm), zr(iinstp), &
                              materPara, BEHInteg)

! - Get output fields
    if (lMatr) then
        call jevech('PMATUNS', 'E', imatr)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivect)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if

! - Compute
    call nmfihm(BEHInteg, &
                ndim, nddl, nno1, nno2, npg, &
                lgpg, iw, zr(iw), zr(ivf1), zr(ivf2), &
                idf2, zr(idf2), option, zr(jvGeom), &
                zr(iddlm), zr(iddld), iu, ip, zr(icontm), &
                zr(icontp), zr(ivect), zr(imatr), zr(ivarim), zr(ivarip), &
                zr(iinstm), zr(iinstp), zr(jvCarcri), compor, typmod, &
                lVect, lMatr, lSigm, codret)

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', jv_codret)
        zi(jv_codret) = codret
    end if
!
end subroutine
