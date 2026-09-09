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
subroutine te0545(option, nomte)
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
#include "asterfort/elrefv.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngfint.h"
#include "asterfort/nglgic.h"
#include "asterfort/ngvlog.h"
#include "asterfort/nmgvmb.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_GRAD_INCO, 3D_GRAD_VARI
!           D_PLAN_GRAD_INCO, D_PLAN_GRAD_VARI
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
    character(len=8), parameter :: fami = 'RIGI'
    character(len=8) :: typmod(2)
    character(len=16) :: defoComp
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: axi, matsym
    integer(kind=8) :: nnoQ, nnoL, npg, ndim, nddl, neps, lgpg
    integer(kind=8) :: jv_poids, jv_vfQ, jv_dfdeQ, jv_vfL, jv_dfdeL
    integer(kind=8) :: jvMaterc, icontm, ivarim, jvInstmr, jvInstpr, ideplm, ideplp
    integer(kind=8) :: ivectu, icontp, ivarip, imatns, imatuu, jvCarcri, ivarix, jvGeom, icoret
    integer(kind=8) :: iret, nnos, jv_ganoQ, jv_ganoL, itab(7)
    integer(kind=8) :: codret
    real(kind=8), allocatable :: b(:, :, :), w(:, :), ni2ldc(:, :)
    aster_logical :: lMatr, lVect, lSigm, lVari
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    matsym = ASTER_FALSE
    ivectu = 1
    icontp = 1
    ivarip = 1
    icoret = 1
    imatns = 1

! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))
    axi = typmod(1) .eq. 'AXIS'

! - Get parameters of element
    call elrefv('RIGI', ndim, nnoL, nnoQ, nnos, &
                npg, jv_poids, jv_vfL, jv_vfQ, jv_dfdeL, &
                jv_dfdeQ, jv_ganoL, jv_ganoQ)

! - PARAMETRES EN ENTREE ET DIMENSION
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)
    call tecach('OOO', 'PDEPLPR', 'L', iret, nval=2, itab=itab)
    nddl = itab(2)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nnoQ, jvGeom, materPara%lcsPara)

! - Get fields for non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Properties of behaviour
    defoComp = compor(DEFO)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              zr(jvInstmr), zr(jvInstpr), &
                              materPara, BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

!   NOMBRE DE VARIABLES INTERNES
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=itab)
    lgpg = max(itab(6), 1)*itab(7)

! - PARAMETRES EN SORTIE
    if (lMatr) then
        matsym = ASTER_FALSE
        call jevech('PMATUNS', 'E', imatns)
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', icoret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        zr(ivarip:ivarip-1+npg*lgpg) = zr(ivarix:ivarix-1+npg*lgpg)
    end if

! - CALCUL DES FORCES INTERIEURES ET MATRICES TANGENTES
    if (defoComp .eq. 'GDEF_LOG') then
        if (lteatt('INCO', 'C5GV')) then
            call nglgic(BEHInteg, &
                        option, typmod, ndim, nnoQ, &
                        nnoL, npg, nddl, jv_poids, zr(jv_vfQ), &
                        zr(jv_vfL), jv_dfdeQ, jv_dfdeL, zr(jvGeom), compor, &
                        lgpg, zr(jvCarcri), zr(jvInstmr), &
                        zr(jvInstpr), matsym, zr(ideplm), zr(ideplp), zr(icontm), &
                        zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), zr(imatns), &
                        lMatr, lVect, lSigm, lVari, codret)
        else
            call ngvlog(BEHInteg, &
                        option, typmod, ndim, nnoQ, &
                        nnoL, npg, nddl, jv_poids, zr(jv_vfQ), &
                        zr(jv_vfL), jv_dfdeQ, jv_dfdeL, zr(jvGeom), compor, &
                        lgpg, zr(jvCarcri), zr(jvInstmr), &
                        zr(jvInstpr), matsym, zr(ideplm), zr(ideplp), zr(icontm), &
                        zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), zr(imatns), &
                        lMatr, lVect, lSigm, lVari, codret)
        end if

    else if (defoComp(1:5) .eq. 'PETIT') then
        call nmgvmb(ndim, nnoQ, nnoL, npg, axi, &
                    zr(jvGeom), zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                    jv_poids, nddl, neps, b, w, &
                    ni2ldc)
        matsym = .not. (nint(zr(jvCarcri-1+CARCRI_MATRSYME)) .gt. 0)
        call ngfint(BEHInteg, &
                    option, typmod, ndim, nddl, neps, &
                    npg, w, b, compor, &
                    lgpg, zr(jvCarcri), zr(jvInstmr), &
                    zr(jvInstpr), zr(ideplm), zr(ideplp), ni2ldc, zr(icontm), &
                    zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), matsym, &
                    zr(imatuu), zr(imatns), lMatr, lVect, lSigm, &
                    codret)
        deallocate (b)
        deallocate (w)
        deallocate (ni2ldc)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (lSigm) then
        zi(icoret) = codret
    end if
!
end subroutine
