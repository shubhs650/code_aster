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
subroutine te0359(option, nomte)

    use Behaviour_type
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/assert.h"
#include "asterfort/eimatb.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngpipe.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
! Elementary computation
!
! Elements: 3D_INTERFACE
!           PLAN_INTERFACE, AXIS_INTERFACE
!
! Options: PILO_PRED_ELAS
!
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter:: ntrou_max = 10
! --------------------------------------------------------------------------------------------------
    aster_logical :: axi
    character(len=8) :: typmod(2), lielrf(ntrou_max), attrib
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: ntrou, ndim_fe, ndim_sp, nno2, nno1, npg, lgpg, jtab(7), nddl, neps
    integer(kind=8) :: jv_poids, jv_vff2, jv_vff1, jv_dfde2, jv_dfde1
    integer(kind=8) :: jv_geom, jv_materc, jv_carcri
    integer(kind=8) :: jv_contm, jv_varim, jv_copil, jv_borne, jv_dtau, jv_typilo
    integer(kind=8) :: jv_deplm, jv_ddepl, jv_depl0, jv_depl1, iret
    real(kind=8) :: instam, instap
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
    real(kind=8), allocatable:: b(:, :, :), w(:, :), ni2ldc(:, :)
! --------------------------------------------------------------------------------------------------

! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), vattr_missing=' ')
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp
    axi = lteatt('AXIS', 'OUI')

! - Get parameters of element
    call elref2(nomte, ntrou_max, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(2), fami=fami, ndim=ndim_fe, nno=nno2, &
                     npg=npg, jpoids=jv_poids, jvf=jv_vff2, jdfde=jv_dfde2)
    call elrefe_info(elrefe=lielrf(1), fami=fami, ndim=ndim_fe, nno=nno1, &
                     npg=npg, jpoids=jv_poids, jvf=jv_vff1, jdfde=jv_dfde1)
    ASSERT(ndim_sp .eq. ndim_fe+1)

! - Option parameters
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PDEPLMR', 'L', jv_deplm)
    call jevech('PCONTMR', 'L', jv_contm)
    call jevech('PVARIMR', 'L', jv_varim)
    call jevech('PDDEPLR', 'L', jv_ddepl)
    call jevech('PDEPL0R', 'L', jv_depl0)
    call jevech('PDEPL1R', 'L', jv_depl1)
    call jevech('PMATERC', 'L', jv_materc)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jv_carcri)
    call jevech('PCDTAU', 'L', jv_dtau)
    call jevech('PBORNPI', 'L', jv_borne)
    call jevech('PCOPILO', 'E', jv_copil)

! Number of internal variables
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jv_materc), materPara)

! - Local coordinate system
    call initLCSPg(ndim_sp, nno2, materPara)

! - Set main parameters for behaviour (on cell)
    instam = r8vide()
    instap = r8vide()
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jv_carcri), &
                              instam, instap, &
                              materPara, BEHInteg)

    ! Kinematics
    nddl = ndim_sp*(2*nno1+nno2)
    neps = 2*ndim_sp
    allocate (b(neps, npg, nddl), w(neps, npg), ni2ldc(neps, npg))
    call eimatb(nomte, ndim_sp, axi, nno1, nno2, npg, &
                zr(jv_poids), zr(jv_vff1), zr(jv_vff2), zr(jv_dfde2), zr(jv_geom), &
                materPara%lcsPara%lcsAnglePg, b, w, ni2ldc)

    ! Computation of path-following coefficients
    call ngpipe(BEHInteg, typmod, compor, ndim_sp, npg, neps, nddl, b, ni2ldc, &
                zr(jv_deplm), zr(jv_ddepl), zr(jv_depl0), zr(jv_depl1), &
                lgpg, zr(jv_contm), zr(jv_varim), &
                zr(jv_dtau), zr(jv_borne+1), zr(jv_borne), zr(jv_copil))

    deallocate (b, w, ni2ldc)

end subroutine
