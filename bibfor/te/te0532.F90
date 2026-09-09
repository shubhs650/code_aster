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
subroutine te0532(option, nomte)

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
#include "asterfort/nmbamb.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ngpipe.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
! Elementary computation PILO_PRED_ELAS
! Elements: BARRE_2D_NL, BARRE_3D_NL
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    character(len=8) :: typmod(2), attrib
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: ndim_sp, nno, npg, lgpg, jtab(7), nddl, neps
    integer(kind=8) :: jv_poids, jv_vff, jv_dfde
    integer(kind=8) :: jv_geom, jv_materc, jv_carcri, jv_sect
    integer(kind=8) :: jv_contm, jv_varim, jv_copil, jv_borne, jv_dtau, jv_typilo
    integer(kind=8) :: jv_deplm, jv_ddepl, jv_depl0, jv_depl1, iret
    real(kind=8) :: instam, instap, aire
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
    real(kind=8), allocatable:: b(:, :, :), w(:, :), ni2ldc(:, :)
! --------------------------------------------------------------------------------------------------

! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), vattr_missing=' ')
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp

! - Get parameters of element
    call elrefe_info(fami=fami, nno=nno, npg=npg, jpoids=jv_poids, jvf=jv_vff, jdfde=jv_dfde)

! - Option parameters
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PCAGNBA', 'L', jv_sect)
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

    ! Cross section area
    aire = zr(jv_sect)

! Number of internal variables
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jv_materc), materPara)

    ! No definition of local coordinate system
    call initLCSNone(materPara)

! - Set main parameters for behaviour (on cell)
    instam = r8vide()
    instap = r8vide()
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jv_carcri), &
                              instam, instap, &
                              materPara, BEHInteg)

    ! Kinematics
    call nmbamb(ndim_sp, nno, npg, zr(jv_geom), aire, &
                zr(jv_dfde), zr(jv_poids), nddl, neps, b, w, ni2ldc)

    ! Computation of path-following coefficients
    call ngpipe(BEHInteg, typmod, compor, ndim_sp, npg, neps, nddl, b, ni2ldc, &
                zr(jv_deplm), zr(jv_ddepl), zr(jv_depl0), zr(jv_depl1), &
                lgpg, zr(jv_contm), zr(jv_varim), &
                zr(jv_dtau), zr(jv_borne+1), zr(jv_borne), zr(jv_copil))

    deallocate (b, w, ni2ldc)

end subroutine
