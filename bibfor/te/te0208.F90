! --------------------------------------------------------------------
! Copyright (C) 2007 - 2026 - EDF - www.code-aster.org
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

subroutine te0208(option, nomte)

    use Behaviour_type
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none

#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/te0208_implement.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16) :: nomte, option
!-----------------------------------------------------------------------
!     PILOTAGE POUR LES ELEMENTS DE JOINT 3D
!     OPTION : PILO_PRED_ELAS
!-----------------------------------------------------------------------
!
    character(len=4), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    character(len=8) :: typmod(2), attrib
    integer(kind=8) :: nddl, npg, jtab(7), iret, nno, ndim_sp, lgpg
    integer(kind=8) :: iddepl, idepl0, idepl1, ictau, icopil, jv_carcri, jv_bornes, jv_contm
    integer(kind=8) :: jv_geom, jv_materc, ideplm, ivarim, jv_poids, jv_vff, jv_dfde
    character(len=16), pointer :: compor(:) => null()
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
! --------------------------------------------------------------------------------------------------

    ! Finite element characteristics (only the boundary face, not the volume)
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))
    call elrefe_info(fami=fami, npg=npg, nno=nno, jpoids=jv_poids, jvf=jv_vff, jdfde=jv_dfde)
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp
    nddl = ndim_sp*(2*nno)
!
! - PARAMETRES EN ENTREE
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_materc)
    call jevech('PCARCRI', 'L', jv_carcri)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PCONTMR', 'L', jv_contm)
    call jevech('PBORNPI', 'L', jv_bornes)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDDEPLR', 'L', iddepl)
    call jevech('PDEPL0R', 'L', idepl0)
    call jevech('PDEPL1R', 'L', idepl1)
    call jevech('PCDTAU', 'L', ictau)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCOPILO', 'E', icopil)

    ! Initializations of material parameters on current cell
    call initParaCell(fami, zi(jv_materc), materPara)
    call initLCSNone(materPara)
    call behaviourSetParaCell(typmod, option, compor, zr(jv_carcri), r8vide(), r8vide(), &
                              materPara, BEHInteg)

    ! recuperation du nombre de variables internes par points de gauss :
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

    ! Compute path-following coefficients
    call te0208_implement(BEHInteg, typmod, compor, ndim_sp, nno, nddl, npg, zr(jv_geom), &
                          zr(jv_poids), zr(jv_vff), zr(jv_dfde), &
                          zr(ideplm), zr(iddepl), zr(idepl0), zr(idepl1), &
                          lgpg, zr(jv_contm), zr(ivarim), zr(jv_bornes+1), zr(jv_bornes), &
                          zr(ictau), zr(icopil))

end subroutine
