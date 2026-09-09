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

subroutine te0543(option, nomte)
!
    use Behaviour_type
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/te0543_implement.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DES COEFFICIENTS A0 ET A1 POUR LE PILOTAGE PAR CRITERE ELASTIQUE
!  POUR LES ELEMENTS A VARIABLES LOCALES 2D et 3D
!
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    character(len=8) :: typmod(2)
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: ndim, nno, npg, lgpg, jtab(7)
    integer(kind=8) :: jv_poids, jv_vff, jv_dfde, jv_geom, jv_materc, jv_carcri
    integer(kind=8) :: jv_contm, jv_varim, jv_copil, jv_borne, jv_ctau, jv_typilo
    integer(kind=8) :: jv_deplm, jv_ddepl, jv_depl0, jv_depl1, iret
    real(kind=8) :: instam, instap
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
! --------------------------------------------------------------------------------------------------
!
! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), vattr_missing=' ')

! - Get parameters of element
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, jpoids=jv_poids, &
                     jvf=jv_vff, jdfde=jv_dfde)

! - Common parameters
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
    call jevech('PCDTAU', 'L', jv_ctau)
    call jevech('PBORNPI', 'L', jv_borne)
    call jevech('PCOPILO', 'E', jv_copil)

! Number of internal variables
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jv_materc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Set main parameters for behaviour (on cell)
    instam = r8vide()
    instap = r8vide()
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jv_carcri), &
                              instam, instap, &
                              materPara, BEHInteg)

    call behaviourPrepESVAGeom(nno, npg, ndim, &
                               jv_poids, jv_vff, jv_dfde, &
                               zr(jv_geom), BEHInteg)

! - Main subroutine to compute coefficients
    call te0543_implement(BEHInteg, &
                          typmod, compor, &
                          ndim, nno, npg, &
                          jv_poids, jv_vff, jv_dfde, zr(jv_geom), &
                          lgpg, zr(jv_deplm), zr(jv_contm), zr(jv_varim), &
                          zr(jv_ddepl), zr(jv_depl0), zr(jv_depl1), &
                          zr(jv_borne), zr(jv_ctau), zr(jv_copil))

end subroutine
