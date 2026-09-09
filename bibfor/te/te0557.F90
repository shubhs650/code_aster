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

subroutine te0557(option, nomte)
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
#include "asterfort/jevech.h"
#include "asterfort/ngfint.h"
#include "asterfort/nmbamb.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"

    character(len=16), intent(in) :: option, nomte

! --------------------------------------------------------------------------------------------------
! Elements: BARRE / 2D_BARRE
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
    aster_logical, parameter:: matsym = ASTER_TRUE
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: ndim_fe, ndim_sp, nno, npg, nddl, neps, iret, lgpg
    integer(kind=8):: jv_poids, jv_dfde
    integer(kind=8) :: jv_materc, jv_geom, jv_sect, jv_instmr
    integer(kind=8) :: jv_instpr, jv_deplm, jv_deplp, jv_contm, jv_varim
    integer(kind=8) :: jv_carcri, jv_matuu, jv_matns
    integer(kind=8) :: jv_vectu, jv_contp, jv_varip, jv_codret
    integer(kind=8) :: codret, itab(7)
    real(kind=8) :: aire
    character(len=8) :: typmod(2), attrib
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: defoComp
    aster_logical :: lVect, lMatr, lVari, lSigm
    real(kind=8), allocatable :: b(:, :, :), w(:, :), ni2ldc(:, :)
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
! --------------------------------------------------------------------------------------------------

! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), vattr_missing=' ')

! - Get parameters of element
    call elrefe_info(fami=fami, ndim=ndim_fe, nno=nno, npg=npg, jpoids=jv_poids, jdfde=jv_dfde)
    call teattr('S', 'DIM_COOR_MODELI', attrib)
    read (attrib, '(I8)') ndim_sp

! - Get input fields
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PCAGNBA', 'L', jv_sect)
    call jevech('PINSTMR', 'L', jv_instmr)
    call jevech('PINSTPR', 'L', jv_instpr)
    call jevech('PDEPLMR', 'L', jv_deplm)
    call jevech('PDEPLPR', 'L', jv_deplp)
    call jevech('PCONTMR', 'L', jv_contm)
    call jevech('PVARIMR', 'L', jv_varim)
    call jevech('PMATERC', 'L', jv_materc)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jv_carcri)

    ! Cross section area
    aire = zr(jv_sect)

    ! Initializations of material parameters on current cell
    call initParaCell(fami, zi(jv_materc), materPara)

    ! No definition of local coordinate system
    call initLCSNone(materPara)

    ! Properties of behaviour
    defoComp = compor(DEFO)
    ASSERT(defoComp .eq. 'PETIT')

    ! Number of internal variables
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=itab)
    lgpg = max(itab(6), 1)*itab(7)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jv_carcri), &
                              zr(jv_instmr), zr(jv_instpr), &
                              materPara, BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Get output fields
    if (lMatr) then
        call jevech('PMATUUR', 'E', jv_matuu)
        jv_matns = jv_matuu
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', jv_vectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', jv_contp)
        call jevech('PCODRET', 'E', jv_codret)
    end if
    ! if (lVari .or. option(1:10) .eq. 'RIGI_MECA_') then
    if (lVari) then
        call jevech('PVARIPR', 'E', jv_varip)
    end if

    ! if ((relaComp .eq. 'ELAS') .or. (relaComp .eq. 'VMIS_ISOT_LINE') .or. &
    !     (relaComp .eq. 'VMIS_ISOT_TRAC') .or. (relaComp .eq. 'CORR_ACIER') .or. &
    !     (relaComp .eq. 'VMIS_CINE_LINE') .or. (relaComp .eq. 'RELAX_ACIER')) then

    ! Calcul de la cinématique
    call nmbamb(ndim_sp, nno, npg, zr(jv_geom), aire, &
                zr(jv_dfde), zr(jv_poids), nddl, neps, b, w, ni2ldc)

    ! Intégration du comportement et calcul des vecteurs et matrices élémentaires
    call ngfint(BEHInteg, &
                option, typmod, ndim_fe, nddl, neps, &
                npg, w, b, compor, &
                lgpg, zr(jv_carcri), zr(jv_instmr), &
                zr(jv_instpr), zr(jv_deplm), zr(jv_deplp), ni2ldc, zr(jv_contm), &
                zr(jv_varim), zr(jv_contp), zr(jv_varip), zr(jv_vectu), matsym, &
                zr(jv_matuu), zr(jv_matns), lMatr, lVect, lSigm, &
                codret)

    deallocate (b, w, ni2ldc)

    if (lSigm) zi(jv_codret) = codret

end subroutine
