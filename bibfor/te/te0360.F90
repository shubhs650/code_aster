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
subroutine te0360(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/eimatb.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngfint.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_INTERFACE
!           PLAN_INTERFACE, AXIS_INTERFACE
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
    character(len=16) :: defoComp
    character(len=8) :: typmod(2), lielrf(10)
    aster_logical :: axi, matsym, lVect, lMatr, lVari, lSigm
    integer(kind=8) :: iret, ntrou, codret, jtab(7)
    integer(kind=8) :: ndim, nno1, nno2, npg, lgpg, nddl, neps
    integer(kind=8) :: jv_w, jv_vff1, jvGeom, jvMaterc, jv_vff2, jv_dff2
    integer(kind=8) :: jv_varim, jv_varip, jvInstmr, jvInstpr, jv_codret
    integer(kind=8) :: jv_ddlm, jv_ddld, jvCarcri
    integer(kind=8) :: jv_vectu, jv_contp, jv_varix, jv_matuu, jv_matns
    character(len=16), pointer :: compor(:) => null()
    real(kind=8), allocatable:: wg(:, :), ni2ldc(:, :), b(:, :, :)
    real(kind=8), allocatable:: sigm(:, :), transp_ktan(:, :)
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    jv_vectu = 1
    jv_contp = 1
    jv_varip = 1
    jv_codret = 1
    jv_matuu = 1
    jv_matns = 1
    codret = 0

! - Get element parameters
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami=fami, ndim=ndim, nno=nno1, npg=npg, &
                     jpoids=jv_w, jvf=jv_vff1)
    call elrefe_info(elrefe=lielrf(2), fami=fami, ndim=ndim, nno=nno2, npg=npg, &
                     jpoids=jv_w, jvf=jv_vff2, jdfde=jv_dff2)
    ndim = ndim+1
    nddl = ndim*(2*nno1+nno2)
    neps = 2*ndim

    allocate (b(neps, npg, nddl), wg(neps, npg), ni2ldc(neps, npg))
    allocate (sigm(neps, npg), transp_ktan(nddl, nddl))

! - Type of finite element
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))
    axi = lteatt('AXIS', 'OUI')

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PVARIMR', 'L', jv_varim)
    call jevech('PDEPLMR', 'L', jv_ddlm)
    call jevech('PDEPLPR', 'L', jv_ddld)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)
    call initParaCell(fami, zi(jvMaterc), materPara)
    call initLCSPg(ndim, nno2, materPara)

! - Get fields for non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Properties of behaviour
    defoComp = compor(DEFO)
    if (defoComp(1:5) .ne. 'PETIT') then
        call utmess('F', 'JOINT1_2', sk=defoComp)
    end if

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Total number of internal state variables on element
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Get time
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)

! - Get output fields
    if (lVect) call jevech('PVECTUR', 'E', jv_vectu)
    if (lSigm) call jevech('PCONTPR', 'E', jv_contp)
    if (lVari) then
        call jevech('PVARIPR', 'E', jv_varip)
        call jevech('PVARIMP', 'L', jv_varix)
        zr(jv_varip:jv_varip+npg*lgpg-1) = zr(jv_varix:jv_varix+npg*lgpg-1)
    end if
    if (lMatr) then
        call jevech('PMATUUR', 'E', jv_matuu)
        call jevech('PMATUNS', 'E', jv_matns)
        matsym = .not. (nint(zr(jvCarcri-1+CARCRI_MATRSYME)) .gt. 0)
    end if

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              zr(jvInstmr), zr(jvInstpr), &
                              materPara, BEHInteg)

! - Calcul de la matrice cinematique
    call eimatb(nomte, ndim, axi, nno1, nno2, npg, &
                zr(jv_w), zr(jv_vff1), zr(jv_vff2), zr(jv_dff2), zr(jvGeom), &
                materPara%lcsPara%lcsAnglePg, b, wg, ni2ldc)

! - Calcul des forces interieures et de la matrice tangente transposee
    sigm = 0
    call ngfint(BEHInteg, &
                option, typmod, ndim, nddl, neps, &
                npg, wg, b, compor, &
                lgpg, zr(jvCarcri), zr(jvInstmr), &
                zr(jvInstpr), zr(jv_ddlm), zr(jv_ddld), ni2ldc, sigm, &
                zr(jv_varim), zr(jv_contp), zr(jv_varip), zr(jv_vectu), &
                matsym, zr(jv_matuu), zr(jv_matns), lMatr, lVect, lSigm, codret)

! - Stockage du code retour
    if (lSigm) then
        call jevech('PCODRET', 'E', jv_codret)
        zi(jv_codret) = codret
    end if

    deallocate (b, wg, ni2ldc, sigm, transp_ktan)
end subroutine
