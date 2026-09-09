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
subroutine te0054(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nmsfin.h"
#include "asterfort/nmsfon.h"
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
! Elements: D_PLAN_MIX_STA et 3D_MIX_STA
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA, FORC_NODA et REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    character(len=8) :: typmod(2), elrefe
    aster_logical :: axi, lNonLine
    integer(kind=8) :: nno, npg, ndim, lgpg, nddl
    integer(kind=8) :: jv_poids, jv_vf, jv_dfde
    integer(kind=8) :: jvMaterc, icontm, ivarim, jvInstmr, jvInstpr, ideplm, ideplp
    integer(kind=8) :: ivectu, icontp, ivarip, imatuu, jvCarcri, ivarix, jvGeom, icoret
    integer(kind=8) :: icont
    integer(kind=8) :: iret, itab(7)
    integer(kind=8) :: codret
    real(kind=8) :: sigref, lagref
    real(kind=8), allocatable:: sref(:)
    aster_logical :: lMatr, lVect, lSigm, lVari, lrefe
    character(len=16), pointer :: compor(:) => null()
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!
    ivectu = 1
    icontp = 1
    ivarip = 1
    icoret = 1
    imatuu = 1

! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    typmod(2) = ' '
    axi = typmod(1) .eq. 'AXIS'
    lrefe = ASTER_FALSE

! - Get parameters of element
    call elref1(elrefe)
    call elrefe_info(elrefe=elrefe, fami=fami, &
                     ndim=ndim, nno=nno, &
                     npg=npg, jpoids=jv_poids, &
                     jvf=jv_vf, jdfde=jv_dfde)
    nddl = 3*nno*ndim

! - PARAMETRES EN ENTREE ET DIMENSION
    lNonLine = ASTER_FALSE
    call jevech('PGEOMER', 'L', jvGeom)
    if (option .eq. "FORC_NODA") then
        call jevech('PSIEFR', 'L', icont)
        call jevech('PVECTUR', 'E', ivectu)
    else if (option .eq. "REFE_FORC_NODA") then
        allocate (sref(4*ndim))
        call jevech('PVECTUR', 'E', ivectu)
    else
        lNonLine = ASTER_TRUE
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PINSTMR', 'L', jvInstmr)
        call jevech('PINSTPR', 'L', jvInstpr)
        call tecach('OOO', 'PDEPLPR', 'L', iret, nval=2, itab=itab)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=itab)
        lgpg = max(itab(6), 1)*itab(7)

    end if

    if (lNonLine) then
! ----- Material parameters
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), materPara)

! ----- Set local coordinate system from user
        call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! ----- Get fields for non-linear behaviour
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', jvCarcri)

! ----- Initialisation of behaviour datastructure
        call behaviourInit(BEHInteg)

! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

! ----- Access to other fields
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
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

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  zr(jvInstmr), zr(jvInstpr), &
                                  materPara, BEHInteg)

    end if

    if (option .eq. "FORC_NODA") then
        call nmsfon(lrefe, ndim, nno, npg, nddl, &
                    zr(jvGeom), zr(jv_vf), jv_dfde, &
                    jv_poids, zr(icont), zr(ivectu))

    else if (option .eq. "REFE_FORC_NODA") then
        lrefe = ASTER_TRUE
        call refe%Init(nomte)
        sigref = refe%GetRef('SIGM')
        lagref = refe%GetRef('LAGR')
        call refe%Check()

        sref(1:2*ndim) = sigref
        sref(2*ndim+1:4*ndim) = lagref
        call nmsfon(lrefe, ndim, nno, npg, nddl, &
                    zr(jvGeom), zr(jv_vf), jv_dfde, &
                    jv_poids, transpose(spread(sref, 1, npg)), zr(ivectu))

    else
        call nmsfin(BEHInteg, &
                    option, typmod, &
                    compor, zr(jvCarcri), &
                    ndim, nno, npg, nddl, &
                    jv_poids, zr(jv_vf), jv_dfde, &
                    zr(jvGeom), zr(ideplm), zr(ideplp), &
                    zr(jvInstmr), zr(jvInstpr), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    zr(ivectu), zr(imatuu), &
                    lMatr, lVect, lSigm, &
                    codret)

    end if

    if (lrefe) then
        deallocate (sref)
    end if

    if (lSigm) then
        zi(icoret) = codret
    end if

end subroutine
