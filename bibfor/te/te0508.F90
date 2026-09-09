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
subroutine te0508(option, nomte)
!
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefv.h"
#include "asterfort/jevech.h"
#include "asterfort/lggvfc.h"
#include "asterfort/lgicfc.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngforc.h"
#include "asterfort/nmgvmb.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
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
! Options: FORC_NODA, REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: defo_comp
    character(len=8) :: typmod(2)
    aster_logical :: axi, grand, inco, lrefe
    integer(kind=8) :: nnoQ, nnoL, npg, ndim, nddl, neps, itab(2)
    integer(kind=8) :: iret, nnos, jv_ganoQ, jv_poids, jv_vfQ, jv_dfdeQ, jv_vfL, jv_dfdeL, jv_ganoL
    integer(kind=8) :: igeom, icont, ivectu, idepl
    character(len=16), pointer :: compor(:) => null()
    real(kind=8) :: sigref, varref, lagref, epsref
    real(kind=8), allocatable:: b(:, :, :), w(:, :), ni2ldc(:, :)
    real(kind=8), allocatable:: sref(:)
    real(kind=8), allocatable:: ddl(:)
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!
    lrefe = option .eq. 'REFE_FORC_NODA'
!
! - Type of modelling
!
    call teattr('S', 'TYPMOD', typmod(1))
    typmod(2) = ' '
    inco = lteatt('INCO', 'C5GV')
    axi = typmod(1) .eq. 'AXIS'
!
! - Get parameters of element
!
    call elrefv('RIGI', ndim, &
                nnoL, nnoQ, nnos, &
                npg, jv_poids, &
                jv_vfL, jv_vfQ, &
                jv_dfdeL, jv_dfdeQ, &
                jv_ganoL, jv_ganoQ)
    neps = merge(3*ndim+4, 3*ndim+2, inco)

! Parametres de l'option et nbr de ddl
    call jevech('PGEOMER', 'L', igeom)

    call tecach('OOO', 'PVECTUR', 'E', iret, nval=2, itab=itab)
    ivectu = itab(1)
    nddl = itab(2)
    if (.not. lrefe) then
        call jevech('PSIEFR', 'L', icont)
        call jevech('PDEPLAR', 'L', idepl)
    else
        allocate (sref(neps))
        ! En attendant de lire le deplacement dans l'option REFE_FORC_NODA
        allocate (ddl(nddl))
        ddl = 0
    end if

!
! - Behaviour
!
    call jevech('PCOMPOR', 'L', vk16=compor)
    defo_comp = compor(DEFO)
    grand = defo_comp(1:8) .eq. 'GDEF_LOG'

! -------------------------!
!   GRAD_INCO + GDEF_LOG   !
! -------------------------!

    if (inco .and. grand) then

        if (lrefe) then
            call refe%Init(nomte)
            sigref = refe%GetRef('SIGM')
            varref = refe%GetRef('VARI')
            lagref = refe%GetRef('LAGR')
            epsref = refe%GetRef('EPSI')
            call refe%Check()

            if (ndim .eq. 2) then
                sref(1:neps) = [sigref, sigref, sigref, sigref, epsref, &
                                sigref, lagref, varref, 0.d0, 0.d0]
            else if (ndim .eq. 3) then
                sref(1:neps) = [sigref, sigref, sigref, sigref, sigref, &
                                sigref, epsref, sigref, lagref, varref, &
                                0.d0, 0.d0, 0.d0]
            end if

            call lgicfc(lrefe, ndim, nnoQ, nnoL, npg, nddl, axi, &
                        zr(igeom), ddl, zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                        jv_poids, transpose(spread(sref, 1, npg)), &
                        zr(ivectu))

        else
            call lgicfc(lrefe, ndim, nnoQ, nnoL, npg, nddl, axi, &
                        zr(igeom), zr(idepl), zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                        jv_poids, zr(icont), zr(ivectu))

        end if

! -------------------------!
!   GRAD_VARI + GDEF_LOG   !
! -------------------------!

    else if (.not. inco .and. grand) then

        if (lrefe) then
            call refe%Init(nomte)
            sigref = refe%GetRef('SIGM')
            varref = refe%GetRef('VARI')
            lagref = refe%GetRef('LAGR')
            call refe%Check()

            if (ndim .eq. 2) then
                sref(1:neps) = [sigref, sigref, sigref, sigref, lagref, &
                                varref, 0.d0, 0.d0]
            else if (ndim .eq. 3) then
                sref(1:neps) = [sigref, sigref, sigref, sigref, sigref, &
                                sigref, lagref, varref, 0.d0, 0.d0, 0.d0]
            end if

            call lggvfc(lrefe, ndim, nnoQ, nnoL, npg, nddl, axi, &
                        zr(igeom), ddl, zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                        jv_poids, transpose(spread(sref, 1, npg)), &
                        zr(ivectu))

        else
            call lggvfc(lrefe, ndim, nnoQ, nnoL, npg, nddl, axi, &
                        zr(igeom), zr(idepl), zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                        jv_poids, zr(icont), zr(ivectu))

        end if

! -------------------------!
!   GRAD_VARI + PETIT      !
! -------------------------!

    else if (.not. inco .and. .not. grand) then
        call nmgvmb(ndim, nnoQ, nnoL, npg, axi, &
                    zr(igeom), zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                    jv_poids, nddl, neps, b, w, ni2ldc)

        if (lrefe) then
            call refe%Init(nomte)
            sigref = refe%GetRef('SIGM')
            varref = refe%GetRef('VARI')
            lagref = refe%GetRef('LAGR')
            call refe%Check()

            if (ndim .eq. 2) then
                sref(1:neps) = [sigref, sigref, sigref, sigref, lagref, &
                                varref, 0.d0, 0.d0]
            else if (ndim .eq. 3) then
                sref(1:neps) = [sigref, sigref, sigref, sigref, sigref, &
                                sigref, lagref, varref, 0.d0, 0.d0, 0.d0]
            end if

            call ngforc(w, abs(b), ni2ldc, transpose(spread(sref(1:neps), 1, npg)), zr(ivectu))

        else
            call ngforc(w, b, ni2ldc, zr(icont), zr(ivectu))

        end if

        deallocate (b, w, ni2ldc)

    else
        call utmess('F', 'ELEMENTS_2')
        ! Combinaison inconnue
    end if

    if (lrefe) then
        deallocate (ddl)
        deallocate (sref)
    end if

end subroutine
