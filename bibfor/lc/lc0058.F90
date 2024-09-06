! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504,C1505
!
subroutine lc0058(BEHinteg, &
                  fami, kpg, ksp, ndim, typmod, &
                  imate, compor, carcri, instam, instap, &
                  neps, epsm, deps, nsig, sigm, &
                  nvi, vim, option, angmas, &
                  sigp, vip, ndsde, dsidep, codret)
!
    use Behaviour_type
    use logging_module, only: DEBUG, LOGLEVEL_MGIS, is_enabled
!
    implicit none
!
#include "asterc/mgis_debug.h"
#include "asterc/mgis_get_number_of_props.h"
#include "asterc/mgis_integrate.h"
#include "asterc/mgis_set_external_state_variables.h"
#include "asterc/mgis_set_gradients.h"
#include "asterc/mgis_set_internal_state_variables.h"
#include "asterc/mgis_set_material_properties.h"
#include "asterc/mgis_set_rotation_matrix.h"
#include "asterc/mgis_set_thermodynamic_forces.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcicma.h"
#include "asterfort/matrot.h"
#include "asterfort/mfront_get_mater_value.h"
#include "asterfort/mfrontPrepareStrain.h"
#include "asterfort/use_orient.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*), intent(in) :: fami
    integer, intent(in) :: kpg, ksp, ndim
    character(len=8), intent(in) :: typmod(*)
    integer, intent(in) :: imate
    character(len=16), intent(in) :: compor(*)
    real(kind=8), intent(in) :: carcri(*)
    real(kind=8), intent(in) :: instam, instap
    integer, intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    integer, intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer, intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: angmas(*)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    integer, intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                        merge(neps, 6, nsig*neps .eq. ndsde))
    integer, intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! MFRONT
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHinteg         : parameters for integration of behaviour
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  ndim             : dimension of problem (2 or 3)
! In  typmod           : type of modelization (TYPMOD2)
! In  imate            : coded material address
! In  compor           : name of comportment definition (field)
! In  carcri           : parameters for comportment
! In  instam           : time at beginning of time step
! In  instap           : time at end of time step
! In  neps             : number of components of strains
! In  epsm             : strains at beginning of current step time
! In  deps             : increment of strains during current step time
! In  nsig             : number of components of stresses
! In  sigm             : stresses at beginning of current step time
! In  nvi              : number of components of internal state variables
! In  vim              : internal state variables at beginning of current step time
! In  option           : name of option to compute
! In  angmas           : nautical angles
! Out sigm             : stresses at end of current step time
! Out vip              : internal state variables at end of current step time
! Out dsidep           : tangent matrix
! Out codret           : code for error
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lMatr, lSigm, lVari
    integer :: nstatv, i, j, nstran
    integer :: strain_model
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer, parameter :: s0 = 0, s1 = 1
    real(kind=8) :: drot(3, 3), dstran(9)
    real(kind=8) :: time(2)
    real(kind=8) :: ddsdde(54)
    real(kind=8) :: stran(9)
    real(kind=8) :: dtime, pnewdt
    character(len=16) :: rela_comp, defo_comp, extern_addr
    aster_logical :: l_greenlag, l_czm, l_pred
    real(kind=8) :: sigp_loc(6), vi_loc(nvi), dsidep_loc(6, 6)
    integer, parameter :: npropmax = 197
    real(kind=8) :: props(npropmax)
    integer :: ntens, ndi, nprops, retcode
    common/tdim/ntens, ndi
    aster_logical :: dbg
!
! --------------------------------------------------------------------------------------------------
!

    ASSERT(neps*nsig .eq. ndsde .or. (ndsde .eq. 36 .and. neps .le. 9 .and. nsig .le. 6))
    ASSERT(nsig .ge. 2*ndim)
    ASSERT(neps .ge. 2*ndim)

    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lMatr = L_MATR(option)

    sigp_loc = 0.d0
    vi_loc = 0.d0
    dsidep_loc = 0.d0
    stran = 0.d0
    dstran = 0.d0
    props = 0.d0

    dbg = is_enabled(LOGLEVEL_MGIS, DEBUG)

    ntens = 2*ndim
    ndi = 3
    codret = 0
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    l_pred = option(1:9) .eq. 'RIGI_MECA'
!
! - Finite element
!
    l_czm = typmod(2) .eq. 'ELEMJOIN'
    ASSERT(.not. l_czm)
!
! - Strain model
!
    strain_model = nint(carcri(EXTE_STRAIN))

    l_greenlag = defo_comp .eq. 'GREEN_LAGRANGE'

    if (ndim==2) then
        nstran = 5
    else if (ndim==3) then
        nstran = 9
    endif
!
! - Pointer to MGISBehaviour
!
    extern_addr = compor(MGIS_ADDR)
!
! - Get and set the material properties
    call mgis_get_number_of_props(extern_addr, nprops)
    ASSERT(nprops <= npropmax)
!
    call mfront_get_mater_value(extern_addr, BEHinteg, rela_comp, fami, kpg, &
                                ksp, imate, props, nprops)
!
! - Prepare strains
!
    call mfrontPrepareStrain(l_greenlag, l_pred, neps, epsm, deps, stran, dstran)
!
! - Number of internal state variables
!
    nstatv = nvi
!
! - Time parameters
!
    time(1) = instap-instam
    time(2) = instam
    dtime = instap-instam
!
! - Anisotropic case
!
    if (use_orient(angmas, 3)) then
        call matrot(angmas, drot)
        call mgis_set_rotation_matrix(extern_addr, drot)
    end if
!
! - Type of matrix for MFront
!
    ddsdde = 0.d0
    if (option .eq. 'RIGI_MECA_TANG') then
        ddsdde(1) = 4.d0
    else if (option .eq. 'RIGI_MECA_ELAS') then
        ddsdde(1) = 1.d0
    else if (option .eq. 'FULL_MECA_ELAS') then
        ddsdde(1) = 2.d0
    else if (option .eq. 'FULL_MECA') then
        ddsdde(1) = 4.d0
    else if (option .eq. 'RAPH_MECA') then
        ddsdde(1) = 0.d0
    end if
!
! - Call MFront
!
    pnewdt = 1.d0
    sigp_loc(1:nsig) = sigm(1:nsig)
    vi_loc(1:nstatv) = vim(1:nstatv)

    ! nstatv must be equal to the value returned by mgis_get_sizeof_isvs
    ! (not increased by kit...)

    if (dbg) then
        write (6, *) "+++ inputs +++ ", option
        write (6, *) "ddsdde", (ddsdde(i), i=1, ntens*ntens)
        write (6, *) "sigp_loc", (sigp_loc(i), i=1, 6)
        write (6, *) "vi_loc", (vi_loc(i), i=1, nstatv)
        write (6, *) "stran:", (stran(i), i=1, neps)
        write (6, *) "dstran:", (dstran(i), i=1, neps)
        write (6, *) "dtime:", dtime
        write (6, *) "predef:", (BEHinteg%exte%predef(i), i=1, BEHinteg%exte%nb_pred)
        write (6, *) "dpred:", (BEHinteg%exte%dpred(i), i=1, BEHinteg%exte%nb_pred)
        write (6, *) "props:", (props(i), i=1, nprops)
        write (6, *) "angl_naut:", (angmas(i), i=1, ndim)
        write (6, *) "ntens/nstatv:", ntens, nstatv
    end if

    call mgis_set_material_properties(extern_addr, s0, props, nprops)
    call mgis_set_gradients(extern_addr, s0, stran, nstran)
    call mgis_set_thermodynamic_forces(extern_addr, s0, sigp_loc, 2*ndim)
    call mgis_set_internal_state_variables(extern_addr, s0, vi_loc, nstatv)
    call mgis_set_external_state_variables(extern_addr, s0, BEHinteg%exte%predef, &
                                           BEHinteg%exte%nb_pred)

    call mgis_set_material_properties(extern_addr, s1, props, nprops)
    call mgis_set_gradients(extern_addr, s1, stran+dstran, nstran)
    call mgis_set_external_state_variables(extern_addr, s1, &
                                           BEHinteg%exte%predef+BEHinteg%exte%dpred, &
                                           BEHinteg%exte%nb_pred)

    ! call mgis_debug(extern_addr, "Before integration:")

    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA' .or. &
        option(1:9) .eq. 'RIGI_MECA') then
        call mgis_integrate(extern_addr, sigp_loc, vi_loc, ddsdde, dtime, &
                            pnewdt, retcode)
        ASSERT(nstatv .le. nvi)
    end if
!
    if (dbg) then
        write (6, *) "+++ outputs +++"
        write (6, *) "sigp_loc", (sigp_loc(i), i=1, 6)
        write (6, *) "vi_loc", (vi_loc(i), i=1, nstatv)
        write (6, *) "ddsdde", (ddsdde(i), i=1, ntens*ntens)
        write (6, *) "pnewdt/retcode:", pnewdt, retcode
        write (6, *) "ntens/nstatv:", ntens, nstatv
    end if
!
! - Convert stresses
!
!   TODO: sqrt(2) should be removed, same convention seems to be in used in MGIS/MFront
    ! sigp_loc(4:6) = sigp_loc(4:6)*rac2
!
! - Convert matrix
!
    if (option(1:9) .eq. 'RIGI_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        call lcicma(ddsdde, ntens, ntens, ntens, ntens, &
                    1, 1, dsidep_loc, 6, 6, 1, 1)
        dsidep_loc(1:6, 4:6) = dsidep_loc(1:6, 4:6)*rac2
        dsidep_loc(4:6, 1:6) = dsidep_loc(4:6, 1:6)*rac2
    end if
!
! - Returned code from mgis_integrate (retcode):
!    -1: integration failed
!     0: integration succeeded but results are unreliable
!     1: integration succeeded and results are reliable
!   Use 'pnewdt' to return state to caller:
    if (retcode .lt. 0) then
        codret = 1
    end if
    if (pnewdt .lt. 0.0d0) then
        if (pnewdt .lt. -0.99d0 .and. pnewdt .gt. -1.01d0) then
            codret = 1
        else if (pnewdt .lt. -1.99d0 .and. pnewdt .gt. -2.01d0) then
            call utmess('F', 'MFRONT_1')
        else if (pnewdt .lt. -2.99d0 .and. pnewdt .gt. -3.01d0) then
            call utmess('F', 'MFRONT_2')
        else if (pnewdt .lt. -3.99d0 .and. pnewdt .gt. -4.01d0) then
            codret = 1
        else
            call utmess('F', 'MFRONT_3')
        end if
    end if

    if (lSigm) then
        sigp = 0.d0
        sigp(1:2*ndim) = sigp_loc(1:2*ndim)
    end if
    if (lVari) then
        vip = 0.d0
        vip(1:nstatv) = vi_loc(1:nstatv)
    end if
    if (lMatr) then
        dsidep = 0.d0
        do i = 1, 2*ndim
            do j = 1, 2*ndim
                dsidep(i, j) = dsidep_loc(j, i)
            end do
        end do
    end if
!
end subroutine
