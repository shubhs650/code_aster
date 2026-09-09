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
! aslint: disable=W1306,W1504,C1505,W0104
!
subroutine lc9058(BEHInteg, &
                  fami, kpg, ksp, ndim, &
                  typmod, jvMaterCode, compor, carcri, instam, &
                  instap, neps, epsm, deps, nsig, &
                  sigm, nvi, vim, option, &
                  sigp, vip, ndsde, dsidep, codret)
!
    use Behaviour_type
    use BehaviourMGIS_module
    use logging_module, only: DEBUG, LOGLEVEL_MGIS, is_enabled
    use MaterialPara_type
    implicit none
!
#include "asterc/mgis_get_number_of_props.h"
#include "asterc/mgis_integrate.h"
#include "asterc/mgis_set_gradients.h"
#include "asterc/mgis_set_internal_state_variables.h"
#include "asterc/mgis_set_material_properties.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/czm_post.h"
#include "asterfort/mfront_get_mater_value.h"
#include "asterfort/mfrontPrepareStrain.h"
#include "asterfort/rcvalb.h"
#include "asterfort/use_orient.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, ndim
    character(len=8), intent(in) :: typmod(2)
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                        merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! MFRONT for CZM with INTERFAC elements (CZM_MFRONT)
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHInteg         : parameters for integration of behaviour
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  ndim             : dimension of problem (2 or 3)
! In  typmod           : type of modelization (TYPMOD2)
! In  jvMaterCode            : coded material address
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
! Out sigp             : stresses at end of current step time
! Out vip              : internal state variables at end of current step time
! Out dsidep           : tangent matrix
! Out codret           : code for error
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lMatr, lSigm, lVari
    integer(kind=8) :: i, j
    integer(kind=8) :: nstran, nforc, nstatv, nmatr
    integer(kind=8), parameter :: s0 = 0, s1 = 1
    real(kind=8) :: dstran(ndim), stran(ndim)
    real(kind=8) :: dsidepMGIS(merge(nsig, 6, nsig*neps .eq. ndsde)* &
                               merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8) :: dtime, pnewdt, rdt
    character(len=16) :: relaComp, defoComp, extern_addr
    aster_logical :: lGreenLagr, lCZM, lGradVari
    real(kind=8) :: sigp_loc(ndim), vi_loc(nvi)
    real(kind=8) :: dsidep_loc(ndim, ndim)
    real(kind=8) :: props(MGIS_MAX_PROPS)
    integer(kind=8) :: nprops, retcode
    aster_logical :: dbg
    integer(kind=8) :: cod(1)
    real(kind=8) :: val(1)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = BEHInteg%materPara
    sigp_loc = 0.d0
    vi_loc = 0.d0
    dsidep_loc = 0.d0
    stran = 0.d0
    dstran = 0.d0
    props = 0.d0
    codret = 0

    ASSERT(nsig .ge. 2*ndim)
    ASSERT(neps .ge. 2*ndim)

! - Flags for quantity to compute
    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lMatr = L_MATR(option)

! - Debug flag
    dbg = is_enabled(LOGLEVEL_MGIS, DEBUG)

! - Get main parameters
    relaComp = compor(RELA_NAME)
    defoComp = compor(DEFO)
    lCZM = typmod(2) .eq. 'INTERFAC'
    lGradVari = typmod(2) .eq. 'GRADVARI'
    lGreenLagr = defoComp .eq. 'GREEN_LAGRANGE'
    ASSERT(lCZM)

! - Pointer to MGISBehaviour
    extern_addr = compor(MGIS_ADDR)

! - Management of dimensions
    call getMGISDime(lGreenLagr, lCZM, lGradVari, ndim, &
                     neps, nsig, nvi, ndsde, &
                     nstran, nforc, nstatv, nmatr)

! - Get and set the material properties
    call mgis_get_number_of_props(extern_addr, nprops)
    ASSERT(nprops <= MGIS_MAX_PROPS)
    call mfront_get_mater_value(extern_addr, BEHInteg, relaComp, fami, kpg, &
                                ksp, jvMaterCode, props, nprops)

! - Prepare strains
    call rcvalb(fami, kpg, ksp, '+', jvMaterCode, ' ', 'MFRONT', 0, ' ', [0.d0], &
                1, 'PENA_LAGR', val, cod, 2)
    call mfrontPrepareStrain(lGreenLagr, ndim, &
                             epsm(ndim+1:2*ndim)+val(1)*epsm(1:ndim), &
                             deps(ndim+1:2*ndim)+val(1)*deps(1:ndim), &
                             stran, dstran)

! - Input internal state variables
    vi_loc(1:nstatv) = vim(1:nstatv)

! - Time increment
    dtime = instap-instam

! - Anisotropic case
    if (use_orient(materPara%lcsPara%lcsAngle, 3)) then
        call utmess('F', 'MGIS1_2', sk=typmod(2))
    end if

! - Type of matrix for MFront
    dsidepMGIS = 0.d0
    if (option .eq. 'RIGI_MECA_TANG') then
        dsidepMGIS(1) = float(MGIS_BV_INTEGRATION_CONSISTENT_TANGENT_OPERATOR)
    else if (option .eq. 'RIGI_MECA_ELAS') then
        dsidepMGIS(1) = float(MGIS_BV_INTEGRATION_ELASTIC_OPERATOR)
    else if (option .eq. 'FULL_MECA_ELAS') then
        dsidepMGIS(1) = float(MGIS_BV_INTEGRATION_SECANT_OPERATOR)
    else if (option .eq. 'FULL_MECA') then
        dsidepMGIS(1) = float(MGIS_BV_INTEGRATION_CONSISTENT_TANGENT_OPERATOR)
    else if (option .eq. 'RAPH_MECA') then
        dsidepMGIS(1) = float(MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR)
    else
        WRITE (6, *) "Option <", option, ">"
    end if

! - Set material properties (begin and end of current time step)
    call mgis_set_material_properties(extern_addr, s0, props, nprops)
    call mgis_set_material_properties(extern_addr, s1, props, nprops)

! - Set strains and increment of strains
    call mgis_set_gradients(extern_addr, s0, stran, nstran)
    call mgis_set_gradients(extern_addr, s1, stran+dstran, nstran)

! - Set internal state variables
    call mgis_set_internal_state_variables(extern_addr, s0, vi_loc, nstatv)

! - Désactivation de l'augmentation du pas de temps dans la LdC
    rdt = 1.d0

! - Call to integrator
    call mgis_integrate(extern_addr, sigp_loc, vi_loc, dsidepMGIS, dtime, rdt, &
                        pnewdt, retcode)

! - Convert matrix
    if (lMatr) then
        do i = 1, ndim
            do j = 1, ndim
                dsidep_loc(i, j) = dsidepMGIS(ndim*(i-1)+j)
            end do
        end do
    end if

! - Returned code from mgis_integrate (retcode):
!    -1: integration failed
!     0: integration succeeded but results are unreliable
!     1: integration succeeded and results are reliable
    if (retcode .eq. -1) then
        codret = 1
    elseif (retcode .eq. 0) then
        codret = 2
    elseif (retcode .eq. 1) then
        codret = 0
    else
        ASSERT(ASTER_FALSE)
    end if

! - Copy output
    call czm_post(ndim, lSigm, lMatr, val(1), &
                  epsm(ndim+1:2*ndim)+deps(ndim+1:2*ndim), &
                  epsm(1:ndim)+deps(1:ndim), &
                  sigp_loc, dsidep_loc, sigp, dsidep)

    if (lVari) then
        vip = 0.d0
        vip(1:nstatv) = vi_loc(1:nstatv)
    end if

end subroutine
