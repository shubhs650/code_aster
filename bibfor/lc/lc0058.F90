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
subroutine lc0058(BEHInteg, &
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
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/lcicma.h"
#include "asterfort/matrot.h"
#include "asterfort/mfront_get_mater_value.h"
#include "asterfort/mfrontPrepareStrain.h"
#include "asterfort/use_orient.h"
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
! MFRONT
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
    real(kind=8) :: drot(3, 3), dstran(neps), stran(neps)
    real(kind=8) :: dsidepMGIS(merge(nsig, 6, nsig*neps .eq. ndsde)* &
                               merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8) :: dtime, pnewdt, rdt
    character(len=16) :: rela_comp, defo_comp, extern_addr
    aster_logical :: lGreenLagr, lCZM, lGradVari
    real(kind=8) :: sigp_loc(nsig), vi_loc(nvi)
    real(kind=8) :: dsidep_loc(merge(nsig, 6, nsig*neps .eq. ndsde), &
                               merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8) :: props(MGIS_MAX_PROPS)
    integer(kind=8) :: nprops, retcode
    aster_logical :: dbg
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

! - Flags for quantity to compute
    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lMatr = L_MATR(option)

! - Debug flag
    dbg = is_enabled(LOGLEVEL_MGIS, DEBUG)

! - Get main parameters
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    lCZM = typmod(2) .eq. 'ELEMJOIN' .or. typmod(2) .eq. 'INTERFAC'
    lGradVari = typmod(2) .eq. 'GRADVARI'
    lGreenLagr = defo_comp .eq. 'GREEN_LAGRANGE'
    ASSERT(.not. lCZM)
    ASSERT(.not. lGradVari)

! - Pointer to MGISBehaviour
    extern_addr = compor(MGIS_ADDR)

! - Management of dimensions
    call getMGISDime(lGreenLagr, lCZM, lGradVari, ndim, &
                     neps, nsig, nvi, ndsde, &
                     nstran, nforc, nstatv, nmatr)

! - Get and set the material properties
    call mgis_get_number_of_props(extern_addr, nprops)
    ASSERT(nprops <= MGIS_MAX_PROPS)
    call mfront_get_mater_value(extern_addr, BEHInteg, rela_comp, fami, kpg, &
                                ksp, jvMaterCode, props, nprops)

! - Prepare strains
    call mfrontPrepareStrain(lGreenLagr, neps, epsm, deps, stran, dstran)

! - Input stresses
    sigp_loc(1:nforc) = sigm(1:nforc)

! - Input internal state variables
    vi_loc(1:nstatv) = vim(1:nstatv)

! - Time increment
    dtime = instap-instam

! - Anisotropic case
    if (use_orient(materPara%lcsPara%lcsAngle, 3)) then
        call matrot(materPara%lcsPara%lcsAngle, drot)
        call mgis_set_rotation_matrix(extern_addr, drot)
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

! - Debug - Inputs
    if (dbg) then
        write (6, *) "+++ inputs +++ ", option
        write (6, *) "sigp_loc", (sigp_loc(i), i=1, nforc)
        write (6, *) "vi_loc", (vi_loc(i), i=1, nstatv)
        write (6, *) "stran:", (stran(i), i=1, nstran)
        write (6, *) "dstran:", (dstran(i), i=1, nstran)
        write (6, *) "dtime:", dtime
        write (6, *) "predef:", (BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(i), &
                                 i=1, BEHInteg%behavESVA%behavESVAExte%nbESVAScal)
        write (6, *) "dpred:", (BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(i), &
                                i=1, BEHInteg%behavESVA%behavESVAExte%nbESVAScal)
        write (6, *) "props:", (props(i), i=1, nprops)
        write (6, *) "angl_naut:", (materPara%lcsPara%lcsAngle(i), i=1, ndim)
        write (6, *) "nforc:", nforc
        write (6, *) "nstatv:", nstatv
        write (6, *) "nmatr:", nmatr
    end if

! - Set material properties (begin and end of current time step)
    call mgis_set_material_properties(extern_addr, s0, props, nprops)
    call mgis_set_material_properties(extern_addr, s1, props, nprops)

! - Set strains and increment of strains
    call mgis_set_gradients(extern_addr, s0, stran, nstran)
    call mgis_set_gradients(extern_addr, s1, stran+dstran, nstran)

! - Set stresses
    call mgis_set_thermodynamic_forces(extern_addr, s0, sigp_loc, nforc)

! - Set internal state variables
    call mgis_set_internal_state_variables(extern_addr, s0, vi_loc, nstatv)

! - Set external state variables (begin and end of current time step)
    call mgis_set_external_state_variables(extern_addr, s0, &
                                           BEHInteg%behavESVA%behavESVAExte%scalESVAPrev, &
                                           BEHInteg%behavESVA%behavESVAExte%nbESVAScal)
    call mgis_set_external_state_variables(extern_addr, s1, &
                                           BEHInteg%behavESVA%behavESVAExte%scalESVAPrev+ &
                                           BEHInteg%behavESVA%behavESVAExte%scalESVAIncr, &
                                           BEHInteg%behavESVA%behavESVAExte%nbESVAScal)

! - Désactivation de l'augmentation du pas de temps dans la LdC
    rdt = 1.d0

! - Call to integrator
    ! call mgis_debug(extern_addr, "Before integration:")
    call mgis_integrate(extern_addr, sigp_loc, vi_loc, dsidepMGIS, dtime, rdt, &
                        pnewdt, retcode)

! - Debug - Outputs
    if (dbg) then
        write (6, *) "+++ outputs +++"
        write (6, *) "sigp_loc", (sigp_loc(i), i=1, nforc)
        write (6, *) "vi_loc", (vi_loc(i), i=1, nstatv)
        write (6, *) "dsidepMGIS", (dsidepMGIS(i), i=1, nmatr*nmatr)
        write (6, *) "pnewdt( pas utilisé par aster)", pnewdt
        write (6, *) "retcode:", retcode
    end if

! - Convert stresses (nothing to do: same convention in code_aster)

! - Convert matrix
    if (lMatr) then
        call lcicma(dsidepMGIS, nmatr, nmatr, nmatr, nmatr, &
                    1, 1, dsidep_loc, 6, 6, &
                    1, 1)
    end if
!
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

! - Copy outputs
    if (lSigm) then
        sigp = 0.d0
        sigp(1:nforc) = sigp_loc(1:nforc)
    end if
    if (lVari) then
        vip = 0.d0
        vip(1:nstatv) = vi_loc(1:nstatv)
    end if
    if (lMatr) then
        dsidep = 0.d0
        do i = 1, nmatr
            do j = 1, nmatr
                dsidep(i, j) = dsidep_loc(j, i)
            end do
        end do
    end if
!
end subroutine
