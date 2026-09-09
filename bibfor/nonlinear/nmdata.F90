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
subroutine nmdata(model, mesh, mater, mateco, cara_elem, ds_constitutive, &
                  listLoad, solver, ds_conv, sddyna, ds_posttimestep, &
                  ds_energy, ds_errorindic, ds_print, ds_algopara, &
                  ds_inout, ds_contact, ds_measure, ds_algorom, &
                  nlDynaDamping)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use Rom_Datastructure_type
    use NonLinearDyna_module
    use listLoad_type
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/GetIOField.h"
#include "asterfort/gettco.h"
#include "asterfort/infdbg.h"
#include "asterfort/ndcrdy.h"
#include "asterfort/ndlect.h"
#include "asterfort/nmdoch.h"
#include "asterfort/nmdocn.h"
#include "asterfort/nmdomt_ls.h"
#include "asterfort/nmdomt.h"
#include "asterfort/nmdopo.h"
#include "asterfort/nmdorc.h"
#include "asterfort/nmlect.h"
#include "asterfort/nonlinDSContactRead.h"
#include "asterfort/nonlinDSEnergyRead.h"
#include "asterfort/nonlinDSErrorIndicRead.h"
#include "asterfort/nonlinDSInOutRead.h"
#include "asterfort/nonlinDSMeasureRead.h"
#include "asterfort/nonlinDSPrintRead.h"
#include "asterfort/nonlinDSPrintSepLine.h"
#include "asterfort/utmess.h"
#include "asterfort/verif_affe.h"
!
    character(len=24), intent(out) :: model
    character(len=*), intent(out) :: mesh, mater, mateco, cara_elem
    type(NL_DS_Constitutive), intent(inout) :: ds_constitutive
    character(len=*), intent(out) :: listLoad, solver
    type(NL_DS_Conv), intent(inout) :: ds_conv
    character(len=19) :: sddyna
    type(NL_DS_PostTimeStep), intent(inout) :: ds_posttimestep
    type(NL_DS_Energy), intent(inout) :: ds_energy
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    type(NL_DS_InOut), intent(inout) :: ds_inout
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
    type(NL_DS_ErrorIndic), intent(inout) :: ds_errorindic
    type(NLDYNA_DAMPING), intent(out) :: nlDynaDamping
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Read parameters
!
! --------------------------------------------------------------------------------------------------
!
! Out mesh             : name of mesh
! Out model            : name of model
! Out mater            : name of material characteristics (field)
! Out mateco           : name of coded material
! Out cara_elem        : name of elementary characteristics (field)
! IO  ds_constitutive  : datastructure for constitutive laws management
! Out listLoad        : name of datastructure for list of loads
! Out solver           : name of datastructure for solver
! IO  ds_conv          : datastructure for convergence management
! IN  SDDYNA : SD DYNAMIQUE
! IO  ds_posttimestep  : datastructure for post-treatment at each time step
! IO  ds_energy        : datastructure for energy management
! IO  ds_errorindic    : datastructure for error indicator
! IO  ds_print         : datastructure for printing parameters
! IO  ds_algopara      : datastructure for algorithm parameters
! IO  ds_inout         : datastructure for input/output management
! IO  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
! IO  ds_algorom       : datastructure for ROM parameters
! Out nlDynaDamping    : damping parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=1), parameter:: jvBase = "V"
    character(len=8) :: result, modelCara
    character(len=16) :: k16bid, nomcmd
    aster_logical :: l_etat_init, l_sigm, lHasPilo, staticOperator
    character(len=24) :: typco
    character(len=8) :: stin_evol, cara_elem_in
    type(ListLoad_Prep) :: listLoadPrep
    integer(kind=8) :: nocc
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call nonlinDSPrintSepLine()
        call utmess('I', 'MECANONLINE12_1')
    end if

    model = " "

! - Get command parameters
    call getres(result, k16bid, nomcmd)
    staticOperator = nomcmd .eq. 'STAT_NON_LINE'

! - Read parameters for input/output management
    call nonlinDSInOutRead('MECA', result, ds_inout)

! - Initial state (complete datastructure or only stresses)
    call GetIOField(ds_inout, 'SIEF_ELGA', l_read_=l_sigm)
    l_etat_init = ((ds_inout%l_stin_evol) .or. (l_sigm))

! - LECTURE DONNEES GENERALES
    call nmlect(result, model, mater, mateco, cara_elem, solver)

! - Continueation method
    lHasPilo = ASTER_FALSE
    if (staticOperator) then
        call getfac('PILOTAGE', nocc)
        lHasPilo = nocc .gt. 0
    end if

! - Get loads/BC and create list of loads datastructure
    listLoad = '&&OP00XX.LIST_LOAD'
    listLoadPrep%model = model(1:8)
    listLoadPrep%lHasPilo = lHasPilo
    listLoadPrep%funcIsCplx = ASTER_FALSE
    listLoadPrep%staticOperator = staticOperator
    call nmdoch(listLoadPrep, listLoad, jvBase)

! - Get mesh (only one !)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
!   ------------------------------------------------------------------------------------------------
!   VERIFICATION OF CARA_ELEM
!
!   COEF_RIGI_DRZ prohibited in non-linear
    if (cara_elem .ne. ' ') then
        call gettco(cara_elem, typco)
        if (typco .eq. 'CARA_ELEM') then
            call verif_affe(model, cara_elem, non_lin=ASTER_TRUE)
        end if
    end if

! - Check consistency of models for CARa_ELEM
    if (cara_elem .ne. ' ') then
        call dismoi("NOM_MODELE", cara_elem, "CARA_ELEM", repk=modelCara)
        if (modelCara .ne. model) then
            call utmess("A", "MECANONLINE5_74")
        end if
    end if

!   IF exist ETAT_INIT/EVOL_NOLI : Is cara_elem in ETAT_INIT same as curent ?
    if (ds_inout%l_stin_evol) then
        stin_evol = ds_inout%stin_evol
        ! Get cara_elem of stin_evol
        call dismoi('CARA_ELEM', stin_evol, 'RESULTAT', repk=cara_elem_in)
        if (cara_elem_in(1:6) .ne. "#AUCUN") then
            if (cara_elem .ne. cara_elem_in) then
                call utmess('A', 'MECANONLINE_1')
            end if
        end if
    end if
!   ------------------------------------------------------------------------------------------------

! - Read parameters for algorithm management
!
    call nmdomt(ds_algopara, ds_algorom)
!
! - Read parameters for algorithm management (line search)
!
    call nmdomt_ls(ds_algopara)

! - Read objects for constitutive laws
    call nmdorc(model, mater, l_etat_init, &
                ds_constitutive%compor, ds_constitutive%carcri, ds_constitutive%mult_comp)
!
! - Read parameters for convergence
!
    call nmdocn(ds_conv, mesh)

! - Create datastructure for dynamic
    call ndcrdy(result, sddyna)

! - Get parameters from command file for dynamic
    call ndlect(model, mater, cara_elem, listLoad, &
                sddyna, nlDynaDamping)
!
! - Read parameters for post-treatment management (CRIT_STAB and MODE_VIBR)
!
    call nmdopo(sddyna, ds_posttimestep)
!
! - Read parameters for contact management
!
    call nonlinDSContactRead(ds_contact)
!
! - Read parameters for energy management
!
    call nonlinDSEnergyRead(ds_energy)
!
! - Read parameters for measure and statistic management
!
    call nonlinDSMeasureRead(ds_measure)
!
! - Read parameters for error indicator
!
    if (staticOperator) then
        call nonlinDSErrorIndicRead(ds_errorindic)
    end if
!
! - Read parameters for printing
!
    call nonlinDSPrintRead(ds_print)
!
end subroutine
