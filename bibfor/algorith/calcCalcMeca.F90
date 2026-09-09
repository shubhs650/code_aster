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
subroutine calcCalcMeca(nb_option, list_option, &
                        l_elem_nonl, &
                        listLoadZ, modelZ, caraElemZ, &
                        ds_constitutive, ds_material, ds_system, &
                        hval_incr, hval_algo, &
                        vediri, vefnod, &
                        nb_obje_maxi, obje_name, obje_sdname, nb_obje, &
                        l_pred)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/knindi.h"
#include "asterfort/medime.h"
#include "asterfort/merimo.h"
#include "asterfort/nmchex.h"
#include "asterfort/utmess.h"
#include "asterfort/vebtla.h"
#include "asterfort/vefnme.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
!
    integer(kind=8), intent(in) :: nb_option
    character(len=16), intent(in) :: list_option(:)
    aster_logical, intent(in) :: l_elem_nonl
    character(len=*), intent(in) :: listLoadZ, modelZ, caraElemZ
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: hval_incr(:), hval_algo(:)
    character(len=19), intent(in) :: vediri, vefnod
    integer(kind=8), intent(in) :: nb_obje_maxi
    character(len=16), intent(inout) :: obje_name(nb_obje_maxi)
    character(len=24), intent(inout) :: obje_sdname(nb_obje_maxi)
    integer(kind=8), intent(out) ::  nb_obje
    aster_logical, intent(in) :: l_pred
!
! --------------------------------------------------------------------------------------------------
!
! Command CALCUL
!
! Compute for mechanics
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_option        : number of options to compute
! In  list_option      : list of options to compute
! In  listLoad         : name of datastructure for list of loads
! In  model            : name of model
! In  caraElem         : name of elementary characteristics (field)
! In  l_elem_nonl      : .true. if all elements can compute non-linear options
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_material      : datastructure for material parameters
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  vediri           : name of elementary for reaction (Lagrange) vector
! In  vefnod           : name of elementary for forces vector (FORC_NODA)
! In  nb_obje_maxi     : maximum number of new objects to add
! IO  obje_name        : name of new objects to add
! IO  obje_sdname      : datastructure name of new objects to add
! Out nb_obje          : number of new objects to add
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_matr, l_nonl, l_forc_noda
    aster_logical :: l_lagr
    character(len=16) :: option
    character(len=24), parameter :: disp = "&&OP0026.DISP"
    character(len=19) :: varc_curr, disp_curr, sigm_curr, vari_curr
    character(len=19) :: vari_prev, disp_prev, sigm_prev, disp_cumu_inst
    integer(kind=8) :: iter_newt, ixfem, nb_subs_stat
    aster_logical :: l_xfem, l_macr_elem
    integer(kind=8) :: ldccvg
    character(len=19) :: ligrmo, caco3d, listLoad
    character(len=1), parameter :: jvBase = "G"
    character(len=24) :: model, caraElem
!
! --------------------------------------------------------------------------------------------------
!
    nb_obje = 0
    caco3d = '&&MERIMO.CARA_ROTAF'
    model = modelZ
    caraElem = caraElemZ
    listLoad = listLoadZ
!
! - Get LIGREL
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
!
! - Specific functionnalities
!
    call exixfe(model, ixfem)
    l_xfem = ixfem .ne. 0
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nb_subs_stat)
    l_macr_elem = nb_subs_stat .gt. 0
!
! - Name of variables
!
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varc_curr)
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
    call nmchex(hval_incr, 'VALINC', 'SIGPLU', sigm_curr)
    call nmchex(hval_incr, 'VALINC', 'VARPLU', vari_curr)
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', disp_prev)
    call nmchex(hval_incr, 'VALINC', 'SIGMOI', sigm_prev)
    call nmchex(hval_incr, 'VALINC', 'VARMOI', vari_prev)
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', disp_cumu_inst)
!
! - What we are computing
!
    l_matr = (knindi(16, 'MATR_TANG_ELEM', list_option, nb_option) .gt. 0)
    l_nonl = (knindi(16, 'MATR_TANG_ELEM', list_option, nb_option) .gt. 0) .or. &
             (knindi(16, 'COMPORTEMENT', list_option, nb_option) .gt. 0) .or. &
             (knindi(16, 'FORC_INTE_ELEM', list_option, nb_option) .gt. 0)
    l_forc_noda = (knindi(16, 'FORC_NODA_ELEM', list_option, nb_option) .gt. 0)
    l_lagr = l_matr
!
! - Some checks
!
    if (l_nonl) then
        if (.not. l_elem_nonl) then
            call utmess('F', 'CALCUL1_8')
        end if
        if (l_xfem) then
            call utmess('F', 'CALCUL1_10')
        end if
        if (l_macr_elem) then
            call utmess('F', 'CALCUL1_11')
        end if
        if (disp_prev .eq. ' ' .or. sigm_prev .eq. ' ' .or. vari_prev .eq. ' ') then
            call utmess('F', 'CALCUL1_12')
        end if
    end if
!
    if (l_forc_noda) then
        if (disp_prev .eq. ' ' .or. sigm_prev .eq. ' ') then
            call utmess('F', 'CALCUL1_13')
        end if
    end if
!
! - How we are computing
!
    option = ' '
    if (l_matr) then
        if (l_pred) then
            option = 'RIGI_MECA_TANG'
            l_nonl = ASTER_FALSE
        else
            option = 'FULL_MECA'
        end if
    else
        option = 'RAPH_MECA'
    end if
!
! - Physical dof computation
!
    if (l_nonl .or. l_matr) then
        iter_newt = 1
        call merimo(jvBase, &
                    l_xfem, l_macr_elem, &
                    model, caraElem, iter_newt, &
                    ds_constitutive, ds_material, ds_system, &
                    hval_incr, hval_algo, &
                    option, ldccvg)
        call detrsd('CHAMP', caco3d)
    end if
!
! - Lagrange dof computation
!
    if (l_lagr) then
        call medime(jvBase, 'CUMU', model, listLoad, ds_system%merigi)
        call vebtla(jvBase, model, disp_curr, listLoad, vediri)
    end if
!
! - Nodal forces
!
    if (l_forc_noda) then
        option = 'FORC_NODA'
        if (.not. l_nonl .and. .not. l_matr) then
! ----- Update displacement
            call copisd('CHAMP_GD', 'V', disp_prev, disp)
            call vtzero(disp)
            call vtaxpy(1.0, disp_prev, disp)
            call vtaxpy(1.0, disp_cumu_inst, disp)
            ! calcul avec sigma init (sans integration de comportement)
            call copisd('CHAMP_GD', 'V', sigm_prev, sigm_curr)
            call vefnme(option, model, ds_material%mateco, caraElem, &
                        ds_constitutive%compor, 0, ligrmo, &
                        varc_curr, sigm_curr, ' ', disp, &
                        jvBase, vefnod)
        else
            ! t(i) => t(i+1) : depplu => depmoi, depdel = 0
            call vefnme(option, model, ds_material%mateco, caraElem, &
                        ds_constitutive%compor, 0, ligrmo, &
                        varc_curr, sigm_curr, ' ', disp_curr, jvBase, vefnod)
        end if
    end if
!
! - New objects in table
!
    nb_obje = 0
    if (l_lagr) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'FORC_DIRI_ELEM'
        obje_sdname(nb_obje) = vediri
    end if
    if (l_nonl) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'FORC_INTE_ELEM'
        obje_sdname(nb_obje) = ds_system%veinte
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'SIEF_ELGA'
        obje_sdname(nb_obje) = sigm_curr
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'VARI_ELGA'
        obje_sdname(nb_obje) = vari_curr
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'CODE_RETOUR_INTE'
        obje_sdname(nb_obje) = ds_constitutive%comp_error
    end if
    if (l_matr) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'MATR_TANG_ELEM'
        obje_sdname(nb_obje) = ds_system%merigi
        if (l_pred) then
            nb_obje = nb_obje+1
            ASSERT(nb_obje .le. nb_obje_maxi)
            obje_name(nb_obje) = 'CODE_RETOUR_INTE'
            obje_sdname(nb_obje) = ds_constitutive%comp_error
        end if
    end if
    if (l_forc_noda) then
        nb_obje = nb_obje+1
        ASSERT(nb_obje .le. nb_obje_maxi)
        obje_name(nb_obje) = 'FORC_NODA_ELEM'
        obje_sdname(nb_obje) = vefnod
    end if
!
end subroutine
