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
subroutine op0030()
!
    use loadMecaDefinition_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/asmpi_comm.h"
#include "asterc/getres.h"
#include "asterfort/adalig.h"
#include "asterfort/aflrch.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/caform.h"
#include "asterfort/calico.h"
#include "asterfort/caliun.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/check_model.h"
#include "asterfort/chveno.h"
#include "asterfort/copisd.h"
#include "asterfort/cormgi.h"
#include "asterfort/defContactCreateObjects.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/infmaj.h"
#include "asterfort/initel.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/lgtlgr.h"
#include "asterfort/utmess.h"
#include "Contact_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! COMMANDE:  DEFI_CONTACT
!
! --------------------------------------------------------------------------------------------------
!
    mpi_int :: nb_proc, mpicou
    integer(kind=8) :: iret, geomDime, iexi
    character(len=24) :: phenomenon
    character(len=4), parameter :: valeType = 'REEL'
    character(len=8) :: mesh, model, sdcont, partit
    character(len=16) :: k16dummy
    character(len=19), parameter :: slavElemLigr = '&&OP0030.LIGRET', ligrelTmp = '&&OP0030.LIGREL'
    character(len=19) :: contLigrel
    integer(kind=8) :: contForm, algo_cont
    aster_logical :: lallv
    character(len=24) :: sdcont_defi
    aster_logical :: lLineRela
    character(len=19) :: listRela
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Initializations
    contForm = CONT_FORM_UNDEF
    call asmpi_comm('GET', mpicou)
    call asmpi_info(mpicou, size=nb_proc)

! - Get output datastructure
    call getres(sdcont, k16dummy, k16dummy)

! - Main datastructure on contact definition
    sdcont_defi = sdcont(1:8)//'.CONTACT'

! - Get model
    call getvid(' ', 'MODELE', scal=model)

! - Mesh
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    if (isParallelMesh(mesh)) then
        call utmess('F', 'CONTACT_30')
    end if

! - Checks
    call dismoi('PHENOMENE', model, 'MODELE', repk=phenomenon)
    if (phenomenon .ne. 'MECANIQUE') then
        call utmess('F', 'CONTACT_64')
    end if

! - Dimension of problem
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)

! - Get contact formulation
    call caform(contForm)

! - Check mesh (for LAC method)
    call check_model(mesh, contForm)

! - Create general datastructure
    call defContactCreateObjects(sdcont)

! - Read and create datastructures
    lLineRela = ASTER_FALSE
    listRela = " "
    if (contForm .eq. CONT_FORM_UNIL) then
        call caliun(sdcont, mesh, model)
    else
        call calico(sdcont, mesh, model, &
                    geomDime, contForm, &
                    slavElemLigr, lLineRela, listRela)
    end if

! - Create linear relations for QUAD8 elements (discrete contact)
    if (lLineRela) then
        call creaLoadObje("G", sdcont, model)
        call aflrch(listRela, sdcont, 'NLIN')
    end if

! - LIGREL of virtual cells
    if (lLineRela) then
        contLigrel = sdcont//'.CHME.LIGRE'
    else
        contLigrel = sdcont//'.CONT.LIGRE'
    end if

! - MPI forbidden for some methods (issue25897)
    if ((contForm .eq. 1) .or. (contForm .eq. 4)) then
        algo_cont = cfdisi(sdcont_defi, 'ALGO_CONT')
        if (nb_proc .gt. 1 .and. algo_cont .ne. 2) then
            call dismoi('PARTITION', model//'.MODELE', 'LIGREL', repk=partit)
            call exisd('PARTITION', partit, iexi)
            if (iexi .ne. 0) then
                call utmess('F', 'CONTACT3_45')
            end if
        end if
    end if

! - New <LIGREL>
    lallv = cfdisl(sdcont_defi, 'ALL_VERIF')
    if (contForm .eq. 2 .or. contForm .eq. 5) then
        if (.not. lallv) then
            call lgtlgr('V', slavElemLigr, ligrelTmp)
            call detrsd('LIGRET', slavElemLigr)
            call copisd('LIGREL', 'G', ligrelTmp, contLigrel)
            call detrsd('LIGREL', ligrelTmp)
        end if
    end if

! - Update <LIGREL>
    call jeexin(contLigrel//'.LGRF', iret)
    if (iret .ne. 0) then
        call adalig(contLigrel)
        call cormgi('G', contLigrel)
        call jeecra(contLigrel//'.LGRF', 'DOCU', cval='MECA')
        call initel(contLigrel)
    end if

! - Check mesh orientation (normals)
    if ((contForm .eq. 1) .or. (contForm .eq. 2) .or. (contForm .eq. 5)) then
        call chveno(valeType, mesh, model)
    end if
!
    call jedema()
!
end subroutine
