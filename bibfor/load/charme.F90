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
subroutine charme(load, valeType)
!
    implicit none
!
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/caarei.h"
#include "asterfort/caddli.h"
#include "asterfort/caddlp.h"
#include "asterfort/caethm.h"
#include "asterfort/caethr.h"
#include "asterfort/cafaci.h"
#include "asterfort/cafond.h"
#include "asterfort/cafono.h"
#include "asterfort/cafthm.h"
#include "asterfort/cagene.h"
#include "asterfort/cagrou.h"
#include "asterfort/caimch.h"
#include "asterfort/caliag.h"
#include "asterfort/caliai.h"
#include "asterfort/calich.h"
#include "asterfort/calicp.h"
#include "asterfort/caliel.h"
#include "asterfort/calimc.h"
#include "asterfort/caliob.h"
#include "asterfort/calipj.h"
#include "asterfort/calirc.h"
#include "asterfort/caliso.h"
#include "asterfort/calyrc.h"
#include "asterfort/caprec.h"
#include "asterfort/carbe3.h"
#include "asterfort/carota.h"
#include "asterfort/caveas.h"
#include "asterfort/caveis.h"
#include "asterfort/cbchei.h"
#include "asterfort/cbonde.h"
#include "asterfort/cbondp.h"
#include "asterfort/cbpesa.h"
#include "asterfort/cbprca.h"
#include "asterfort/cbpres.h"
#include "asterfort/cbsint.h"
#include "asterfort/cbvite.h"
#include "asterfort/char_crea_neum.h"
#include "asterfort/chveno.h"
#include "asterfort/cormgi.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/initel.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnl.h"
#include "asterfort/utmess.h"
#include "asterfort/verif_affe.h"
!
    character(len=8), intent(in) :: load
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads for AFFE_CHAR_MECA_*
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: phenom = 'MECANIQUE'
    character(len=4), parameter :: phenomS = 'MECA'
    character(len=16), parameter :: command = 'AFFE_CHAR_MECA'
    character(len=16), parameter :: keywFactEnforceDOF = 'DDL_IMPO'
    integer(kind=8) :: geomDime, iret
    character(len=3) :: answer
    character(len=8) :: mesh, model
    character(len=13) :: loadDescBase
    character(len=19) :: loadLigrel
    character(len=8), pointer :: loadLigrelLgrf(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Mesh, Ligrel for model, dimension of model
    call cagene(load, command, model, mesh, geomDime)
    if (geomDime .gt. 3) then
        call utmess('A', 'CHARGES2_4')
    end if
! - Get Ligrel for load
    call lisnnl(phenom, load, loadDescBase)
    loadLigrel = loadDescBase//'.LIGRE'

! - Others loadings
    if (valeType .eq. 'REEL') then

! ----- LIAISON_INTERF
        call calimc(load)

! ----- RELA_CINE_BP
        call caprec(load, loadLigrel, mesh, model, valeType)

! ----- VITE_FACE
        call cbvite(phenom, load, mesh, valeType)

! ----- ONDE_FLUI
        call cbonde(load, mesh, valeType)

! ----- FLUX_THM_REP
        call cafthm(load, mesh, model, valeType)

! ----- ECHA_THM
        call caethm(load, mesh, model, valeType)

! ----- ECHA_THM_HR
        call caethr(load, mesh, model, valeType)

! ----- FORCE_SOL
        call caveis(load)

    else if (valeType .eq. 'COMP') then

    else if (valeType .eq. 'FONC') then

! ----- VITE_FACE
        call cbvite(phenom, load, mesh, valeType)

! ----- ONDE_PLANE
        call cbondp(load, mesh, geomDime, valeType)

! ----- FLUX_THM_REP
        call cafthm(load, mesh, model, valeType)

! ----- ECHA_THM
        call caethm(load, mesh, model, valeType)

! ----- ECHA_THM_HR
        call caethr(load, mesh, model, valeType)

    else
        ASSERT(ASTER_FALSE)

    end if
!
! --------------------------------------------------------------------------------------------------
!
!   Neumann loadings
!
! --------------------------------------------------------------------------------------------------
!
    if (valeType .eq. 'REEL') then

! ----- PRES_REP/FORCE_TUYAU
        call cbpres(load, mesh, model, geomDime, valeType)

! ----- PRE_EPSI
        call cbchei(load, mesh, model, valeType)

! ----- PRE_SIGM
        call cbsint(load, mesh)

! ----- EFFE_FOND
        call cafond(load, mesh, model, geomDime, valeType)

! ----- EVOL_CHAR
        call cbprca(phenom, load)

! ----- PESANTEUR
        call cbpesa(load, mesh, valeType)

! ----- ROTATION
        call carota(load, mesh, valeType)

! ----- VECT_ASSE
        call caveas(load)

! ----- FORCE_NODALE
        call cafono(load, loadLigrel, mesh, model, valeType)

! ----- FORCE_CONTOUR/FORCE_INTERNE/FORCE_ARETE/FORCE_FACE/FORCE_POUTRE/FORCE_COQUE
        call char_crea_neum(load, model, mesh, geomDime, valeType)

    else if (valeType .eq. 'COMP') then

! ----- FORCE_CONTOUR/FORCE_INTERNE/FORCE_ARETE/FORCE_FACE/FORCE_POUTRE/FORCE_COQUE
        call char_crea_neum(load, model, mesh, geomDime, valeType)

    else if (valeType .eq. 'FONC') then

! ----- PRES_REP/FORCE_TUYAU
        call cbpres(load, mesh, model, geomDime, valeType)

! ----- PRE_EPSI
        call cbchei(load, mesh, model, valeType)

! ----- PRE_SIGM
        call cbsint(load, mesh)

! ----- EFFE_FOND
        call cafond(load, mesh, model, geomDime, valeType)

! ----- FORCE_NODALE
        call cafono(load, loadLigrel, mesh, model, valeType)

! ----- FORCE_CONTOUR/FORCE_INTERNE/FORCE_ARETE/FORCE_FACE/FORCE_POUTRE/FORCE_COQUE
        call char_crea_neum(load, model, mesh, geomDime, valeType)

    else
        ASSERT(ASTER_FALSE)

    end if
!
! --------------------------------------------------------------------------------------------------
!
!   Kinematic conditions
!
! --------------------------------------------------------------------------------------------------
!
    if (valeType .eq. 'REEL') then

! ----- DDL_POUTRE
        call caddlp(load, mesh, model, valeType)

! ----- DDL_IMPO
        call caddli(keywFactEnforceDOF, load, mesh, model, valeType)

! ----- ARETE_IMPO
        call caarei(load, mesh, model, valeType)

! ----- FACE_IMPO
        call cafaci(load, mesh, model, valeType)

! ----- LIAISON_DDL
        call caliai(valeType, load, phenomS)

! ----- LIAISON_MAIL
        call calirc(phenom, load, model)

! ----- LIAISON_PROJ
        call calipj(load, model)

! ----- LIAISON_CYCL
        call calyrc(load, mesh, model, geomDime)

! ----- LIAISON_ELEM
        call caliel(load, model)

! ----- LIAISON_CHAMNO
        call calich(load, phenomS)

! ----- CHAMNO_IMPO
        call caimch(load)

! ----- LIAISON_RBE3
        call carbe3(load)

! ----- LIAISON_OBLIQUE
        call caliob(load, mesh, model, valeType)

! ----- LIAISON_GROUP
        call caliag(valeType, load, phenomS)

! ----- LIAISON_UNIF
        call cagrou(load, mesh, valeType, phenomS)

! ----- LIAISON_SOLIDE
        call caliso(load, mesh, model, valeType)

! ----- LIAISON_COQUE
        call calicp(load, mesh, model, valeType)

    else if (valeType .eq. 'COMP') then

! ----- DDL_IMPO
        call caddli(keywFactEnforceDOF, load, mesh, model, valeType)

! ----- LIAISON_DDL
        call caliai(valeType, load, phenomS)

    else if (valeType .eq. 'FONC') then

! ----- DDL_IMPO
        call caddli(keywFactEnforceDOF, load, mesh, model, valeType)

! ----- FACE_IMPO
        call cafaci(load, mesh, model, valeType)

! ----- LIAISON_DDL
        call caliai(valeType, load, phenomS)

! ----- LIAISON_OBLIQUE
        call caliob(load, mesh, model, valeType)

! ----- LIAISON_GROUP
        call caliag(valeType, load, phenomS)

! ----- LIAISON_UNIF
        call cagrou(load, mesh, valeType, phenomS)

! ----- LIAISON_COQUE
        call calicp(load, mesh, model, valeType)

    else
        ASSERT(ASTER_FALSE)

    end if

! - Update loads <LIGREL>
    call jeexin(loadLigrel//'.LGRF', iret)
    if (iret .ne. 0) then
        call adalig(loadLigrel)
        call cormgi('G', loadLigrel)
        call jeecra(loadLigrel//'.LGRF', 'DOCU', cval=phenomS)
        call initel(loadLigrel)
        call jeveuo(loadLigrel//'.LGRF', 'E', vk8=loadLigrelLgrf)
        call dismoi('PARTITION', model, "MODELE", repk=loadLigrelLgrf(2))
        call getvtx(' ', 'DOUBLE_LAGRANGE', scal=answer)
        if (answer .eq. 'NON') then
            loadLigrelLgrf(3) = 'LAG1'
        end if
    end if

! - Check mesh orientation (normals)
    if (valeType .ne. 'COMP') then
        call chveno(valeType, mesh, model)
    end if

! - Audit assignments
    call verif_affe(modele=model, sd=load)
!
end subroutine
