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
subroutine op0154()
!
    use mesh_module, only: checkInclude
    use mesh_modification_module, only: meshOperModiGetPara, meshOperModiDelPara, &
                                        meshOrieShell, meshOrieInterfPou
    use mesh_operators_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/abscur.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/chgref.h"
#include "asterfort/chpver.h"
#include "asterfort/conori.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/echell.h"
#include "asterfort/exipat.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/mai2a3.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/modi_alea.h"
#include "asterfort/momaba.h"
#include "asterfort/orilgm.h"
#include "asterfort/rotama.h"
#include "asterfort/symema.h"
#include "asterfort/tranma.h"
#include "asterfort/utmess.h"
#include "asterfort/vtgpld.h"
!
! --------------------------------------------------------------------------------------------------
!
! MODI_MAILLAGE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n1, n2, nbocc, iOcc, nbDime, ier, iOrieShell, iOrieInterfPou
    aster_logical :: bidim, lModiTopo
    character(len=8) :: mesh, meshReuse, dispMesh
    character(len=16) :: kbi1, kbi2, option
    character(len=19) :: geomInit, geomModi, disp
    real(kind=8) :: ltchar, pt(3), pt2(3), dir(3), angl
    real(kind=8) :: axe1(3), axe2(3), perp(3), alea
    type(MESH_OPER_MODI_PARA) :: meshOperModiPara
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
! - Check Includes
!
    call checkInclude()
!
! - Main datastructure
!
    call getvid(' ', 'MAILLAGE', scal=mesh, nbret=nbocc)
    ASSERT(nbocc .eq. 1)
    call getres(meshReuse, kbi1, kbi2)
    if (mesh .ne. meshReuse) then
        call utmess('F', 'MESH1_15')
    end if
!
! - Flag when mesh topology has been changed (new nodes, new elements)
!
    lModiTopo = ASTER_TRUE
!
! - For "ORIE_FISSURE"
!
    call getfac('ORIE_FISSURE', nbocc)
    if (nbocc .ne. 0) then
        lModiTopo = ASTER_TRUE
        call conori(mesh)
    end if
!
! - For "MODI_MAILLE"
!
    call getfac('MODI_MAILLE', nbocc)
    if (nbocc .ne. 0) then
        lModiTopo = ASTER_TRUE
        call getvtx('MODI_MAILLE', 'OPTION', iocc=1, scal=option)
        if (option .eq. 'NOEUD_QUART') then
            call momaba(mesh)
        end if
    end if
!
! - For "DEFORME"
!
    call getfac('DEFORME', nbocc)
    if (nbocc .ne. 0) then
        lModiTopo = ASTER_FALSE
        call getvr8('DEFORME', 'ALEA', iocc=1, scal=alea, nbret=n1)
        if (n1 .eq. 1) then
            call modi_alea(mesh, alea)
        else
            call getvid('DEFORME', 'DEPL', iocc=1, scal=disp, nbret=n1)
            ASSERT(n1 .eq. 1)
            call dismoi('NOM_MAILLA', disp, 'CHAM_NO', repk=dispMesh)
            if (dispMesh .ne. mesh) then
                call utmess('F', 'MESH1_1')
            end if
            call chpver('F', disp, 'NOEU', 'DEPL_R', ier)
            geomInit = mesh//'.COORDO'
            geomModi = mesh//'.COORD2'
            call vtgpld('CUMU', 1.d0, geomInit, disp, 'V', geomModi)
            call detrsd('CHAMP_GD', geomInit)
            call copisd('CHAMP_GD', 'G', geomModi, geomInit)
            call detrsd('CHAMP_GD', geomModi)
        end if
    end if
!
! - For "TRANSLATION"
!
    call getvr8(' ', 'TRANSLATION', nbval=0, nbret=n1)
    if (n1 .ne. 0) then
        lModiTopo = ASTER_FALSE
        geomInit = mesh//'.COORDO'
        bidim = ASTER_FALSE
        call getvr8(' ', 'TRANSLATION', nbval=0, nbret=nbDime)
        nbDime = -nbDime
        if (nbDime .eq. 2) then
            call getvr8(' ', 'TRANSLATION', nbval=2, vect=dir, nbret=n1)
            dir(3) = 0.d0
            bidim = ASTER_TRUE
        elseif (nbDime .eq. 3) then
            call getvr8(' ', 'TRANSLATION', nbval=3, vect=dir, nbret=n1)
            bidim = ASTER_FALSE
        else
            ASSERT(ASTER_FALSE)
        end if
        call tranma(geomInit, dir, bidim)
    end if
!
! - For "MODI_BASE"
!
    call getfac('MODI_BASE', nbocc)
    if (nbocc .ne. 0) then
        lModiTopo = ASTER_TRUE
        geomInit = mesh//'.COORDO'
        bidim = ASTER_FALSE
        call getvr8('MODI_BASE', 'VECT_X', iocc=1, nbval=0, nbret=nbDime)
        nbDime = -nbDime
        if (nbDime .eq. 2) then
            call getvr8('MODI_BASE', 'VECT_X', iocc=1, nbval=2, vect=pt, nbret=n1)
            pt(3) = 0.d0
            pt2(:) = 0.d0
            bidim = ASTER_TRUE
        elseif (nbDime .eq. 3) then
            call getvr8('MODI_BASE', 'VECT_X', iocc=1, nbval=3, vect=pt, nbret=n1)
            call getvr8('MODI_BASE', 'VECT_Y', iocc=1, nbval=3, vect=pt2, nbret=n1)
        else
            ASSERT(ASTER_FALSE)
        end if
        call chgref(geomInit, pt, pt2, bidim)
    end if
!
! - For "ROTATION"
!
    call getfac('ROTATION', nbocc)
    if (nbocc .ne. 0) then
        lModiTopo = ASTER_TRUE
        geomInit = mesh//'.COORDO'
        bidim = ASTER_FALSE
        do iOcc = 1, nbocc
            call getvr8('ROTATION', 'POIN_1', iocc=iOcc, nbval=0, nbret=nbDime)
            call getvr8('ROTATION', 'ANGLE', iocc=iOcc, scal=angl, nbret=n1)
            call getvr8('ROTATION', 'POIN_2', iocc=iOcc, nbval=0, nbret=n2)
            nbDime = -nbDime
            if (nbDime .eq. 2) then
                call getvr8('ROTATION', 'POIN_1', iocc=iOcc, nbval=2, vect=pt, nbret=n1)
                pt(3) = 0.d0
                pt2(:) = 0.d0
                dir(:) = 0.d0
                bidim = ASTER_TRUE
            elseif (nbDime .eq. 3) then
                call getvr8('ROTATION', 'POIN_1', iocc=iOcc, nbval=3, vect=pt, nbret=n1)
                if (n2 .ne. 0) then
                    call getvr8('ROTATION', 'POIN_2', iocc=iOcc, nbval=3, vect=pt2, nbret=n1)
                    dir = pt2-pt
                else
                    call getvr8('ROTATION', 'DIR', iocc=iOcc, nbval=3, vect=dir, nbret=n1)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
            call rotama(geomInit, pt, dir, angl, bidim)
        end do
    end if
!
! - For "SYMETRIE"
!
    call getfac('SYMETRIE', nbocc)
    if (nbocc .ne. 0) then
        lModiTopo = ASTER_TRUE
        geomInit = mesh//'.COORDO'
        do iOcc = 1, nbocc
            call getvr8('SYMETRIE', 'POINT', iocc=iOcc, nbval=0, nbret=nbDime)
            call getvr8('SYMETRIE', 'AXE_1', iocc=iOcc, nbval=0, nbret=n1)
            call getvr8('SYMETRIE', 'AXE_2', iocc=iOcc, nbval=0, nbret=n2)
!
!           DIM, N1, N2 = 2 OU 3 ==> IMPOSE PAR LES CATALOGUES
!           EN 2D : DIM=N1=2    , AXE_2 N'EXISTE PAS N2=0
!           EN 3D : DIM=N1=N2=3
            if (nbDime .eq. -2) then
                if (n1 .ne. nbDime) then
                    call utmess('F', 'MESH1_62')
                end if
                if (n2 .ne. 0) then
                    call utmess('A', 'MESH1_63')
                end if
                call getvr8('SYMETRIE', 'POINT', iocc=iOcc, nbval=2, vect=pt, nbret=nbDime)
                call getvr8('SYMETRIE', 'AXE_1', iocc=iOcc, nbval=2, vect=axe1, nbret=n1)
!              CONSTRUCTION DU VECTEUR PERPENDICULAIRE A Z ET AXE1
                perp(1) = -axe1(2)
                perp(2) = axe1(1)
                perp(3) = 0.0d0
            else
                if (n1 .ne. nbDime) then
                    call utmess('F', 'MESH1_62')
                end if
                if (n2 .ne. nbDime) then
                    call utmess('F', 'MESH1_64')
                end if
                call getvr8('SYMETRIE', 'POINT', iocc=iOcc, nbval=3, vect=pt, nbret=nbDime)
                call getvr8('SYMETRIE', 'AXE_1', iocc=iOcc, nbval=3, vect=axe1, nbret=n1)
                call getvr8('SYMETRIE', 'AXE_2', iocc=iOcc, nbval=3, vect=axe2, nbret=n2)
!              CONSTRUCTION DU VECTEUR PERPENDICULAIRE A AXE1 ET AXE2
                perp(1) = axe1(2)*axe2(3)-axe1(3)*axe2(2)
                perp(2) = axe1(3)*axe2(1)-axe1(1)*axe2(3)
                perp(3) = axe1(1)*axe2(2)-axe1(2)*axe2(1)
            end if
            call symema(geomInit, perp, pt)
        end do
    end if
!
! - For "ECHELLE"
!
    call getvr8(' ', 'ECHELLE', nbval=0, nbret=n1)
    if (n1 .ne. 0) then
        lModiTopo = ASTER_TRUE
        geomInit = mesh//'.COORDO'
        call getvr8(' ', 'ECHELLE', scal=ltchar, nbret=n2)
        call echell(geomInit, ltchar)
    end if

! - For "ORIE_PEAU" , "ORIE_LIGNE", "ORIE_NORM_COQUE" and "ORIE_INTERF_POU"
    call meshOperModiGetPara(mesh, meshOperModiPara)
    if (meshOperModiPara%orieShell .gt. 0) then
        do iOrieShell = 1, meshOperModiPara%orieShell
            ! "ORIE_NORM_COQUE"
            call meshOrieShell(mesh, meshOperModiPara%meshOperOrieShell(iOrieShell))
        end do
    else
        ! "ORIE_LIGNE" and "ORIE_PEAU"
        call orilgm(mesh)
    end if
    if (meshOperModiPara%orieInterfPou .gt. 0) then
        do iOrieInterfPou = 1, meshOperModiPara%OrieInterfPou
            ! "ORIE_INTERF_POU"
            call meshOrieInterfPou(mesh, meshOperModiPara%meshOperOrieInterfPou(iOrieInterfPou))
        end do
    end if
    call meshOperModiDelPara(meshOperModiPara)
!
! - For "ABSC_CURV"
!
    call getfac('ABSC_CURV', nbocc)
    if (nbocc .eq. 1) then
        lModiTopo = ASTER_TRUE
        call abscur(mesh)
    end if
!
! - On interdit MODI_MAILLAGE après DECOUPE_LAC (excepte DEFORME et TRANSLATION)
!
    if (lModiTopo) then
        call exipat(mesh, ier)
        if (ier == 1) then
            call utmess('F', 'MESH1_12')
        end if
    end if
!
! - Update parameters for modified mesh (bounding box and dimensions)
!
    call cargeo(mesh)
!
! - Check topological size of mesh
!
    call mai2a3(mesh)
!
    call jedema()
!
end subroutine
