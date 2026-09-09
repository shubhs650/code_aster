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
subroutine raco3d(iocc, listRelaZ, loadZ)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/alchml.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rco3d_apco3d.h"
#include "asterfort/rco3d_clcrela.h"
#include "asterfort/rco3d_crch.h"
#include "asterfort/rco3d_crealigrel.h"
#include "asterfort/rco3d_crep.h"
#include "asterfort/getelem.h"
#include "asterfort/getnode.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    integer(kind=8), intent(in) :: iocc
    character(len=*), intent(in) :: listRelaZ, loadZ
!
! --------------------------------------------------------------------------------------------------
!
! LIAISON_ELEM
!
! For Shell/Solid (3D
!
! --------------------------------------------------------------------------------------------------
!
! In  iocc             : index of factor keyword
! In  listRela         : name of object for linear relations
! In  load             : load
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = "LIAISON_ELEM"
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    character(len=19) :: modelLigrel
    character(len=19), parameter :: chmlrac = '&&RACO3D.PCACOQU.CM'
    character(len=24), parameter :: lismavo = '&&RACO3D.LMAILLES.VOL'
    character(len=24), parameter :: lisnoco = '&&RACO3D.LNOEUDS.COQ'
    character(len=24), parameter :: lismaco = '&&RACO3D.LMAILLES.COQ'
    character(len=19), parameter :: ligrel = '&&RACO3D'
    character(len=8) :: model, mesh
    integer(kind=8) :: nbmavo, nbmaco, nt_nodes
    integer(kind=8) :: nb_pairs, iret
    integer(kind=8) :: i, n1
    real(kind=8) :: epai, crig
    integer(kind=8), pointer :: list_pairs(:) => null()
    integer(kind=8) :: nbnocot, jlisnoco
    integer(kind=8), allocatable :: map_noco_pair(:, :, :)
    integer(kind=8), allocatable :: map_noco_nbnoco(:, :, :)
    integer(kind=8), allocatable :: map_noco_nbelem(:, :)
    real(kind=8), pointer ::  v_epai(:) => null()
    integer(kind=8), pointer :: list_total_no_co(:) => null()
    character(len=8) :: caraElem
    character(len=8) :: load
    character(len=19) :: listRela
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    load = loadZ
    listRela = listRelaZ
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Main parameters
    call dismoi('NOM_MODELE', load, 'CHARGE', repk=model)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', modelLigrel, 'LIGREL', repk=mesh)

! - RECUPERER COEF_RIGI_DRZ
    call getvr8(factorKeyword, 'COEF_RIGI_DRZ', iocc=iocc, scal=crig, nbret=n1)
    if (n1 .eq. 0) then
        crig = 1.d0-5
    end if

! - RECUPERER LA LISTE DES MAILLES
    call getelem(mesh, factorKeyword, iocc, 'F', lismaco, &
                 nbmaco, '_COQUE')

    call getelem(mesh, factorKeyword, iocc, 'F', lismavo, &
                 nbmavo, '_MASSIF')

! - RECUPERER LA LISTE DES NOEUDS DU BORD DE LA COQUE
    call getnode(mesh, factorKeyword, iocc, 'V', lisnoco, &
                 nbnocot, ' ', '_COQUE')
    call jeveuo(lisnoco, 'L', jlisnoco)
    !
    AS_ALLOCATE(vi=list_total_no_co, size=nbnocot)
    !
    do i = 1, nbnocot
        list_total_no_co(i) = zi(jlisnoco-1+i)
    end do

!-- RECUPERER LES EPAISSEURS

    AS_ALLOCATE(vr=v_epai, size=nbmaco)

    call getvid(factorKeyword, 'CARA_ELEM', iocc=iocc, scal=caraElem, nbret=n1)
    call rco3d_crep(caraElem, mesh, lismaco, nbmaco, v_epai)
    ! RECUPERER LE MAX POUR L APPARIEMMENT
    epai = maxval(v_epai)
    !

!-- RECUPERER LA LISTE DES PAIRES
    nb_pairs = 0
    nt_nodes = 0

    call rco3d_apco3d(mesh, lismavo, lismaco, nbmavo, nbmaco, epai, &
                      list_pairs, nb_pairs, nt_nodes)
!

!-- CONSTRUCTION DU LIGREL
!
!   2D ET 3D ARRAYs POUR ACCELERER L ACCES AUX DONNEES AU MOMENT
!   DE L ASSEMBLAGE  DES MATRICES

    allocate (map_noco_pair(9, nbnocot, nb_pairs))
    allocate (map_noco_nbnoco(9, nbnocot, nb_pairs))
    allocate (map_noco_nbelem(9, nbnocot))
    !
    call rco3d_crealigrel(ligrel, mesh, model, list_pairs, &
                          nb_pairs, nt_nodes, &
                          list_total_no_co, nbnocot, map_noco_pair, &
                          map_noco_nbelem, map_noco_nbnoco)

!   CREATION DU CHAMP D ENTREE

    call alchml(ligrel, 'LIAI_CO_3D', 'PCACOQU', 'V', chmlrac, iret, ' ')
    call rco3d_crch(ligrel, mesh, chmlrac, lismaco, nbmaco, crig, v_epai)

! - Set input fields
    nbFieldIn = 1
    lpain(nbFieldIn) = 'PGEOMER'
    lchin(nbFieldIn) = mesh//'.COORDO'

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem, &
                       cacoqueZ_=chmlrac)

! - Set output field
    lpaout(1) = 'PMATUNS'
    lchout(1) = '&&RACO3D.PMATUNS'

!-- Compute elementary matrices
    call calcul('S', 'LIAI_CO_3D', ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

!-- add the linear relations
    call rco3d_clcrela(ligrel, mesh, nb_pairs, nbnocot, &
                       list_total_no_co, map_noco_pair, map_noco_nbelem, &
                       map_noco_nbnoco, lchout(1) (1:19), listRela)

! - Clean
    call detrsd('LIGREL', ligrel)
    call detrsd('CHAM_ELEM', chmlrac)
    call detrsd('RESUELEM', '&&RACO3D.PMATUNS')
    call jedetr('&&RACO3D.LMAILLES.VOL')
    call jedetr('&&RACO3D.LMAILLES.COQ')
    call jedetr('&&RACO3D.LNOEUDS.COQ')
    AS_DEALLOCATE(vi=list_pairs)
    AS_DEALLOCATE(vi=list_total_no_co)
    AS_DEALLOCATE(vr=v_epai)
    deallocate (map_noco_pair)
    deallocate (map_noco_nbelem)
    deallocate (map_noco_nbnoco)

    call jedema()

end subroutine
