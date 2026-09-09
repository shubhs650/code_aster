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
subroutine caarle(numeDofZ, iocc, listRelaZ, loadZ)
!
    implicit none
#include "jeveux.h"
#include "asterfort/arlcou.h"
#include "asterfort/arllec.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=*), intent(in) :: numeDofZ
    integer(kind=8), intent(in) :: iocc
    character(len=*), intent(in) :: listRelaZ, loadZ
!
! --------------------------------------------------------------------------------------------------
!
! LIAISON_ELEM
!
! For ARLEQUIN
!
! --------------------------------------------------------------------------------------------------
!
! In  numeDof          : name of numbering object (NUME_DDL)
! In  iocc             : index of factor keyword
! In  listRela         : name of object for linear relations
! In  load             : load
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = "LIAISON_ELEM"
    character(len=8) :: load
    character(len=14) :: numeDof
    character(len=19) :: listRela
    character(len=24), parameter :: typmai = '&&CAARLE.NOMTM'
    character(len=10) :: noma, nomb, nom1, nom2
    character(len=8) :: model, mesh, partModel(2), partKine(2)
    character(len=8) :: k8bid
    integer(kind=8) :: dime
    integer(kind=8) :: nbtyp
    integer(kind=8) :: ibid, i
    integer(kind=8) :: jtypm
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    numeDof = numeDofZ
    load = loadZ
    listRela = listRelaZ

! - Get model
    call getvid(' ', 'MODELE', iocc=0, nbval=1, scal=model, nbret=ibid)

! - Get mesh
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - STRUCTURES DE DONNEES
    noma = load(1:8)//'.A'
    nomb = load(1:8)//'.B'
    nom1 = load(1:8)//'.1'
    nom2 = load(1:8)//'.2'

! --- CREATION D'UN VECTEUR CONTENANT LE NOM DES TYPES DE MAILLES
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtyp, k8bid)
    call wkvect(typmai, 'V V K8', nbtyp, jtypm)
    do i = 1, nbtyp
        call jenuno(jexnum('&CATA.TM.NOMTM', i), zk8(jtypm-1+i))
    end do

! - LECTURE ET VERIFICATION DES MAILLES DES MODELES
    call arllec(factorKeyword, iocc, model, noma, nomb, &
                partModel, partKine, dime)

! - CALCUL DES EQUATIONS DE COUPLAGE
    call arlcou(mesh, iocc, model, typmai, noma, &
                nomb, partKine, dime, listRela, load)

! - Cleaning
    call jedetr(nom1)
    call jedetr(nom2)
    call jedetr(noma//'.GROUPEMA')
    call jedetr(nomb//'.GROUPEMA')
    call jedetr(typmai)
!
    call jedema()
!
end subroutine
