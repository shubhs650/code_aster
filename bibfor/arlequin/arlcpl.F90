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
subroutine arlcpl(zocc, nbma1, nbma2, mail, nomo, &
                  typmai, nom1, nom2, ndim, lisrel, &
                  charge)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/arlcp2.h"
#include "asterfort/arlcp3.h"
#include "asterfort/arlmai.h"
#include "asterfort/arlmod.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "MeshTypes_type.h"
!
    character(len=24) :: typmai
    character(len=8) :: mail, nomo
    character(len=8) :: charge
    character(len=19) :: lisrel
    character(len=10) :: nom1, nom2
    integer(kind=8) :: ndim, zocc
    integer(kind=8) :: nbma1, nbma2
!
! ----------------------------------------------------------------------
!
! ROUTINE ARLEQUIN
!
! CALCUL DES MATRICES DE COUPLAGE ELEMENTAIRE ARLEQUIN
! ASSEMBLAGE DANS LES MATRICES ARLEQUIN MORSES
!
! ----------------------------------------------------------------------
!
! IN  MAIL   : NOM DU MAILLAGE
! IN  NOMO   : NOM DU MODELE
! IN  TYPMAI : SD CONTENANT NOM DES TYPES ELEMENTS (&&CATA.NOMTM)
! IN  NOM1   : NOM DE LA SD DE STOCKAGE PREMIER GROUPE
! IN  NOM2   : NOM DE LA SD DE STOCKAGE SECOND GROUPE
! IN  NDIM   : DIMENSION DE L'ESPACE GLOBAL (2 OU 3)
!
!
    integer(kind=8) :: nbnomx
    parameter(nbnomx=MT_NNOMAX)
    integer(kind=8) :: nliai, nddl
    parameter(nliai=12, nddl=nliai*nliai)
    aster_logical :: proj
    character(len=8) :: marlel, modarl, mailar
    character(len=24) :: tabcor
    integer(kind=8) :: iaux, jaux
    integer(kind=8) :: imatu1, imatu2, iexi
    integer(kind=8) :: jma1, jma2
    integer(kind=8) :: nbnoc1, nbnoc2
    integer(kind=8) :: chtest, i, j, jj, k
    character(len=19) :: ligarl, arlmt1, arlmt2
    real(kind=8) :: m1de(nliai, nliai)
    real(kind=8) :: m3de(nliai, 3*nbnomx)
    character(len=5) :: ch1, ch2
!
    integer(kind=8) :: len1, len2, iproj
    character(len=5), dimension(2+nbnomx, nbma1) :: numno1
    character(len=5), dimension(2, nbma2) :: numno2
    character(len=5), dimension(nbnomx*nbma1) :: numn1t
    character(len=5), dimension(2*nbma2) :: numn2t
    real(kind=8) :: m3dea(12, 3*nbnomx, nbma1), m1dea(12, 12, nbma2)
!
! ----------------------------------------------------------------------
    call jemarq()
!
! --- INITIALISATIONS
!
    numno1 = '0'
    numno2 = '0'
    numn1t = '0'
    numn2t = '0'
    m1dea = 0.0
    m3dea = 0.0
!
! --------------------------------------------------------------------
!
    marlel = '&&MARLEL'
!
! --- CREATION PSEUDO-MAILLAGE
!
    call arlmai(mail, mailar, ndim, nom1, nom2, &
                tabcor, nbma1, nbma2)
!
! --- CREATION PSEUDO-MODELE
!
    call arlmod(nomo, mailar, modarl, tabcor)
!
! --- MATRICES DE COUPLAGE ELEMENTAIRES
!
    iproj = 0
    proj = .false.
    do jma2 = 1, nbma2
        do jma1 = 1, nbma1
            call arlcp2(zocc, mail, nomo, typmai, nom1, &
                        nom2, marlel, modarl, jma1, jma2, &
                        tabcor, mailar, proj)
            if (proj) then
                iproj = iproj+1
                call dismoi('NOM_LIGREL', modarl, 'MODELE', repk=ligarl)
                arlmt2 = marlel(1:8)//'.ARLMT2'
                call jeexin(jexnum(arlmt2(1:19)//'.RESL', 2), iexi)
                if (iexi == 0) goto 210
                call jeveuo(jexnum(arlmt2(1:19)//'.RESL', 2), 'L', imatu2)
                nbnoc1 = nint(zr(imatu2+nddl))
                nbnoc2 = nint(zr(imatu2+nddl+1))
                do j = 1, nbnoc2
                    call codent(nint(zr(imatu2-1+(nddl+4)+nbnoc1+j)), 'G', ch2)
                    numno2(j, jma2) = ch2
                    numno1(j, iproj) = ch2
                end do
                do j = 1, nbnoc1
                    call codent(nint(zr(imatu2-1+(nddl+4)+j)), 'G', ch1)
                    numno1(nbnoc2+j, iproj) = ch1
                end do
!
! --- RECUPERATION DE LA MATRICE DE COUPLAGE 1D-1D
!
                do iaux = 1, 6*nbnoc2
                    do jaux = 1, 6*nbnoc2
                        m1de(iaux, jaux) = zr(imatu2-1+(6*nbnoc2)*(iaux-1)+jaux)
                        m1dea(iaux, jaux, jma2) = m1de(iaux, jaux)
                    end do
                end do
!
! --- RECUPERATION DE LA MATRICE DE COUPLAGE 1D-3D
!
                arlmt1 = marlel(1:8)//'.ARLMT1'
                call jeexin(jexnum(arlmt1(1:19)//'.RESL', 1), iexi)
                if (iexi == 0) goto 210
                call jeveuo(jexnum(arlmt1(1:19)//'.RESL', 1), 'L', imatu1)
                do iaux = 1, 6*nbnoc2
                    do jaux = 1, 3*nbnoc1
                        m3de(iaux, jaux) = zr(imatu1-1+(3*nbnoc1)*(iaux-1)+jaux)
                        m3dea(iaux, jaux, iproj) = m3de(iaux, jaux)
                    end do
                end do
            end if
        end do
    end do
!
! --- CREATION DES VECTEURS NUMEROS NOEUDS POUR L'ASSEMBLAGE
!
    jj = 0
    do i = 1, nbma2
        do j = 1, nbnoc2
            chtest = 0
            do k = 1, nbnoc2*nbma2
                if (numno2(j, i) == numn2t(k)) then
                    chtest = 1
                end if
            end do
            if (chtest == 0) then
                jj = jj+1
                numn2t(jj) = numno2(j, i)
            end if
        end do
    end do
    len2 = jj
    jj = 0
    do i = 1, nbma1
        do j = 1, nbnoc1
            chtest = 0
            do k = 1, nbnoc1*nbma1
                if ((numno1(2+j, i) == numn1t(k)) .or. (numno1(2+j, i) == '0')) then
                    chtest = 1
                end if
            end do
            if (chtest == 0) then
                jj = jj+1
                numn1t(jj) = numno1(2+j, i)
            end if
        end do
    end do
    len1 = jj
!
! --- ASSEMBLAGE DES MATRICES DE COUPLAGE ELEMENTAIRES
! --- ET AFFECTATION DES RELATIONS CINEMATIQUES
!
    call arlcp3(nbma1, nbma2, numno1, numno2, m3dea, &
                m1dea, numn1t, numn2t, len1, len2, &
                lisrel, charge)
!
210 continue
!
    call jedetr(tabcor)
!
    call jedema()
!
end subroutine
