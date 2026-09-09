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
subroutine dxeffi(plateCara, plateOrie, &
                  option, nomte, cont, nbEfgeNd, &
                  effint)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dxdmul.h"
#include "asterfort/dxmate.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8) :: cont(*), effint(*)
    character(len=16) :: nomte
    character(len=*) :: option
    integer(kind=8) :: nbEfgeNd
!
! --------------------------------------------------------------------------------------------------
!
!     IN  NOMTE  : NOM DE L'ELEMENT TRAITE
!     IN  XYZL   : COORDONNEES DES NOEUDS
!     IN  UL     : DEPLACEMENT A L'INSTANT T
!     IN  nbEfgeNd    : =6 : 6 CMP D'EFFORT PAR NOEUD
!     IN  nbEfgeNd    : =8 : 8 CMP D'EFFORT PAR NOEUD
!     OUT EFFINT : EFFORTS INTERNES
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0, deux = 2.d0
    integer(kind=8) :: npg
    integer(kind=8) :: nbcon, nbLayer, npgh, k, ipg, iLayer, igauh, icpg
    real(kind=8) :: hLayer, h, zic, zmin, coef, distn, coehsd
    real(kind=8) :: n(3), m(3), t(2)
    integer(kind=8) :: multic, iniv
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: hm(3, 3)
    real(kind=8) :: d1i(2, 2), d2i(2, 4)
    aster_logical :: coupmf
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', npg=npg)
!
!     RECUPERATION DES OBJETS &INEL ET DES CHAMPS PARAMETRES :
!     --------------------------------------------------------
    if (nomte .ne. 'MEDKTR3 ' .and. nomte .ne. 'MEDSTR3 ' .and. nomte .ne. 'MEDKQU4 ' .and. &
        nomte .ne. 'MEDSQU4 ' .and. nomte .ne. 'MEQ4QU4 ' .and. nomte .ne. 'MET3TR3 ') then
        call utmess('F', 'ELEMENTS_34', sk=nomte)
    end if
!
    nbcon = 6
    nbLayer = plateCara%nbLayer
    if (nbLayer .le. 0) then
        call utmess('F', 'PLATE1_10')
    end if

! - Multi-layers or not ?
    multic = 0
    if (option .eq. 'FORC_NODA') then
        call dxmate(plateCara, plateOrie, &
                    'RIGI', df, dm, dmf, dc, &
                    dci, dmc, dfc, &
                    multic, coupmf)
    end if
!
!     -- GRANDEURS GEOMETRIQUES :
!     ---------------------------
    npgh = 3
    if (multic .eq. 0) then
        h = plateCara%thick
        hLayer = h/nbLayer
        distn = plateCara%offset
        zmin = -h/deux+distn
    end if
!
    call r8inir(32, zero, effint, 1)
!
!===============================================================
!     -- BOUCLE SUR LES POINTS DE GAUSS DE LA SURFACE:
!     -------------------------------------------------
    do ipg = 1, npg
        call r8inir(3, zero, n, 1)
        call r8inir(3, zero, m, 1)
        call r8inir(2, zero, t, 1)
!
        do iLayer = 1, nbLayer
            do igauh = 1, npgh
                icpg = nbcon*npgh*nbLayer*(ipg-1)+ &
                       nbcon*npgh*(iLayer-1)+ &
                       nbcon*(igauh-1)
!
                if (igauh .eq. 1) then
                    zic = zmin+(iLayer-1)*hLayer
                    coef = 1.d0/3.d0
                else if (igauh .eq. 2) then
                    zic = zmin+hLayer/2.d0+(iLayer-1)*hLayer
                    coef = 4.d0/3.d0
                else
                    zic = zmin+hLayer+(iLayer-1)*hLayer
                    coef = 1.d0/3.d0
                end if
                if (multic .gt. 0) then
                    iniv = igauh-2
                    call dxdmul(plateCara, plateOrie, &
                                .false._1, iLayer, iniv, &
                                hm, d1i, d2i, zic, hLayer)
                end if
!
!         -- CALCUL DES EFFORTS GENERALISES DANS L'EPAISSEUR (N, M ET T)
!         --------------------------------------------------------------
                coehsd = coef*hLayer/2.d0
                n(1) = n(1)+coehsd*cont(icpg+1)
                n(2) = n(2)+coehsd*cont(icpg+2)
                n(3) = n(3)+coehsd*cont(icpg+4)
                m(1) = m(1)+coehsd*zic*cont(icpg+1)
                m(2) = m(2)+coehsd*zic*cont(icpg+2)
                m(3) = m(3)+coehsd*zic*cont(icpg+4)
                t(1) = t(1)+coehsd*cont(icpg+5)
                t(2) = t(2)+coehsd*cont(icpg+6)
            end do
        end do
        do k = 1, 3
            effint((ipg-1)*nbEfgeNd+k) = n(k)
            effint((ipg-1)*nbEfgeNd+k+3) = m(k)
        end do
        if (nbEfgeNd .gt. 6) then
            effint((ipg-1)*nbEfgeNd+7) = t(1)
            effint((ipg-1)*nbEfgeNd+8) = t(2)
        end if
    end do
!
end subroutine
