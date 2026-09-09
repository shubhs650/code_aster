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
subroutine te0451(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/excent.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/jevech.h"
#include "asterfort/plate_type.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT, DST, Q4G, COQUE_AXIS, COQUE_3D
!
! Options: EFGE_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: elasKeyword = 'ELAS_COQMU'
    integer(kind=8) :: nbLayer, npgh, jvSiefIn, idec, jvEfgeElga, npg, itab(7), iret
    integer(kind=8) :: nbsp, kpg, nbsig, nbeff, iLayer
    integer(kind=8) :: jvMaterc, elasID
    real(kind=8) :: nxx, nyy, mxx, myy, nxy, mxy, qx, qy, excen
    real(kind=8) :: cb, cm, ch, h, hb, hm, hh
    real(kind=8) :: siyyb, siyym, siyyh, sixxb, sixxm, sixxh, sixyb, sixym
    real(kind=8) :: sixyh
    real(kind=8) :: siyzb, siyzm, siyzh, sixzb, sixzm, sixzh, epcou(100), hLayer
    character(len=3) :: iLayerStr
    character(len=2) :: oneStr
    aster_logical :: lComposite, lreel
    integer(kind=8), parameter :: nbProp = 1
    character(len=16) :: propName(1)
    integer(kind=8) :: propCode(1)
    real(kind=8) :: propVale(1)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'EFGE_ELGA')

! - Get plate properties
    call getCara(plateCara, plateOrie)
    call compCoorSystNone(plateOrie)

! - Size of stress/force
    nbsig = 6
    nbeff = 8
    if (plateCara%type .eq. PLATE_COAX) then
        nbsig = 4
        nbeff = 6
    end if

! - Get shell properties
    h = plateCara%thick
    nbLayer = plateCara%nbLayer
    ASSERT(nbLayer .le. 100)

! - Detect ELAS_COQMU
    call jevech('PMATERC', 'L', jvMaterc)
    call get_elas_id(zi(jvMaterc), elasID)
    lComposite = elasID .eq. ELAS_COMPOSITE
    if (lComposite) then
        ASSERT(nbLayer .le. 100)
        call codent(1, 'G', oneStr)
        do iLayer = 1, nbLayer
! --------- Get thickness of current layer
            call codent(iLayer, 'G', iLayerStr)
            propName(1) = 'C'//iLayerStr//'_V'//oneStr
            call rcvala(zi(jvMaterc), ' ', elasKeyword, &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            ASSERT(propCode(1) .eq. 0)
            hLayer = propVale(1)
            ASSERT(hLayer .ge. 0.d0)
            epcou(iLayer) = hLayer
        end do
    end if

! - Acces to input stress in layers
    call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, itab=itab)
    jvSiefIn = itab(1)
    npg = itab(3)
    nbsp = itab(7)
    npgh = 3
    ASSERT(nbsp .eq. nbLayer*npgh)
    ASSERT(itab(2) .eq. nbsig*npg)

! - Output field
    call tecach('OOO', 'PEFGER', 'E', iret, nval=7, itab=itab)
    jvEfgeElga = itab(1)
    ASSERT(itab(2) .eq. nbeff*npg)

    do kpg = 1, npg
        nxx = 0.d0
        nyy = 0.d0
        nxy = 0.d0
        mxx = 0.d0
        myy = 0.d0
        mxy = 0.d0
        qx = 0.d0
        qy = 0.d0

        hb = -h/2
        do iLayer = 1, nbLayer
            idec = ((kpg-1)*nbLayer+ &
                    (iLayer-1))*npgh*nbsig
!
!         -- HB, HM, HH : "HAUTEUR" DES SOUS-POINTS :
            if (lComposite) then
                hLayer = epcou(iLayer)
            else
                hLayer = h/nbLayer
            end if
            hm = hb+hLayer/2.d0
            hh = hm+hLayer/2.d0
!
!         -- SIXXB, SIYYB, ... : CONTRAINTES AU BAS DE LA COUCHE
            sixxb = zr(jvSiefIn-1+idec+1)
            siyyb = zr(jvSiefIn-1+idec+2)
            sixyb = zr(jvSiefIn-1+idec+4)
            if (nbsig .eq. 6) then
                sixzb = zr(jvSiefIn-1+idec+5)
                siyzb = zr(jvSiefIn-1+idec+6)
            end if
!         -- SIXXM, SIYYM, ... : CONTRAINTES AU MILIEU DE LA COUCHE
            sixxm = zr(jvSiefIn-1+idec+1+nbsig)
            siyym = zr(jvSiefIn-1+idec+2+nbsig)
            sixym = zr(jvSiefIn-1+idec+4+nbsig)
            if (nbsig .eq. 6) then
                sixzm = zr(jvSiefIn-1+idec+5+nbsig)
                siyzm = zr(jvSiefIn-1+idec+6+nbsig)
            end if
!
!         -- SIXXH, SIYYH, ... : CONTRAINTES EN HAUT DE LA COUCHE
            sixxh = zr(jvSiefIn-1+idec+1+2*nbsig)
            siyyh = zr(jvSiefIn-1+idec+2+2*nbsig)
            sixyh = zr(jvSiefIn-1+idec+4+2*nbsig)
            if (nbsig .eq. 6) then
                sixzh = zr(jvSiefIn-1+idec+5+2*nbsig)
                siyzh = zr(jvSiefIn-1+idec+6+2*nbsig)
            end if
!
!         -- ON INTEGRE DANS L'EPAISSEUR DE CHAQUE COUCHE
!            AVEC UNE FORRMULE DE NEWTON-COTES A 3 POINTS
!            LES COEFFICIENTS SONT 1/6, 4/6 ET 1/6
            cb = hLayer/6
            cm = 4.d0*hLayer/6
            ch = hLayer/6
!
!         -- NXX, NYY, NXY = SOMME DE SIXX, SIYY, SIXY :
            nxx = nxx+cb*sixxb+cm*sixxm+ch*sixxh
            nyy = nyy+cb*siyyb+cm*siyym+ch*siyyh
            nxy = nxy+cb*sixyb+cm*sixym+ch*sixyh
!
            if (nbeff .eq. 8) then
!           -- QX, QY = SOMME DE SIXZ, SIYZ
                qx = qx+cb*sixzb+cm*sixzm+ch*sixzh
                qy = qy+cb*siyzb+cm*siyzm+ch*siyzh
            end if
!
!         -- MXX, MYY, MXY = MOMENTS DE SIXX, SIYY, SIXY :
            mxx = mxx+cb*sixxb*hb+cm*sixxm*hm+ch*sixxh*hh
            myy = myy+cb*siyyb*hb+cm*siyym*hm+ch*siyyh*hh
            mxy = mxy+cb*sixyb*hb+cm*sixym*hm+ch*sixyh*hh
!
!         -- MISE A JOUR DE HB POUR LA COUCHE SUIVANTE :
            hb = hb+hLayer
        end do
!
        zr(jvEfgeElga-1+(kpg-1)*nbeff+1) = nxx
        zr(jvEfgeElga-1+(kpg-1)*nbeff+2) = nyy
        zr(jvEfgeElga-1+(kpg-1)*nbeff+4) = mxx
        zr(jvEfgeElga-1+(kpg-1)*nbeff+5) = myy
        if (nbeff .eq. 8) then
            zr(jvEfgeElga-1+(kpg-1)*nbeff+3) = nxy
            zr(jvEfgeElga-1+(kpg-1)*nbeff+6) = mxy
            zr(jvEfgeElga-1+(kpg-1)*nbeff+7) = qx
            zr(jvEfgeElga-1+(kpg-1)*nbeff+8) = qy
        end if
    end do
!
!
!     -- POUR LES COQUES EXCENTREES, LES EFFORTS CALCULES SONT
!        DANS LE PLAN 'MOYEN'. IL FAUT LES CALCULER DANS LE PLAN 'MAIL'
!     -----------------------------------------------------------------
    if (plateCara%type .eq. PLATE_DKT .or. plateCara%type .eq. PLATE_DST) then
        excen = plateCara%offset
        lreel = .true.
        call excent('MAIL', excen, npg, nbeff, lreel, &
                    zr(jvEfgeElga), zr(jvEfgeElga), zc(jvEfgeElga), zc(jvEfgeElga))
    end if
!
end subroutine
