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
subroutine te0408(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxtpif.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/fonbpa.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/provec.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT, DST, Q4G, COQUE_3D, COQUE_AXIS
!
! Options: PREP_VRC
!
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: itabp(8), jvTempR, iPara, nbLayer, npgh, itemps, ino
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano, isp
    integer(kind=8) :: ier, igauh, iLayer, iret, jvTempF, jvTempCR, jvGeom
    integer(kind=8) :: ptd
    real(kind=8) :: tpinf, tpmoy, tpsup, cp1, cp2, cp3, tpc, zic, zmin
    real(kind=8) :: hLayer, epais, excent, norm2
    aster_logical :: lTempReal, CasIsOk
    integer(kind=8), parameter :: nbParaMaxi = 4
    real(kind=8) :: paraVale(nbParaMaxi)
    character(len=8) :: paraName(nbParaMaxi)
    character(len=16) :: funcParaName(nbParaMaxi)
    integer(kind=8) :: jvFuncProl, nbPara, casfct
    real(kind=8) :: xyz(3), cdg(3), vect(3), vectab(3), vectcd(3)
    character(len=19) :: funcName
    character(len=24) :: funcProl
    character(len=8) ::  k8bid
    blas_int :: b_incx, b_incy, b_n
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'PREP_VRC')
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    call compCoorSystNone(plateOrie)

! - Get plate properties
    nbLayer = plateCara%nbLayer
    epais = plateCara%thick
    excent = plateCara%offset

    call jevech('PTEMPCR', 'E', jvTempCR)

    call tecach('ONN', 'PTEMPER', 'L', iret, nval=8, itab=itabp)
    if (iret .eq. 0 .or. iret .eq. 3) then
        lTempReal = ASTER_TRUE
        jvTempR = itabp(1)
        cdg = 0.d0

! ----- Temperature on layers INF, SUP, MOY
        tpinf = 0.d0; tpmoy = 0.d0; tpsup = 0.d0
        do ino = 1, nno
            call dxtpif(zr(jvTempR+3*(ino-1)), zl(itabp(8)+3*(ino-1)))
            tpmoy = tpmoy+zr(jvTempR-1+3*(ino-1)+1)/dble(nno)
            tpinf = tpinf+zr(jvTempR-1+3*(ino-1)+2)/dble(nno)
            tpsup = tpsup+zr(jvTempR-1+3*(ino-1)+3)/dble(nno)
        end do

! ----- Coefficients des polynomes de degré 2
        cp1 = tpmoy
        cp2 = (tpsup-tpinf)/epais
        cp3 = 2.d0*(tpinf+tpsup-2.d0*tpmoy)/(epais*epais)
! Dans ce cas pas de prise en compte de l'excentrement
        excent = 0.0

    else
        call tecach('ONO', 'PTEMPEF', 'L', iret, iad=jvTempF)
        ASSERT(iret .eq. 0)

! ----- Les paramètres de la fonction
        funcName = zk8(jvTempF)
        funcProl = funcName//'.PROL'
        call jeveuo(funcProl, 'L', jvFuncProl)
        call fonbpa(funcName, zk24(jvFuncProl), k8bid, nbParaMaxi, nbPara, funcParaName)

! C'est soit :
!   INST    EPAIS
!   INST    EXCEN
!   INST    X   Y   Z
        !
!   INST  EPAIS EXCEN  X   Y   Z   !!! Vive les entiers codés
!   1     2     4      8  16  32
        casfct = 0
        do iPara = 1, nbPara
            select case (funcParaName(iPara))
            case ('INST')
                casfct = casfct+1
            case ('EPAIS')
                casfct = casfct+2
            case ('EXCENT')
                casfct = casfct+4
            case ('X')
                casfct = casfct+8
            case ('Y')
                casfct = casfct+16
            case ('Z')
                casfct = casfct+32
            end select
        end do
        CasIsOk = (casfct .eq. (1+2)) .or. (casfct .eq. (1+4)) .or. (casfct .eq. (1+8+16+32))
        if (.not. CasIsOk) then
            call utmess('F', 'FONCT0_80', sk=funcName)
        end if
        lTempReal = ASTER_FALSE
        call jevech('PINST_R', 'L', itemps)
        paraName(1) = 'INST'; paraVale(1) = zr(itemps)

! ----- Compute barycenter and normal
        call jevech('PGEOMER', 'L', jvGeom)
        if (nnos .eq. 3) then
            ptd = jvGeom
            cdg(1:3) = (zr(jvGeom:jvGeom+2)+ &
                        zr(jvGeom+6:jvGeom+6+2)+ &
                        zr(jvGeom+3:jvGeom+3+2))/3.0
        else if (nnos .eq. 4) then
            ptd = jvGeom+9
            cdg(1:3) = (zr(jvGeom:jvGeom+2)+ &
                        zr(jvGeom+6:jvGeom+6+2)+ &
                        zr(jvGeom+3:jvGeom+3+2)+ &
                        zr(ptd:ptd+2))/4.0
        else
            ASSERT(ASTER_FALSE)
        end if

        vectab(1:3) = zr(jvGeom+6:jvGeom+6+2)-zr(jvGeom:jvGeom+2)
        vectcd(1:3) = zr(ptd:ptd+2)-zr(jvGeom+3:jvGeom+3+2)
        call provec(vectab, vectcd, vect)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norm2 = ddot(b_n, vect, b_incx, vect, b_incy)
        vect = vect(1:3)/sqrt(norm2)
    end if

! COQUES :
!     En thermique les coques ne sont pas excentrées
!     En mécanique elles peuvent être excentrées
!     Si calcul thermique puis mécanique ( lTempReal )
!         Les températures INF, SUP, MOY ne doivent pas tenir compte de EXCENT
!     Si champ de fonction thermique puis calcul mécanique
!         EPAIS  c'est SANS la prise en compte de EXCENT
!         EXCENT c'est AVEC la prise en compte de EXCENT
    !
! hLayer    : épaisseur d'une couche
! npgh   : nombre de points par couche
    !
    hLayer = epais/nbLayer
    npgh = 3
    !
! CALCUL DE LA TEMPERATURE SUR LES COUCHES
    zmin = -epais/2.d0
    !
    do iLayer = 1, nbLayer
        do igauh = 1, npgh
            isp = (iLayer-1)*npgh+igauh
! zic : [-epais/2 ; epais/2 ]
            if (igauh .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
            else if (igauh .eq. 2) then
                zic = zmin+(iLayer-1)*hLayer+hLayer/2.0
            else
                zic = zmin+(iLayer-1)*hLayer+hLayer
            end if
            !
            if (lTempReal) then
                tpc = cp3*zic*zic+cp2*zic+cp1
            else
                if (casfct .eq. (1+2)) then
! Avec EPAIS donc SANS excentrement
                    nbPara = 2
                    paraName(2) = 'EPAIS'; paraVale(2) = zic
                else if (casfct .eq. (1+4)) then
! Avec EXCEN donc AVEC excentrement
                    nbPara = 2
                    paraName(2) = 'EXCENT'; paraVale(2) = excent+zic
                else if (casfct .eq. (1+8+16+32)) then
! Avec X Y Z
                    nbPara = 4
                    xyz(1:3) = cdg(1:3)+vect(1:3)*(excent+zic)
                    paraName(2) = 'X'; paraVale(2) = xyz(1)
                    paraName(3) = 'Y'; paraVale(3) = xyz(2)
                    paraName(4) = 'Z'; paraVale(4) = xyz(3)
                end if
                call fointe(' ', zk8(jvTempF), nbPara, paraName, paraVale, tpc, ier)
                if (ier .ne. 0) then
                    call utmess('F', 'FONCT0_80', sk=zk8(jvTempF))
                end if
            end if
            zr(jvTempCR-1+isp) = tpc
        end do
    end do
    !
end subroutine
