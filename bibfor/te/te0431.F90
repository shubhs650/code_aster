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
subroutine te0431(option, nomte)
!
    use Behaviour_type
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystGrid
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/cargri.h"
#include "asterfort/codere.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmco1d.h"
#include "asterfort/nmgrib.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: GRILLE_MEMBRANE / GRILLE_EXCENTRE
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8), parameter :: nbProp = 1
    character(len=16), parameter :: propName(nbProp) = (/'E'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    integer(kind=8) :: nddl, nno, npg, i, j, j1, n, m, kpg, kk, kkd, lgpg
    integer(kind=8) :: cod(9)
    integer(kind=8) :: imatuu, ipoids, idfde, jvGeom, icontm, ivarim
    integer(kind=8) :: jtab(7), jcret, ideplm, ideplp, iret
    integer(kind=8) :: ivectu, icontp, ivarip, ivarix, icontx
    real(kind=8) :: dff(2, 8), b(6, 8), p(3, 6), jac
    real(kind=8) :: dir11(3), densit, pgl(3, 3), distn
    real(kind=8) :: deps, sigm, sig, tmp, rig
    integer(kind=8) :: iinstm, iinstp
    aster_logical :: lexc, lNonLine, lLine
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=16) :: relaCpla, relaComp
    integer(kind=8) :: jvMaterc, jvCarcri
    type(Behaviour_Integ) :: BEHInteg
    type(Material_Para) :: materPara
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: typmod(2)
    blas_int :: b_incx, b_incy, b_n
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), iret)
    if (iret .eq. 1) typmod(2) = ' '
    lexc = (lteatt('MODELI', 'GRC'))

    lNonLine = (option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA') .or. &
               (option(1:10) .eq. 'RIGI_MECA_')
    lLine = option .eq. 'RIGI_MECA'

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, &
                     nno=nno, npg=npg, &
                     jpoids=ipoids, jdfde=idfde)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    if (lexc) then
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        nddl = 6
    else
        nddl = 3
    end if
    call compCoorSystGrid(pgl, plateCara, plateOrie)

! - Get input fields
    if (lLine) then
        lVect = ASTER_FALSE
        lVari = ASTER_FALSE
        lSigm = ASTER_FALSE
        lMatr = ASTER_TRUE
    else if (lNonLine) then
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PDEPLMR', 'L', ideplm)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
        lgpg = max(jtab(6), 1)*jtab(7)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PVARIMP', 'L', ivarix)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Get parameters for non-linear behaviour
    relaComp = ' '
    relaCpla = ' '
    if (lNonLine) then
! ----- Get fields for non-linear behaviour
        call jevech('PCARCRI', 'L', jvCarcri)
        call jevech('PCOMPOR', 'L', vk16=compor)

! ----- Get parameters for non-linear behaviour
        relaComp = compor(RELA_NAME)
        relaCpla = compor(PLANESTRESS)

! ----- Initialisation of behaviour datastructure
        call behaviourInit(BEHInteg)

! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm)

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  zr(iinstm), zr(iinstp), &
                                  materPara, BEHInteg)
    end if

! - Get output fields
    ivarip = 1
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if

! - PARAMETRES EN SORTIE SUPPLEMENTAIE POUR LA METHODE IMPLEX
    if (option .eq. 'RIGI_MECA_IMPLEX') then
        call jevech('PCONTXR', 'E', icontx)
        b_n = to_blas_int(npg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(icontm), b_incx, zr(icontx), b_incy)
    end if
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    cod = 0

! - LECTURE DES CARACTERISTIQUES DE GRILLE ET CALCUL DE LA DIRECTION D'ARMATURE
    call cargri(plateCara, plateOrie, &
                densit, distn, dir11)

    do kpg = 1, npg
! ----- Set main parameters for material (on point)
        call initParaPoin(kpg, ksp, materPara)

        do n = 1, nno
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

! ----- CALCUL DE LA MATRICE "B" : DEPL NODAL --> EPS11 ET DU JACOBIEN
        call nmgrib(nno, zr(jvGeom), dff, dir11, lexc, &
                    plateOrie%gridNorm, b, jac, p)

        if (lLine) then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, &
                        materPara%schemePara%ksp, &
                        '+', &
                        materPara%jvMaterCode, &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            rig = propVale(1)

        else if (lNonLine) then
            sigm = zr(icontm+kpg-1)

! --------- CALCUL DE LA DEFORMATION DEPS11
            deps = 0.d0
            do i = 1, nno
                do j = 1, nddl
                    deps = deps+b(j, i)*zr(ideplp+(i-1)*nddl+j-1)
                end do
            end do

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! --------- Integrator
            call nmco1d(BEHInteg, &
                        relaComp, &
                        option, deps, sigm, &
                        zr(ivarim+(kpg-1)*lgpg), sig, &
                        zr(ivarip+(kpg-1)*lgpg), rig, cod(kpg))

            if (cod(kpg) .eq. 1) goto 900
            if (lSigm) then
                zr(icontp+kpg-1) = sig
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! --- RANGEMENT DES RESULTATS
        if (lVect) then
            do n = 1, nno
                do i = 1, nddl
                    zr(ivectu+(n-1)*nddl+i-1) = zr(ivectu+(n-1)*nddl+i-1)+ &
                                                b(i, n)*sig*zr(ipoids+kpg-1)*jac*densit
                end do
            end do
        end if
!
        if (lMatr) then
            do n = 1, nno
                do i = 1, nddl
                    kkd = (nddl*(n-1)+i-1)*(nddl*(n-1)+i)/2
                    do j = 1, nddl
                        do m = 1, n
                            if (m .eq. n) then
                                j1 = i
                            else
                                j1 = nddl
                            end if
                            tmp = b(i, n)*rig*b(j, m)*zr(ipoids+kpg-1)*jac*densit
                            if (j .le. j1) then
                                kk = kkd+nddl*(m-1)+j
                                zr(imatuu+kk-1) = zr(imatuu+kk-1)+tmp
                            end if
                        end do
                    end do
                end do
            end do
        end if
    end do
!
900 continue
    if (lSigm) then
        call codere(cod, npg, zi(jcret))
    end if
!
end subroutine
