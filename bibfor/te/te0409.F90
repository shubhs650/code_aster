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
! aslint: disable=W1501
!
subroutine te0409(option, nomte)
!
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate, &
                                isPlateTria, isPlateQuad, isPlateQ4GG, isPlateDKTG
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/coqgth.h"
#include "asterfort/coqrep.h"
#include "asterfort/dhrc_lc.h"
#include "asterfort/dhrc_recup_mate.h"
#include "asterfort/dkqbf.h"
#include "asterfort/dktbf.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/dxqloc.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/glrc_lc.h"
#include "asterfort/glrc_recup_mate.h"
#include "asterfort/glrcmm.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/kit_glrc_dm_vmis.h"
#include "asterfort/maglrc.h"
#include "asterfort/pmrvec.h"
#include "asterfort/q4gbc.h"
#include "asterfort/r8inir.h"
#include "asterfort/t3gbc.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: Q4GG / DKTG
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
    integer(kind=8) :: codret
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: ul(6, 4), dul(6, 4)
    real(kind=8) :: matloc((6*4)*(6*4+1)/2), vecloc(6, 4)
!  ----------------
!  CMPS DE DEPLACEMENT :
!   - MEMBRANE : DX(N1), DY(N1), DX(N2), ..., DY(NNO)
!   - FLEXION  : DZ(N1), BETAX(N1), BETAY(N1), DZ(N2), ..., BETAY(NNO)
!  CMPS DE DEFORMATION COQUE :
!   - MEMBRANE : EPSIXX,EPSIYY,2*EPSIXY
!   - FLEXION  : KHIXX,KHIYY,2*KHIXY
!  CMPS D' EFFORTS COQUE :
!   - MEMBRANE : NXX,NYY,NXY
!   - FLEXION  : MXX,MYY,MXY
!   - CISAILLEMENT : QX , QY
! --------------------------------------------------------------------
!            NPG:    NOMBRE DE POINTS DE GAUSS PAR ELEMENT
!            NC :    NOMBRE DE COTES DE L'ELEMENT
!
    real(kind=8) :: poids
!            POIDS:  POIDS DE GAUSS (Y COMPRIS LE JACOBIEN)
!            AIRE:   SURFACE DE L'ELEMENT
    real(kind=8) :: um(2, 4), uf(3, 4), dum(2, 4), duf(3, 4)
!            UM:     DEPLACEMENT (MEMBRANE) "-"
!            UF:     DEPLACEMENT (FLEXION)  "-"
!           DUM:     INCREMENT DEPLACEMENT (MEMBRANE)
!           DUF:     INCREMENT DEPLACEMENT (FLEXION)
    real(kind=8) :: eps(3), khi(3), gam(2), deps(6), dkhi(3), dgam(2), n(3)
    real(kind=8) :: m(3), q(2)
!            EPS:    DEFORMATION DE MEMBRANE "-"
!            DEPS:   INCREMENT DE DEFORMATION DE MEMBRANE
!            KHI:    DEFORMATION DE FLEXION  "-" (COURBURE)
!            DKHI:   INCREMENT DE DEFORMATION DE FLEXION (COURBURE)
!            GAM :   DISTORSIONS TRANSVERSES
!            DGAM :  INCREMENT DES DISTORSIONS TRANSVERSES
!            N  :    EFFORT NORMAL "+"
!            M  :    MOMENT FLECHISSANT "+"
!            Q  :    CISAILLEMENT TRANSVERSE
!
    real(kind=8) :: sigmEnerm(32), efform(32), effint(32), efforp(32)
!            SIGMAM : EFFORTS DANS LE REPERE UTILISATEUR  A T "-"
!            EFFORM : EFFORTS DANS LE REPERE DE L'ELEMENT A T "-"
!            EFFINT : EFFORTS DANS LE REPERE DE L'ELEMENT A T "+"
!            EFFORP : EFFORTS DANS LE REPERE UTILISATEUR  A T "+"
!
    real(kind=8) :: df(9), dm(9), dmf(9), dc(2, 2)
!            DF :    MATRICE DE RIGIDITE TANGENTE MATERIELLE  (FLEXION)
!            DM :    MATRICE DE RIGIDITE TANGENTE MATERIELLE  (MEMBRANE)
!            DMF:    MATRICE DE RIGIDITE TANGENTE MATERIELLE  (COUPLAGE)
!            DC:     MATRICE DE RIGIDITE ELASTIQUE MATERIELLE
!                                                         (CISAILLEMENT)
!
    real(kind=8) :: dff(3, 3), dmm(3, 3), dmff(3, 3)
    real(kind=8) :: dcc(2, 2), dci(2, 2), dmc(3, 2), dfc(3, 2)
!
    real(kind=8) :: bf(3, 3*4), bm(3, 2*4), bmq(2, 3), bc(2, 3*4)
!            BF :    MATRICE "B" (FLEXION)
!            BM :    MATRICE "B" (MEMBRANE)
!            BC :    MATRICE "B" (CISAILLEMENT)
    real(kind=8) :: flex(3*4, 3*4), memb(2*4, 2*4), flexi(3*4, 3*4)
    real(kind=8) :: mefl(2*4, 3*4), work(3, 3*4)
!           MEMB:    MATRICE DE RIGIDITE DE MEMBRANE
!           FLEX:    MATRICE DE RIGIDITE DE FLEXION
!           WORK:    TABLEAU DE TRAVAIL
!           MEFL:    MATRICE DE COUPLAGE MEMBRANE-FLEXION
!           LE MATERIAU EST SUPPOSE HOMOGENE
!
    aster_logical :: leul, lLinear, lNonLine
    aster_logical :: lVect, lMatr, lVari, lSigm
    aster_logical :: coupmf
    aster_logical :: is_param_opt(2)
!
    integer(kind=8) :: ndim, nno, npg, ipoids, icoopg
    integer(kind=8) :: iret, icontm, ivarim, jvGeom, jvCarcri, jvDispM, jvDispP
    integer(kind=8) :: icontp, ivarip, ino, nbcont, ivectu, jcret, imatuu
    integer(kind=8) :: nbvari, kpg
    integer(kind=8) :: i, i1, i2, j, k, l
    integer(kind=8) :: icpg, icpv
    integer(kind=8) :: jtab(7), nbsig
    integer(kind=8) :: multic
!
    real(kind=8) :: matrElas(6, 6), dsidep(6, 6)
    real(kind=8) :: lambda, deuxmu, deumuf, lamf, gt, gc, gf, seuil, alphaf
    real(kind=8) :: excen

!   -- variables pour dhrc
    real(kind=8) :: a0(6, 6), c0(2, 2, 2)
    real(kind=8) :: aa_t(6, 6, 2), ab_(6, 2, 2), ac_(2, 2, 2), aa_c(6, 6, 2)
    real(kind=8) :: ga_t(6, 6, 2), gb_(6, 2, 2), gc_(2, 2, 2), ga_c(6, 6, 2)
    real(kind=8) :: cstseu(6)
!   -- attention la taille de ecp depend du nombre de variable interne
!   -- lors de l ajout de variable interne il faut incrementer ecr et ecrp
    integer(kind=8), parameter :: nbVariMax = 36
    real(kind=8) :: epst(6), ep, surfgp, sig(8), dsig(8), ecr(nbVariMax), ecrp(nbVariMax)
    real(kind=8) :: epsm(6), qsi, eta, ctor
    real(kind=8) :: carat3(21), jacob(5), caraq4(25)
    real(kind=8) :: matr(50), sigm(8), alfmc
    real(kind=8) :: epsi_c, epsi_els, epsi_lim, val_param_opt(10)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relaPlas, relaComp, defoComp
    type(Material_Para) :: materPara
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    pgl = 0.0
    xyzl = 0.0d0
    dff = 0.0d0
    dmm = 0.0d0
    dmff = 0.0d0
    dcc = 0.0d0
    dci = 0.0d0
    dmc = 0.0d0
    dfc = 0.0d0
    bf = 0.d0
    bm = 0.d0
    bmq = 0.d0
    bc = 0.d0
    flex = 0.d0
    memb = 0.d0
    flexi = 0.d0
    mefl = 0.d0
    work = 0.d0
    multic = -1
    codret = 0

    coupmf = ASTER_FALSE
!
    lLinear = option .eq. 'RIGI_MECA'
    lNonLine = .not. lLinear
!
    call elrefe_info(fami=fami, &
                     ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    ASSERT(isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara))

    nbsig = 6
    if (isPlateQ4GG(plateCara)) then
        nbsig = 8
    end if

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Initializations of material parameters on current cell
    call setMaterPara(fami, materPara)

! - No definition of coordinate system for anisotropic material
    call initLCSNone(materPara)

! - Behaviour
    leul = ASTER_FALSE
    lMatr = ASTER_TRUE
    lVect = ASTER_FALSE
    lVari = ASTER_FALSE
    lSigm = ASTER_FALSE
    codret = 0
    if (lNonLine) then
! ----- Get fields for non-linear behaviour
        call jevech('PVARIMR', 'L', ivarim)
        call tecach('OOO', 'PCONTMR', 'L', iret, nval=7, itab=jtab)
        icontm = jtab(1)
        ASSERT(npg .eq. jtab(3))
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', jvCarcri)

! ----- Properties of behaviour
        relaComp = compor(RELA_NAME)
        defoComp = compor(DEFO)
        leul = defoComp .eq. 'GROT_GDEP'
        read (compor(NVAR), '(I16)') nbvari
        ASSERT(nbvari .le. nbVariMax)
        if (defoComp(6:10) .eq. '_REAC') then
            call utmess('F', 'PLATE1_6')
        end if

! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

    else
        nbvari = 0
        jvCarcri = 1
        ivarim = 1
        icontm = 1
        relaComp = materPara%elasKeyword
    end if

! - Update geometry
    if (lNonLine) then
        call jevech('PDEPLMR', 'L', jvDispM)
        call jevech('PDEPLPR', 'L', jvDispP)
        if (leul) then
            do i = 1, nno
                i1 = 3*(i-1)
                i2 = 6*(i-1)
                zr(jvGeom+i1) = zr(jvGeom+i1)+zr(jvDispM+i2)+zr(jvDispP+i2)
                zr(jvGeom+i1+1) = zr(jvGeom+i1+1)+zr(jvDispM+i2+1)+zr(jvDispP+i2+1)
                zr(jvGeom+i1+2) = zr(jvGeom+i1+2)+zr(jvDispM+i2+2)+zr(jvDispP+i2+2)
            end do
        end if
    end if

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - Change coordinates of displacements
    ul = 0.d0
    dul = 0.d0
    if (lNonLine) then
        call utpvgl(nno, 6, pgl, zr(jvDispM), ul)
        call utpvgl(nno, 6, pgl, zr(jvDispP), dul)
    end if

! - COQUE HOMOGENEISEE ?
    if (lNonLine) then
        if ((materPara%elasID .eq. ELAS_SHELL) .or. &
            (materPara%elasID .eq. ELAS_COMPOSITE) .or. &
            (materPara%elasID .eq. ELAS_ORTH)) then
            call utmess('F', 'PLATE1_5')
        end if
    end if

    call r8inir(8, 0.d0, sig, 1)
    call r8inir(8, 0.d0, dsig, 1)

! - Get plate properties
    ep = plateCara%thick
    ctor = plateCara%coefRigiDRZ
    excen = plateCara%offset
!
    icontp = icontm
    ivarip = ivarim
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
!
!   -- GRANDEURS GEOMETRIQUES :
    if (nno .eq. 3) then
        call gtria3(xyzl, carat3)
    else if (nno .eq. 4) then
        call gquad4(xyzl, caraq4)
    end if
!
!   -- MISES A ZERO :
    if (lMatr) then
        call r8inir((3*nno)*(3*nno), 0.d0, flexi, 1)
        call r8inir((3*nno)*(3*nno), 0.d0, flex, 1)
        call r8inir((2*nno)*(2*nno), 0.d0, memb, 1)
        call r8inir((2*nno)*(3*nno), 0.d0, mefl, 1)
        call r8inir((6*4)*(6*4+1)/2, 0.d0, matloc, 1)
    end if
!
    if (lVect) then
        call r8inir(6*nno, 0.d0, vecloc, 1)
    end if
    call r8inir(32, 0.d0, effint, 1)
    call r8inir(32, 0.d0, efforp, 1)
!
    call r8inir(36, 0.d0, matrElas, 1)
    call r8inir(32, 0.d0, sigmEnerm, 1)
    call r8inir(32, 0.d0, efform, 1)
    call r8inir(2*4, 0.d0, um, 1)
    call r8inir(3*4, 0.d0, uf, 1)
    call r8inir(2*4, 0.d0, dum, 1)
    call r8inir(3*4, 0.d0, duf, 1)
!
!   -- PARTITION DU DEPLACEMENT EN MEMBRANE/FLEXION :
    do ino = 1, nno
        um(1, ino) = ul(1, ino)
        um(2, ino) = ul(2, ino)
        uf(1, ino) = ul(3, ino)
        uf(2, ino) = ul(5, ino)
        uf(3, ino) = -ul(4, ino)
        dum(1, ino) = dul(1, ino)
        dum(2, ino) = dul(2, ino)
        duf(1, ino) = dul(3, ino)
        duf(2, ino) = dul(5, ino)
        duf(3, ino) = -dul(4, ino)
    end do
!
!       -- CONTRAINTE 2D : NXX,NYY,NXY,MXX,MYY,MXY,QX,QY
    nbcont = 8
!
    if (lNonLine) then
! ----- Stress from global to parametric frame
        do i = 1, nbcont*npg
            sigmEnerm(i) = zr(icontm-1+i)
        end do
        call dxefro(npg, plateOrie%t2ui, sigmEnerm, efform)
    end if

! - Loop on Gauss point
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        if (lNonLine) then
            call initParaPoin(kpg, ksp, materPara)
        end if

        call r8inir(nbVariMax, 0.d0, ecrp, 1)
        call r8inir(3, 0.d0, n, 1)
        call r8inir(3, 0.d0, m, 1)
        call r8inir(2, 0.d0, q, 1)
        call r8inir(9, 0.d0, df, 1)
        call r8inir(9, 0.d0, dm, 1)
        call r8inir(9, 0.d0, dmf, 1)
        call r8inir(4, 0.d0, dc, 1)
!
        qsi = zr(icoopg-1+ndim*(kpg-1)+1)
        eta = zr(icoopg-1+ndim*(kpg-1)+2)
!
        icpg = (kpg-1)*nbcont
        icpv = (kpg-1)*nbvari
!
        if (isPlateDKTG(plateCara)) then
            if (isPlateTria(plateCara)) then
                call dxtbm(carat3(9), bm)
                call dktbf(qsi, eta, carat3, bf)
                poids = zr(ipoids+kpg-1)*carat3(7)
            elseif (isPlateQuad(plateCara)) then
                call jquad4(xyzl, qsi, eta, jacob)
                call dxqbm(qsi, eta, jacob(2), bm)
                call dkqbf(qsi, eta, jacob(2), caraq4, bf)
                poids = zr(ipoids+kpg-1)*jacob(1)
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (isPlateQ4GG(plateCara)) then
            if (isPlateTria(plateCara)) then
                call dxtbm(carat3(9), bm)
                call dstbfb(carat3(9), bf)
                call t3gbc(xyzl, qsi, eta, bc)
                poids = carat3(8)
            elseif (isPlateQuad(plateCara)) then
                call jquad4(xyzl, qsi, eta, jacob)
                call dxqbm(qsi, eta, jacob(2), bm)
                call dsqbfb(qsi, eta, jacob(2), bf)
                call q4gbc(qsi, eta, jacob(2), caraq4, bc)
                poids = zr(ipoids+kpg-1)*jacob(1)
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

        call r8inir(3, 0.d0, eps, 1)
        call r8inir(6, 0.d0, deps, 1)
        call r8inir(6, 0.d0, epsm, 1)
        call r8inir(3, 0.d0, khi, 1)
        call r8inir(3, 0.d0, dkhi, 1)
        call pmrvec('ZERO', 3, 2*nno, bm, um, eps)
        call pmrvec('ZERO', 3, 2*nno, bm, dum, deps)
        call pmrvec('ZERO', 3, 3*nno, bf, uf, khi)
        call pmrvec('ZERO', 3, 3*nno, bf, duf, dkhi)
        if (isPlateQ4GG(plateCara)) then
            call r8inir(2, 0.d0, gam, 1)
            call r8inir(2, 0.d0, dgam, 1)
            call pmrvec('ZERO', 2, 3*nno, bc, uf, gam)
            call pmrvec('ZERO', 2, 3*nno, bc, duf, dgam)
        end if

!           -- EULER_ALMANSI - TERMES QUADRATIQUES
        if (leul) then
            call r8inir(6, 0.d0, bmq, 1)
            do i = 1, 2
                do k = 1, nno
                    do j = 1, 2
                        bmq(i, j) = bmq(i, j)+bm(i, 2*(k-1)+i)*dum(j, k)
                    end do
                    bmq(i, 3) = bmq(i, 3)+bm(i, 2*(k-1)+i)*duf(1, k)
                end do
            end do
            do k = 1, 3
                do i = 1, 2
                    deps(i) = deps(i)-0.5d0*bmq(i, k)*bmq(i, k)
                end do
                deps(3) = deps(3)-bmq(1, k)*bmq(2, k)
            end do
        end if
!
        call r8inir(50, 0.d0, matr, 1)
!
        if (lNonLine) then
            do i = 1, 3
                deps(i+3) = dkhi(i)
                epsm(i) = eps(i)
                epsm(i+3) = khi(i)
            end do
            do i = 1, nbsig
                sig(i) = efform(icpg+i)
                sigm(i) = sig(i)
            end do
        end if
!
!EXCENTREMENT RAJOUTE UN COUPLAGE MEMBRANE_FLEXION : EPSI = EPSI+EXCENT*KHI
        if (relaComp(1:4) .ne. 'ELAS' .and. excen .gt. 0.0d0) then
            do i = 1, 3
                epsm(i) = epsm(i)+excen*khi(i)
                deps(i) = deps(i)+excen*khi(i)
            end do
        end if

        if (relaComp(1:4) .eq. 'ELAS') then
            call r8inir(3*3, 0.d0, dff, 1)
            call r8inir(3*3, 0.d0, dmm, 1)
            call r8inir(3*3, 0.d0, dmff, 1)
            call r8inir(2*2, 0.d0, dcc, 1)
            call r8inir(2*2, 0.d0, dci, 1)
            call r8inir(3*2, 0.d0, dmc, 1)
            call r8inir(3*2, 0.d0, dfc, 1)
            call dxmate(plateCara, plateOrie, &
                        fami, dff, dmm, dmff, dcc, &
                        dci, dmc, dfc, &
                        multic, coupmf)

            call r8inir(36, 0.d0, dsidep, 1)
            do i = 1, 3
                do j = 1, 3
!                       -- MEMBRANE
                    dsidep(i, j) = dmm(i, j)
!                       -- FLEXION
                    dsidep(i+3, j+3) = dff(i, j)
!                       -- MEMBRANE-FLEXION
                    dsidep(i, j+3) = dmff(i, j)
                    dsidep(j+3, i) = dmff(i, j)
                end do
            end do
!
!EXCENTREMENT RAJOUTE UN COUPLAGE MEMBRANE_FLEXION : EPSI = EPSI+EXCENT*KHI
! ON NE FAIT RIEN CAR LE COUPLAGE MEMBRANE_FLEXION INDUIT EST PRIS EN COMPTE DANS DXMATE
!               --  prise en compte de la dilatation thermique
            if (lNonLine) then
                call coqgth(materPara, ep, epsm, deps)
            end if
!
!               -- calcul de l'accroissement de contrainte
            do i = 1, 6
                dsig(i) = 0.d0
                do j = 1, 6
                    dsig(i) = dsig(i)+dsidep(i, j)*deps(j)
                end do
            end do
!               -- calcul de l'accoissement effort cisaillement
            if (isPlateQ4GG(plateCara)) then
                dsig(7) = dcc(1, 1)*dgam(1)
                dsig(8) = dcc(2, 2)*dgam(2)
            end if
!
            do i = 1, nbsig
                sig(i) = sig(i)+dsig(i)
            end do
!
        else if (relaComp(1:11) .eq. 'GLRC_DAMAGE') then
            do i = 1, nbvari
                ecr(i) = zr(ivarim-1+icpv+i)
            end do
            call maglrc(plateCara, plateOrie, materPara%jvMaterCode, &
                        matr, matrElas, ecr)
            if (isPlateQ4GG(plateCara)) then
                dcc(1, 1) = matr(14)
                dcc(2, 2) = matr(15)
                dcc(1, 2) = 0.d0
                dcc(2, 1) = 0.d0
            end if
!
!               -- aire de surface appartenant au point de g.
            surfgp = poids
!
!               --  prise en compte de la dilatation thermique
            call coqgth(materPara, ep, epsm, deps)
            do i = 1, 6
                epst(i) = epsm(i)+deps(i)
            end do
!
            call r8inir(36, 0.d0, dsidep, 1)
            call glrcmm(materPara%jvMaterCode, &
                        matr, ep, surfgp, pgl, &
                        epst, deps, dsig, ecr, matrElas, &
                        dsidep, zr(jvCarcri), codret)
!
            do i = 1, 3
                dsig(i) = dsig(i)*ep
                dsig(3+i) = dsig(3+i)*ep*ep/6.d0
            end do
!
            do i = 1, nbvari
                ecrp(i) = ecr(i)
            end do
!
            if (isPlateQ4GG(plateCara)) then
                dsig(7) = dcc(1, 1)*dgam(1)
                dsig(8) = dcc(2, 2)*dgam(2)
            end if
!
            do i = 1, nbsig
                sig(i) = sig(i)+dsig(i)
            end do
!
        else if (relaComp(1:7) .eq. 'GLRC_DM') then
            if (lNonLine) then
                do i = 1, nbvari
                    ecr(i) = zr(ivarim-1+icpv+i)
                end do
            end if
            call glrc_recup_mate(materPara%jvMaterCode, relaComp, lLinear, ep, lambda=lambda, &
                                 deuxmu=deuxmu, lamf=lamf, deumuf=deumuf, &
                                 gt=gt, gc=gc, gf=gf, seuil=seuil, &
                                 alpha=alphaf, alfmc=alfmc, epsic=epsi_c, &
                                 epsiels=epsi_els, epsilim=epsi_lim, &
                                 is_param_opt_=is_param_opt, val_param_opt_=val_param_opt)

!
!               --  prise en compte de la dilatation thermique
            call coqgth(materPara, ep, epsm, deps)
!
!               -- endommagement seulement
            call r8inir(36, 0.d0, dsidep, 1)
            call glrc_lc(plateOrie, &
                         epsm, deps, ecr, option, sig, &
                         ecrp, dsidep, lambda, deuxmu, lamf, &
                         deumuf, gt, gc, gf, seuil, &
                         alphaf, alfmc, zr(jvCarcri), &
                         epsi_c, epsi_els, epsi_lim, codret, &
                         ep, is_param_opt, val_param_opt)
!
        else if (relaComp(1:4) .eq. 'DHRC') then
!
            if (lNonLine) then
                do i = 1, nbvari
                    ecr(i) = zr(ivarim-1+icpv+i)
                end do
            end if
!
            call dhrc_recup_mate(materPara%jvMaterCode, relaComp, a0, c0, &
                                 aa_t, ga_t, ab_, gb_, ac_, &
                                 gc_, aa_c, ga_c, cstseu)

!
!               --  prise en compte de la dilatation thermique
            call coqgth(materPara, ep, epsm, deps)
!
!               -- endommagement couple glissement acier beton
            call r8inir(36, 0.d0, dsidep, 1)
            call dhrc_lc(plateOrie, epsm, deps, ecr, option, &
                         sig, ecrp, a0, c0, aa_t, &
                         ga_t, ab_, gb_, ac_, gc_, &
                         aa_c, ga_c, cstseu, zr(jvCarcri), codret, &
                         dsidep, ASTER_FALSE)
!
        else if (relaComp(1:7) .eq. 'KIT_DDI') then
!
            if (lNonLine) then
                do i = 1, nbvari
                    ecr(i) = zr(ivarim-1+icpv+i)
                end do
            end if

!
!               --  Prise en compte de la dilatation thermique
!                   En realite, on ne sait pas encore faire ...
!                   La routine coqgth nous arretera s'il y a de la dilatation
            call coqgth(materPara, ep, epsm, deps)
!               -- endommagement plus plasticite
            relaPlas = compor(PLAS_NAME)
            call kit_glrc_dm_vmis(plateCara, plateOrie, &
                                  materPara%jvMaterCode, relaPlas, epsm, deps, ecr, &
                                  option, sigm, sig, ecrp, dsidep, &
                                  zr(jvCarcri), codret)
        else
            call utmess('F', 'PLATE1_7', sk=relaComp)
        end if
!
        if (lNonLine) then
            do i = 1, nbvari
                zr(ivarip-1+icpv+i) = ecrp(i)
            end do
        end if
!
        if (lVect) then
!
!               -- EFFORTS RESULTANTS (N, M ET Q)

            do i = 1, 3
                n(i) = sig(i)
                m(i) = sig(i+3)
            end do
!               -- PRISE EN COMPTE DE L'EXCENTREMENT DANS LE TERME DE FLEXION (VOIR AUSSI DKTNLI)

            if (relaComp(1:4) .ne. 'ELAS' .and. excen .gt. 0.d0) then
                do i = 1, 3
                    m(i) = m(i)+excen*n(i)
                end do
            end if
            do i = 1, 2
                q(i) = sig(i+6)
            end do
!
!               -- CALCUL DE DIV(SIGMA) ET RECOPIE DE N ET M DANS 'PCONTPR'
!                  BTSIG = BTSIG + BFT*M + BMT*N + BCT*Q
!
            do k = 1, 3
                effint((kpg-1)*8+k) = n(k)
                effint((kpg-1)*8+3+k) = m(k)
            end do
!
            if (isPlateQ4GG(plateCara)) then
                do k = 1, 2
                    effint((kpg-1)*8+6+k) = q(k)
                end do
            end if
!
            do ino = 1, nno
                do k = 1, 3
                    vecloc(1, ino) = vecloc(1, ino)+bm(k, 2*(ino-1)+1)*n(k)*poids
                    vecloc(2, ino) = vecloc(2, ino)+bm(k, 2*(ino-1)+2)*n(k)*poids
                    vecloc(3, ino) = vecloc(3, ino)+bf(k, 3*(ino-1)+1)*m(k)*poids
                    vecloc(5, ino) = vecloc(5, ino)+bf(k, 3*(ino-1)+2)*m(k)*poids
                    vecloc(4, ino) = vecloc(4, ino)-bf(k, 3*(ino-1)+3)*m(k)*poids
                end do
!
!                   -- PRISE EN COMPTE DU CISAILLEMENT
                if (isPlateQ4GG(plateCara)) then
                    do k = 1, 2
                        vecloc(3, ino) = vecloc(3, ino)+bc(k, 3*(ino-1)+1)*q(k)*poids
                        vecloc(5, ino) = vecloc(5, ino)+bc(k, 3*(ino-1)+2)*q(k)*poids
                        vecloc(4, ino) = vecloc(4, ino)-bc(k, 3*(ino-1)+3)*q(k)*poids
                    end do
                end if
            end do
        end if
!
        if (lMatr) then
!           -- CALCUL DES MATRICES TANGENTES MATERIELLES (DM,DF,DMF)
            l = 0
            do i = 1, 3
                do j = 1, 3
                    l = l+1
                    dm(l) = dm(l)+poids*dsidep(j, i)
                    dmf(l) = dmf(l)+poids*dsidep(j, i+3)
                    df(l) = df(l)+poids*dsidep(j+3, i+3)
                end do
            end do
!
            if (isPlateQ4GG(plateCara)) then
                dc(1, 1) = poids*dcc(1, 1)
                dc(2, 2) = dc(1, 1)
                dc(1, 2) = 0.d0
                dc(2, 1) = 0.d0
            end if
!
!               CALCUL DE LA MATRICE TANGENTE :
!               KTANG = KTANG + BFT*DF*BF + BMT*DM*BM + BMT*DMF*BF
!                             + BCT*DC*BC
!           -- MEMBRANE
            call utbtab('CUMUL', 3, 2*nno, dm, bm, &
                        work, memb)
!
!           -- FLEXION
            call utbtab('CUMUL', 3, 3*nno, df, bf, &
                        work, flex)
!
!           -- CISAILLEMENT
            if (isPlateQ4GG(plateCara)) then
                call utbtab('CUMUL', 2, 3*nno, dc, bc, &
                            work, flex)
            end if

!           -- COUPLAGE
            call utctab('CUMUL', 3, 3*nno, 2*nno, dmf, &
                        bf, bm, work, mefl)
        end if
    end do
!
    if (lNonLine) then
!           -- PASSAGE DES EFFORTS GENERALISES AUX POINTS D'INTEGRATION
!              DU REPERE INTRINSEQUE AU REPERE LOCAL
!              STOCKAGE DES EFFORTS GENERALISES
!
        call dxefro(npg, plateOrie%t2iu, effint, efforp)
        do i = 1, nbcont*npg
            zr(icontp-1+i) = efforp(i)
        end do
    end if
!
    if (lVect) then
        call utpvlg(nno, 6, pgl, vecloc, zr(ivectu))
    end if
!
    if (lSigm) then
        zi(jcret) = codret
    end if
!
    if (lMatr) then
        if (isPlateTria(plateCara)) then
            call dxtloc(flex, memb, mefl, ctor, matloc)
        else if (isPlateQuad(plateCara)) then
            call dxqloc(flex, memb, mefl, ctor, matloc)
        end if
        call jevech('PMATUUR', 'E', imatuu)
        call utpslg(nno, 6, pgl, matloc, zr(imatuu))
    end if
!
end subroutine
