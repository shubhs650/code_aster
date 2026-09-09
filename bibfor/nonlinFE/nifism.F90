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
! aslint: disable=W1306,W1504,W0413
!
subroutine nifism(BEHInteg, &
                  ndim, nnod, nnog, nnop, npg, &
                  iw, vffd, vffg, vffp, idffd, &
                  idffg, vu, vg, vp, geomi, &
                  typmod, option, compor, lgpg, &
                  carcri, instm, instp, ddlm, ddld, &
                  sigm, vim, sigp, vip, &
                  lMatr, lVect, lMatrPred, vect, matr, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nirela.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmepsi.h"
#include "asterfort/nmmalu.h"
#include "MeshTypes_type.h"
#include "asterfort/rcvala.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8) :: ndim, nnod, nnog, nnop, npg, iw, idffd, idffg, lgpg
    integer(kind=8) :: vu(3, 27), vg(27), vp(27)
    integer(kind=8) :: codret
    real(kind=8) :: vffd(nnod, npg), vffg(nnog, npg), vffp(nnop, npg)
    real(kind=8) :: instm, instp
    real(kind=8) :: geomi(ndim, nnod), ddlm(*), ddld(*)
    real(kind=8) :: sigm(2*ndim+1, npg), sigp(2*ndim+1, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    real(kind=8) :: vect(*), matr(*)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: compor(COMPOR_SIZE), option
    aster_logical, intent(in) :: lMatr, lVect, lMatrPred
!
! --------------------------------------------------------------------------------------------------
!
!          CALCUL DES FORCES INTERNES POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES GRANDES DEFORMATIONS
!          3D/D_PLAN/AXIS
!
! --------------------------------------------------------------------------------------------------
!
! IN  MATSYM  : MATRICE TANGENTE SYMETRIQUE OU NON
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  nnod    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  nnog    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AU GONFLEMENT
! IN  nnop    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  vffd    : VALEUR  DES FONCTIONS DE FORME LIES AUX DEPLACEMENTS
! IN  vffg    : VALEUR  DES FONCTIONS DE FORME LIES AU GONFLEMENT
! IN  vffp    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
! IN  IDFF1   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  IDFF2   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  VU      : TABLEAU DES INDICES DES DDL DE DEPLACEMENTS
! IN  VG      : TABLEAU DES INDICES DES DDL DE GONFLEMENT
! IN  VP      : TABLEAU DES INDICES DES DDL DE PRESSION
! IN  GEOMI   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTM   : INSTANT PRECEDENT
! IN  INSTP   : INSTANT DE CALCUL
! IN  DDLM    : DEGRES DE LIBERTE A L'INSTANT PRECEDENT
! IN  DDLD    : INCREMENT DES DEGRES DE LIBERTE
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT VECT    : FORCES INTERNES
! OUT MATR    : MATRICE DE RIGIDITE (RIGI_MECA_TANG ET FULL_MECA)
! OUT CODRET  : CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1
    aster_logical, parameter :: grand = ASTER_TRUE
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8), parameter :: ndimBizarre = 3
    aster_logical :: axi, nonloc
    integer(kind=8) :: kpg, nddl, ndu, iret
    integer(kind=8) :: ia, na, ra, sa, ib, nb, rb, sb, ja, jb
    integer(kind=8) :: k2ret(1), lij(3, 3), os, kk
    integer(kind=8) :: viaja
    integer(kind=8) :: cod(27)
    real(kind=8) :: geomm(3*MT_NNOMAX3D), geomp(3*MT_NNOMAX3D)
    real(kind=8) :: deplm(3*MT_NNOMAX3D), depld(3*MT_NNOMAX3D)
    real(kind=8) :: r, w, wm, wp, dffd(nnod, 4), dff2(nnog, 3)
    real(kind=8) :: presm(27), presd(27)
    real(kind=8) :: gonfm(27), gonfd(27)
    real(kind=8) :: sigm_ldc(6)
    real(kind=8) :: gm, gd, gp, pm, pd, pp
    real(kind=8) :: fPrev(3, 3), jm, ftm(3, 3), corm
    real(kind=8) :: fIncr(3, 3), jd, jp, ftd(3, 3), cord
    real(kind=8) :: taup(6), taudv(6), tauhy, tauldc(6)
    real(kind=8) :: dsidep(6, 3, 3)
    real(kind=8) :: d(6, 3, 3), hdv(3, 3), dhy(6), h(3, 3), hhy
    real(kind=8) :: gradgp(3), c(1)
    real(kind=8) :: t1, t2
    real(kind=8) :: am, ap, bp, boa, aa, bb, daa, dbb, dboa, d2boa
    blas_int :: b_incx, b_incy, b_n
    real(kind=8), parameter :: kr(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: id(3, 3) = reshape((/1.d0, 0.d0, 0.d0, &
                                                    0.d0, 1.d0, 0.d0, &
                                                    0.d0, 0.d0, 1.d0/), (/3, 3/))
    integer(kind=8), parameter :: vij(3, 3) = reshape((/1, 4, 5, 4, 2, 6, 5, 6, 3/), (/3, 3/))
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    nddl = nnod*ndim+nnop+nnog
    ndu = ndim
    if (axi) then
        ndu = 3
    end if
    cod = 0

! - Output quantities
    if (lVect) then
        vect(1:nddl) = 0.d0
    end if
    if (lMatr) then
        matr(1:nddl*(nddl+1)/2) = 0.d0
    end if

! - Extract for fields
    do na = 1, nnod
        do ia = 1, ndim
            geomm(ia+ndim*(na-1)) = geomi(ia, na)+ddlm(vu(ia, na))
            geomp(ia+ndim*(na-1)) = geomm(ia+ndim*(na-1))+ddld(vu(ia, na))
            deplm(ia+ndim*(na-1)) = ddlm(vu(ia, na))
            depld(ia+ndim*(na-1)) = ddld(vu(ia, na))
        end do
    end do
    do sa = 1, nnop
        presm(sa) = ddlm(vp(sa))
        presd(sa) = ddld(vp(sa))
    end do
    do ra = 1, nnog
        gonfm(ra) = ddlm(vg(ra))
        gonfd(ra) = ddld(vg(ra))
    end do

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Kinematic - Previous strains
        call dfdmip(ndim, nnod, axi, geomi, kpg, &
                    iw, vffd(1, kpg), idffd, r, w, &
                    dffd)
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dffd, deplm, fPrev)

! ----- Kinematic - Increment of strains
        call dfdmip(ndim, nnod, axi, geomm, kpg, &
                    iw, vffd(1, kpg), idffd, r, wm, &
                    dffd)
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dffd, depld, fIncr)

        call dfdmip(ndim, nnod, axi, geomp, kpg, &
                    iw, vffd(1, kpg), idffd, r, wp, &
                    dffd)
        call nmmalu(nnod, axi, r, vffd(1, kpg), dffd, lij)

! ----- Gradient
        jm = fPrev(1, 1)*(fPrev(2, 2)*fPrev(3, 3)-fPrev(2, 3)*fPrev(3, 2))- &
             fPrev(2, 1)*(fPrev(1, 2)*fPrev(3, 3)-fPrev(1, 3)*fPrev(3, 2))+ &
             fPrev(3, 1)*(fPrev(1, 2)*fPrev(2, 3)-fPrev(1, 3)*fPrev(2, 2))
        jd = fIncr(1, 1)*(fIncr(2, 2)*fIncr(3, 3)-fIncr(2, 3)*fIncr(3, 2))- &
             fIncr(2, 1)*(fIncr(1, 2)*fIncr(3, 3)-fIncr(1, 3)*fIncr(3, 2))+ &
             fIncr(3, 1)*(fIncr(1, 2)*fIncr(2, 3)-fIncr(1, 3)*fIncr(2, 2))
        jp = jm*jd

! ----- LONGUEUR CARACTERISTIQUE -> PARAMETRE C
        c(1) = 0.d0
        call rcvala(BEHInteg%materPara%jvMaterCode, ' ', 'NON_LOCAL', &
                    0, ' ', [0.d0], &
                    1, 'C_GONF', c(1), k2ret(1), &
                    0)
        nonloc = k2ret(1) .eq. 0 .and. c(1) .ne. 0.d0

! ----- Gonflement
        b_n = to_blas_int(nnog)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        gm = ddot(b_n, vffg(1, kpg), b_incx, gonfm, b_incy)
        b_n = to_blas_int(nnog)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        gd = ddot(b_n, vffg(1, kpg), b_incx, gonfd, b_incy)
        gp = gm+gd

! ----- Pressure
        b_n = to_blas_int(nnop)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pm = ddot(b_n, vffp(1, kpg), b_incx, presm, b_incy)
        b_n = to_blas_int(nnop)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pd = ddot(b_n, vffp(1, kpg), b_incx, presd, b_incy)
        pp = pm+pd

! - CALCUL DES FONCTIONS A, B,... DETERMINANT LA RELATION LIANT G ET J
        call nirela(1, jp, gm, gp, am, &
                    ap, bp, boa, aa, bb, &
                    daa, dbb, dboa, d2boa, iret)

! - PERTINENCE DES GRANDEURS
        if (iret .ne. 0) then
            codret = 1
            goto 999
        end if
        if (jd .le. 1.d-2 .or. jd .gt. 1.d2) then
            codret = 1
            goto 999
        end if
        if (abs(ap/am) .gt. 1.d2) then
            codret = 1
            goto 999
        end if

! - CALCUL DU GRADIENT DU GONFLEMENT POUR LA REGULARISATION
        if (nonloc) then
            call dfdmip(ndim, nnog, axi, geomi, kpg, &
                        iw, vffg(1, kpg), idffg, r, w, &
                        dff2)
            do ia = 1, ndim
                b_n = to_blas_int(nnog)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                gradgp(ia) = ddot(b_n, dff2(1, ia), b_incx, gonfm, b_incy)+ &
                             ddot(b_n, dff2(1, ia), b_incx, gonfd, b_incy)
            end do
        end if

! ----- CALCUL DES DEFORMATIONS ENRICHIES
        corm = (am/jm)**(1.d0/3.d0)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fPrev, b_incx, ftm, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, corm, ftm, b_incx)
        cord = (ap/am/jd)**(1.d0/3.d0)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fIncr, b_incx, ftd, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, cord, ftd, b_incx)

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        cod(kpg) = 0
        dsidep = 0.d0
        sigm_ldc = 0.d0
        taup = 0.d0
        do ia = 1, 3
            sigm_ldc(ia) = sigm(ia, kpg)+sigm(2*ndim+1, kpg)
        end do
        do ia = 4, 2*ndim
            sigm_ldc(ia) = sigm(ia, kpg)*rac2
        end do
        taup = 0.d0
        call nmcomp(BEHInteg, &
                    ndimBizarre, option, typmod, &
                    instm, instp, &
                    compor, carcri, multComp, &
                    9, ftm, ftd, &
                    6, sigm_ldc, &
                    vim(1, kpg), &
                    taup, vip(1, kpg), &
                    54, dsidep, cod(kpg))
!
        if (cod(kpg) .eq. 1) then
            codret = 1
            ASSERT(lVect)
            goto 999
        end if

! - SUPPRESSION DES RACINES DE 2
        if (lVect) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1/rac2, taup(4), b_incx)
        end if

! - MATRICE TANGENTE SANS LES RACINES DE 2
        if (lMatr) then
            b_n = to_blas_int(9)
            b_incx = to_blas_int(6)
            call dscal(b_n, 1/rac2, dsidep(4, 1, 1), b_incx)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(6)
            call dscal(b_n, 1/rac2, dsidep(5, 1, 1), b_incx)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(6)
            call dscal(b_n, 1/rac2, dsidep(6, 1, 1), b_incx)
        end if

! - CALCUL DE LA FORCE INTERIEURE ET DES CONTRAINTES DE CAUCHY
        if (lVect) then
! --------- CONTRAINTE HYDROSTATIQUE ET DEVIATEUR
            tauhy = (taup(1)+taup(2)+taup(3))/3.d0
            do ia = 1, 6
                taudv(ia) = taup(ia)-tauhy*kr(ia)
            end do
            do ia = 1, 2*ndim
                sigp(ia, kpg) = (taudv(ia)+pp*bb*kr(ia))/jp
            end do
            sigp(2*ndim+1, kpg) = (tauhy-pp*bb)/jp

! --------- Internal forces
            do na = 1, nnod
                do ia = 1, ndu
                    kk = vu(ia, na)
                    t1 = 0.d0
                    do ja = 1, ndu
                        t2 = taudv(vij(ia, ja))+pp*bb*id(ia, ja)
                        t1 = t1+t2*dffd(na, lij(ia, ja))
                    end do
                    vect(kk) = vect(kk)+w*t1
                end do
            end do
!
! - VECTEUR FINT:G
            t2 = tauhy*aa-pp*dboa
            do ra = 1, nnog
                kk = vg(ra)
                t1 = vffg(ra, kpg)*t2
                vect(kk) = vect(kk)+w*t1
            end do
!
            if (nonloc) then
                do ra = 1, nnog
                    kk = vg(ra)
                    b_n = to_blas_int(ndim)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(nnog)
                    t1 = c(1)*ddot(b_n, gradgp, b_incx, dff2(ra, 1), b_incy)
                    vect(kk) = vect(kk)+w*t1
                end do
            end if
!
! - VECTEUR FINT:P
            t2 = bp-boa
            do sa = 1, nnop
                kk = vp(sa)
                t1 = vffp(sa, kpg)*t2
                vect(kk) = vect(kk)+w*t1
            end do
        end if
!
! - MATRICE TANGENTE
        if (lMatr) then
            if (lMatrPred) then
                do ia = 1, 3
                    taup(ia) = (sigm(ia, kpg)+sigm(2*ndim+1, kpg))*jm
                end do
                do ia = 4, 2*ndim
                    taup(ia) = sigm(ia, kpg)*jm
                end do
            end if
!
! - CALCUL DU TENSEUR DE CONTRAINTE : TRACE ET PARTIE DEVIATORIQUE
            tauhy = (taup(1)+taup(2)+taup(3))/3.d0
            do ia = 1, 6
                tauldc(ia) = taup(ia)+(pp*bb-tauhy)*kr(ia)
            end do
!
! - PROJECTIONS DEVIATORIQUES ET HYDROSTATIQUES DE LA MAT. TANG.
            do ia = 1, 3
                do ja = 1, 3
                    h(ia, ja) = (dsidep(1, ia, ja)+dsidep(2, ia, ja)+dsidep(3, ia, ja))/3.d0
                    do na = 1, 6
                        d(na, ia, ja) = dsidep(na, ia, ja)-kr(na)*h(ia, ja)
                    end do
                end do
            end do
!
            hhy = (h(1, 1)+h(2, 2)+h(3, 3))/3.d0
            do ia = 1, 3
                do ja = 1, 3
                    hdv(ia, ja) = h(ia, ja)-hhy*id(ia, ja)
                end do
            end do
            do ia = 1, 6
                dhy(ia) = (d(ia, 1, 1)+d(ia, 2, 2)+d(ia, 3, 3))/3.d0
            end do
!
            do na = 1, nnod
                do ia = 1, ndu
                    os = (vu(ia, na)-1)*nddl
!
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
                    do nb = 1, nnod
                        do ib = 1, ndu
                            kk = os+vu(ib, nb)
                            t1 = 0.d0
                            do ja = 1, ndu
                                do jb = 1, ndu
                                    viaja = vij(ia, ja)
                                    t2 = d(viaja, ib, jb)-dhy(viaja)*id(ib, jb)
                                    t1 = t1+dffd(na, lij(ia, ja))*t2*dffd(nb, lij(ib, jb))
                                end do
                            end do
!
                            t2 = pp*jp*dbb
                            t1 = t1+dffd(na, lij(ia, ia))*t2*dffd(nb, lij(ib, ib))
!
! - RIGIDITE GEOMETRIQUE
                            do jb = 1, ndu
                                t1 = t1-dffd( &
                                     na, lij(ia, ib))*dffd(nb, lij(ib, jb))*tauldc(vij(ia, jb))
                            end do
                            matr(kk) = matr(kk)+w*t1
                        end do
                    end do
!
!
! - TERME K:UG      KUG(NDIM,NNO1,NNO2)
                    do rb = 1, nnog
                        kk = os+vg(rb)
                        t1 = 0.d0
                        do ja = 1, ndu
                            t2 = dhy(vij(ia, ja))*aa
                            t1 = t1+dffd(na, lij(ia, ja))*t2*vffg(rb, kpg)
                        end do
                        matr(kk) = matr(kk)+w*t1
                    end do
!
! - TERME K:UP      KUP(NDIM,NNO1,NNO3)
                    do sb = 1, nnop
                        kk = os+vp(sb)
                        t1 = dffd(na, lij(ia, ia))*bb*vffp(sb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end do
                end do
            end do
!
            do ra = 1, nnog
                os = (vg(ra)-1)*nddl
!
! - TERME K:GU      KGU(NDIM,NNO2,NNO1)
                do nb = 1, nnod
                    do ib = 1, ndu
                        kk = os+vu(ib, nb)
                        t1 = 0.d0
                        do jb = 1, ndu
                            t1 = t1+vffg(ra, kpg)*aa*hdv(ib, jb)*dffd(nb, lij(ib, jb))
                        end do
                        matr(kk) = matr(kk)+w*t1
                    end do
                end do
!
! - TERME K:GG      KGG(NNO2,NNO2)
                do rb = 1, nnog
                    kk = os+vg(rb)
                    t2 = hhy*aa**2-pp*d2boa+tauhy*daa
                    t1 = vffg(ra, kpg)*t2*vffg(rb, kpg)
                    matr(kk) = matr(kk)+w*t1
                end do
!
                if (nonloc) then
                    do rb = 1, nnog
                        kk = os+vg(rb)
                        b_n = to_blas_int(ndim)
                        b_incx = to_blas_int(nnog)
                        b_incy = to_blas_int(nnog)
                        t1 = c(1)*ddot(b_n, dff2(ra, 1), b_incx, dff2(rb, 1), b_incy)
                        matr(kk) = matr(kk)+w*t1
                    end do
                end if
!
! - TERME K:GP      KGP(NNO2,NNO3)
                do sb = 1, nnop
                    kk = os+vp(sb)
                    t1 = -vffg(ra, kpg)*dboa*vffp(sb, kpg)
                    matr(kk) = matr(kk)+w*t1
                end do
            end do
!
            do sa = 1, nnop
                os = (vp(sa)-1)*nddl
!
! - TERME K:PU      KPU(NDIM,NNO3,NNO1)
                do nb = 1, nnod
                    do ib = 1, ndu
                        kk = os+vu(ib, nb)
                        t1 = vffp(sa, kpg)*bb*dffd(nb, lij(ib, ib))
                        matr(kk) = matr(kk)+w*t1
                    end do
                end do
!
! - TERME K:PG      KPG(NNO3,NNO2)
                do rb = 1, nnog
                    kk = os+vg(rb)
                    t1 = -vffp(sa, kpg)*dboa*vffg(rb, kpg)
                    matr(kk) = matr(kk)+w*t1
                end do
!
! - TERME K:PP = 0.D0      KPP(NNO3,NNO3)
!
            end do
        end if
    end do

! - Return code summary
    call codere(cod, npg, codret)
!
999 continue
end subroutine
