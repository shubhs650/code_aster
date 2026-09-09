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
! aslint: disable=W1306,W1504,W1501
!
subroutine nifilg(BEHInteg, &
                  ndim, nnod, nnog, nnop, npg, &
                  iw, vffd, vffg, vffp, idffd, &
                  vu, vg, vp, geomi, typmod, &
                  option, compor, lgpg, carcri, &
                  instm, instp, ddlm, ddld, &
                  sigm, vim, sigp, vip, lMatr, &
                  lVect, lSigm, lVari, vect, matr, &
                  matsym, codret)
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
#include "asterfort/dsde2d.h"
#include "asterfort/nirela.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmepsi.h"
#include "asterfort/nmmalu.h"
#include "MeshTypes_type.h"
#include "asterfort/poslog.h"
#include "asterfort/prelog.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "MeshTypes_type.h"
#include "blas/dscal.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    aster_logical :: matsym
    integer(kind=8) :: ndim, nnod, nnog, nnop, npg, iw, idffd, lgpg
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
    aster_logical, intent(in) :: lMatr, lVect, lSigm, lVari
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
    aster_logical, parameter :: grand = ASTER_TRUE, lCplan = ASTER_FALSE
    aster_logical :: axi
    aster_logical :: lCorr
    integer(kind=8) :: kpg, nddl, ndu
    integer(kind=8) :: ia, na, ra, sa, ib, nb, rb, sb, ja, jb
    integer(kind=8) :: lij(3, 3), os, kk
    integer(kind=8) :: viaja, vibjb, vuiana, vgra, vpsa, iret
    integer(kind=8) :: cod(npg)
    real(kind=8) :: geomm(3*MT_NNOMAX3D), geomp(3*MT_NNOMAX3D)
    real(kind=8) ::  deplm(3*MT_NNOMAX3D), deplp(3*MT_NNOMAX3D)
    real(kind=8) :: r, w, wp, dffd(nnod, 4)
    real(kind=8) :: presm(27), presd(27)
    real(kind=8) :: gonfm(27), gonfd(27)
    real(kind=8) :: sigm_ldc(2*ndim), sigp_ldc(2*ndim)
    real(kind=8) :: gm, gd, gp, pm, pd, pp
    real(kind=8) :: fPrev(3, 3), jm, ftm(3, 3), corm, epslPrev(6)
    real(kind=8) :: fCurr(3, 3), jp, ftp(3, 3), corp, epslIncr(6)
    real(kind=8) :: gn(3, 3), lamb(3), logl(3)
    real(kind=8) :: tlogPrev(6), tlogCurr(6), dtde(6, 6)
    real(kind=8) :: pk2Curr(6), pk2Prev(6)
    real(kind=8) :: taup(6), taudv(6), tauhy, tauldc(6)
    real(kind=8) :: dsidep(6, 6)
    real(kind=8) :: d(6, 6), ddev(6, 6), devd(6, 6), dddev(6, 6)
    real(kind=8) :: iddid, devdi(6), iddev(6)
    real(kind=8) :: ftr(3, 3), t1, t2
    real(kind=8) :: am, ap, bp, boa, aa, bb, daa, dbb, dboa, d2boa
    real(kind=8), parameter :: kr(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: id(3, 3) = reshape((/1.d0, 0.d0, 0.d0, &
                                                    0.d0, 1.d0, 0.d0, &
                                                    0.d0, 0.d0, 1.d0/), (/3, 3/))
    integer(kind=8), parameter :: vij(3, 3) = reshape((/1, 4, 5, 4, 2, 6, 5, 6, 3/), (/3, 3/))
    real(kind=8), parameter :: idev(6, 6) = reshape((/2.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, &
                                                      -1.d0, 2.d0, -1.d0, 0.d0, 0.d0, 0.d0, &
                                                      -1.d0, -1.d0, 2.d0, 0.d0, 0.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 3.d0/), &
                                                    (/6, 6/))
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    lCorr = L_CORR(option)
    axi = typmod(1) .eq. 'AXIS'
    nddl = nnod*ndim+nnop+nnog
    ndu = ndim
    if (axi) then
        ndu = 3
    end if
    dsidep = 0.d0
    codret = 0
    cod = 0

! - Output quantities
    if (lVect) then
        vect(1:nddl) = 0.d0
    end if
    if (lMatr) then
        if (matsym) then
            matr(1:nddl*(nddl+1)/2) = 0.d0
        else
            matr(1:nddl*nddl) = 0.d0
        end if
    end if

! - Extract for fields
    do na = 1, nnod
        do ia = 1, ndim
            geomm(ia+ndim*(na-1)) = geomi(ia, na)+ddlm(vu(ia, na))
            geomp(ia+ndim*(na-1)) = geomm(ia+ndim*(na-1))+ddld(vu(ia, na))
            deplm(ia+ndim*(na-1)) = ddlm(vu(ia, na))
            deplp(ia+ndim*(na-1)) = ddlm(vu(ia, na))+ddld(vu(ia, na))
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

! ----- Kinematic - Current strains
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dffd, deplp, fCurr)
        call dfdmip(ndim, nnod, axi, geomp, kpg, &
                    iw, vffd(1, kpg), idffd, r, wp, &
                    dffd)
        call nmmalu(nnod, axi, r, vffd(1, kpg), dffd, lij)

! ----- Gradient
        jm = fPrev(1, 1)*(fPrev(2, 2)*fPrev(3, 3)-fPrev(2, 3)*fPrev(3, 2))- &
             fPrev(2, 1)*(fPrev(1, 2)*fPrev(3, 3)-fPrev(1, 3)*fPrev(3, 2))+ &
             fPrev(3, 1)*(fPrev(1, 2)*fPrev(2, 3)-fPrev(1, 3)*fPrev(2, 2))
        jp = fCurr(1, 1)*(fCurr(2, 2)*fCurr(3, 3)-fCurr(2, 3)*fCurr(3, 2))- &
             fCurr(2, 1)*(fCurr(1, 2)*fCurr(3, 3)-fCurr(1, 3)*fCurr(3, 2))+ &
             fCurr(3, 1)*(fCurr(1, 2)*fCurr(2, 3)-fCurr(1, 3)*fCurr(2, 2))
        if (jp .le. 0.d0) then
            cod(kpg) = 1
            goto 999
        end if

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

! ----- CALCUL DES FONCTIONS A, B,... DETERMINANT LA RELATION LIANT G ET J
        call nirela(2, jp, gm, gp, am, &
                    ap, bp, boa, aa, bb, &
                    daa, dbb, dboa, d2boa, iret)
        if (iret .ne. 0) then
            cod(kpg) = 1
            goto 999
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
        corp = (ap/jp)**(1.d0/3.d0)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fCurr, b_incx, ftp, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, corp, ftp, b_incx)

! ----- Pre-treatment of kinematic quantities
        call prelog(ndim, lgpg, vim(1, kpg), gn, lamb, &
                    logl, ftm, ftp, epslPrev, epslIncr, &
                    tlogPrev, lCorr, cod(kpg))
        if (cod(kpg) .ne. 0) then
            goto 999
        end if

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        cod(kpg) = 0
        dtde = 0.d0
        tlogCurr = 0.d0
        taup = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instm, instp, &
                    compor, carcri, multComp, &
                    6, epslPrev, epslIncr, &
                    6, tlogPrev, &
                    vim(1, kpg), &
                    tlogCurr, vip(1, kpg), &
                    36, dtde, cod(kpg))
        if (cod(kpg) .eq. 1) then
            goto 999
        end if
        do ia = 1, 3
            sigm_ldc(ia) = sigm(ia, kpg)+sigm(2*ndim+1, kpg)
        end do
        do ia = 4, 2*ndim
            sigm_ldc(ia) = sigm(ia, kpg)
        end do

! ----- Post-treatment of sthenic quantities
        call poslog(BEHInteg, &
                    lCorr, lMatr, lSigm, lVari, &
                    tlogPrev, tlogCurr, ftm, &
                    lgpg, vip(1, kpg), ndim, &
                    ftp, dtde, sigm_ldc, lCplan, &
                    instp, gn, &
                    lamb, logl, sigp_ldc, dsidep, pk2Prev, &
                    pk2Curr, iret)
        if (iret .eq. 1) then
            cod(kpg) = 1
            goto 999
        end if

! ----- CONTRAINTE HYDROSTATIQUE ET DEVIATEUR
        b_n = to_blas_int(2*ndim)
        b_incx = to_blas_int(1)
        call dscal(b_n, exp(gp), sigp_ldc, b_incx)
        b_n = to_blas_int(2*ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigp_ldc, b_incx, taup, b_incy)
        tauhy = (taup(1)+taup(2)+taup(3))/3.d0
        do ia = 1, 6
            taudv(ia) = taup(ia)-tauhy*kr(ia)
        end do

! ----- Cauchy stresses
        if (lSigm) then
            do ia = 1, 2*ndim
                sigp(ia, kpg) = (taudv(ia)+pp*bb*kr(ia))/jp
            end do
            sigp(2*ndim+1, kpg) = (tauhy-pp*bb)/jp
        end if

! ----- Internal forces
        if (lVect) then
            ASSERT(lSigm)
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
            t2 = tauhy*aa-pp*dboa
            do ra = 1, nnog
                kk = vg(ra)
                t1 = vffg(ra, kpg)*t2
                vect(kk) = vect(kk)+w*t1
            end do
            t2 = bp-boa
            do sa = 1, nnop
                kk = vp(sa)
                t1 = vffp(sa, kpg)*t2
                vect(kk) = vect(kk)+w*t1
            end do
        end if

! ----- Rigidity matrix
        if (lMatr) then
! Contraintes generalisees EF (bloc mecanique pour la rigidite geometrique)
            if (lCorr) then
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, ftp, b_incx, ftr, b_incy)
            else
                b_n = to_blas_int(2*ndim)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, sigm_ldc, b_incx, taup, b_incy)
                b_n = to_blas_int(2*ndim)
                b_incx = to_blas_int(1)
                call dscal(b_n, jm, taup, b_incx)
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, ftm, b_incx, ftr, b_incy)
            end if
!
! - CALCUL DE L'OPERATEUR TANGENT SYMÉTRISÉ D
            call dsde2d(3, ftr, dsidep, d)
            devd = matmul(idev/3.d0, d)
            ddev = matmul(d, idev/3.d0)
            dddev = matmul(devd, idev/3.d0)
!
! - CALCUL DU TENSEUR DE CONTRAINTE : TRACE ET PARTIE DEVIATORIQUE
            tauhy = (taup(1)+taup(2)+taup(3))/3.d0
!
! - CALCUL DE D^DEV:ID ET ID:D^DEV ET ID:D:ID
            iddid = 0.d0
            do ia = 1, 6
                devdi(ia) = devd(ia, 1)+devd(ia, 2)+devd(ia, 3)
                iddev(ia) = ddev(1, ia)+ddev(2, ia)+ddev(3, ia)
                taudv(ia) = taup(ia)-tauhy*kr(ia)
                tauldc(ia) = taup(ia)+(pp*bb-tauhy)*kr(ia)
                do ja = 1, 3
                    iddid = iddid+kr(ia)*d(ia, ja)
                end do
            end do
!
            if (matsym) then
! - MATRICE SYMETRIQUE
! - TERME K:UX
                do na = 1, nnod
                    do ia = 1, ndu
                        vuiana = vu(ia, na)
                        os = (vuiana-1)*vuiana/2
!
! - TERME K:UU      KUU(NDIM,nnod,NDIM,nnod)
                        do nb = 1, nnod
                            do ib = 1, ndu
                                if (vu(ib, nb) .le. vuiana) then
                                    kk = os+vu(ib, nb)
                                    t1 = 0.d0
! - RIGIDITE DE COMPORTEMENT
                                    do ja = 1, ndu
                                        viaja = vij(ia, ja)
                                        do jb = 1, ndu
                                            vibjb = vij(ib, jb)
                                            t2 = dddev(viaja, vibjb)
                                            t2 = t2+taup(vij(ia, jb))*kr(vij(ib, ja))
                                            t2 = t2+taup(vij(jb, ja))*kr(vij(ia, ib))
                                            t2 = t2-2.d0/3.d0*( &
                                                 taup(viaja)*kr(vibjb)+taup(vibjb)*kr(viaja))
                                            t2 = t2+2.d0/3.d0*tauhy*kr(viaja)*kr(vibjb)
                                            t1 = t1+dffd(na, lij(ia, ja))*t2*dffd(nb, lij(ib, jb) &
                                                                                  )
                                        end do
                                    end do
                                    t2 = pp*jp*dbb
                                    t1 = t1+dffd(na, lij(ia, ia))*t2*dffd(nb, lij(ib, ib))
!
! - RIGIDITE GEOMETRIQUE
                                    do jb = 1, ndu
                                        t1 = t1-dffd( &
                                             na, lij(ia, ib))*dffd(nb, &
                                                                   lij(ib, jb))*tauldc(vij(ia, jb) &
                                                                                       )
                                    end do
                                    matr(kk) = matr(kk)+w*t1
                                end if
                            end do
                        end do
!
! - TERME K:UG      KUG(NDIM,nnod,nnog)
                        t1 = 0.d0
                        do ja = 1, ndu
                            viaja = vij(ia, ja)
                            t2 = (devdi(viaja)+2.d0*taudv(viaja))
                            t1 = t1+dffd(na, lij(ia, ja))*t2
                        end do
                        t1 = t1*aa/3.d0
!
                        do rb = 1, nnog
                            if (vg(rb) .lt. vuiana) then
                                kk = os+vg(rb)
                                matr(kk) = matr(kk)+w*t1*vffg(rb, kpg)
                            end if
                        end do
!
! - TERME K:UP      KUP(NDIM,nnod,nnop)
                        do sb = 1, nnop
                            if (vp(sb) .lt. vuiana) then
                                kk = os+vp(sb)
                                t1 = dffd(na, lij(ia, ia))*bb*vffp(sb, kpg)
                                matr(kk) = matr(kk)+w*t1
                            end if
                        end do
                    end do
                end do
!
! - TERME K:GX
                do ra = 1, nnog
                    vgra = vg(ra)
                    os = (vgra-1)*vgra/2
!
! - TERME K:GU      KGU(NDIM,nnog,nnod)
                    do nb = 1, nnod
                        do ib = 1, ndu
                            if (vu(ib, nb) .lt. vgra) then
                                kk = os+vu(ib, nb)
                                t1 = 0.d0
                                do jb = 1, ndu
                                    vibjb = vij(ib, jb)
                                    t2 = (iddev(vibjb)+2.d0*taudv(vibjb))
                                    t1 = t1+t2*dffd(nb, lij(ib, jb))
                                end do
                                matr(kk) = matr(kk)+w*t1*aa*vffg(ra, kpg)/3.d0
                            end if
                        end do
                    end do
!
! - TERME K:GG      KGG(nnog,nnog)
                    t2 = (iddid/9.d0+2.d0*tauhy/3.d0)*aa**2-pp*d2boa+tauhy*daa
                    do rb = 1, nnog
                        if (vg(rb) .le. vgra) then
                            kk = os+vg(rb)
                            t1 = vffg(ra, kpg)*t2*vffg(rb, kpg)
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
!
! - TERME K:GP      KGP(nnog,nnop)
                    do sb = 1, nnop
                        if (vp(sb) .lt. vgra) then
                            kk = os+vp(sb)
                            t1 = -vffg(ra, kpg)*dboa*vffp(sb, kpg)
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
                end do
!
! - TERME K:PX
                do sa = 1, nnop
                    vpsa = vp(sa)
                    os = (vpsa-1)*vpsa/2
!
! - TERME K:PU      KPU(NDIM,nnop,nnod)
                    do nb = 1, nnod
                        do ib = 1, ndu
                            if (vu(ib, nb) .lt. vpsa) then
                                kk = os+vu(ib, nb)
                                t1 = vffp(sa, kpg)*bb*dffd(nb, lij(ib, ib))
                                matr(kk) = matr(kk)+w*t1
                            end if
                        end do
                    end do
!
! - TERME K:PG      KPG(nnop,nnog)
                    do rb = 1, nnog
                        if (vg(rb) .lt. vpsa) then
                            kk = os+vg(rb)
                            t1 = -vffp(sa, kpg)*dboa*vffg(rb, kpg)
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
!
! - TERME K:PP = 0.D0      KPP(nnop,nnop)
                end do
!
            else
! - MATRICE NON SYMETRIQUE
! - TERME K:UX
                do na = 1, nnod
                    do ia = 1, ndu
                        os = (vu(ia, na)-1)*nddl
!
! - TERME K:UU      KUU(NDIM,nnod,NDIM,nnod)
                        do nb = 1, nnod
                            do ib = 1, ndu
                                kk = os+vu(ib, nb)
                                t1 = 0.d0
! - RIGIDITE DE COMPORTEMENT
                                do ja = 1, ndu
                                    viaja = vij(ia, ja)
                                    do jb = 1, ndu
                                        vibjb = vij(ib, jb)
                                        t2 = dddev(viaja, vibjb)
                                        t2 = t2+taup(vij(ia, jb))*kr(vij(ib, ja))
                                        t2 = t2+taup(vij(jb, ja))*kr(vij(ia, ib))
                                        t2 = t2-2.d0/3.d0*( &
                                             taup(viaja)*kr(vibjb)+kr(viaja)*taup(vibjb))
                                        t2 = t2+2.d0*kr(viaja)*kr(vibjb)*tauhy/3.d0
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
                                         na, lij(ia, ib))*dffd(nb, &
                                                               lij(ib, jb))*tauldc(vij(ia, jb) &
                                                                                   )
                                end do
                                matr(kk) = matr(kk)+w*t1
                            end do
                        end do
!
! - TERME K:UG      KUG(NDIM,nnod,nnog)
                        t1 = 0.d0
                        do ja = 1, ndu
                            viaja = vij(ia, ja)
                            t2 = (devdi(viaja)+2.d0*taudv(viaja))
                            t1 = t1+dffd(na, lij(ia, ja))*t2
                        end do
                        t1 = t1*aa/3.d0
!
                        do rb = 1, nnog
                            kk = os+vg(rb)
                            matr(kk) = matr(kk)+w*t1*vffg(rb, kpg)
                        end do
!
! - TERME K:UP      KUP(NDIM,nnod,nnop)
                        do sb = 1, nnop
                            kk = os+vp(sb)
                            t1 = dffd(na, lij(ia, ia))*bb*vffp(sb, kpg)
                            matr(kk) = matr(kk)+w*t1
                        end do
                    end do
                end do
!
! - TERME K:GX
                do ra = 1, nnog
                    os = (vg(ra)-1)*nddl
!
! - TERME K:GU      KGU(NDIM,nnog,nnod)
                    do nb = 1, nnod
                        do ib = 1, ndu
                            kk = os+vu(ib, nb)
                            t1 = 0.d0
                            do jb = 1, ndu
                                vibjb = vij(ib, jb)
                                t2 = (iddev(vibjb)+2.d0*taudv(vibjb))
                                t1 = t1+t2*dffd(nb, lij(ib, jb))
                            end do
                            matr(kk) = matr(kk)+w*t1*aa*vffg(ra, kpg)/3.d0
                        end do
                    end do
!
! - TERME K:GG      KGG(nnog,nnog)
                    t2 = (iddid/9.d0+2.d0*tauhy/3.d0)*aa**2-pp*d2boa+tauhy*daa
                    do rb = 1, nnog
                        kk = os+vg(rb)
                        t1 = vffg(ra, kpg)*t2*vffg(rb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end do
!
! - TERME K:GP      KGP(nnog,nnop)
                    do sb = 1, nnop
                        kk = os+vp(sb)
                        t1 = -vffg(ra, kpg)*dboa*vffp(sb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end do
                end do
!
! - TERME K:PX
                do sa = 1, nnop
                    os = (vp(sa)-1)*nddl
!
! - TERME K:PU      KPU(NDIM,nnop,nnod)
                    do nb = 1, nnod
                        do ib = 1, ndu
                            kk = os+vu(ib, nb)
                            t1 = vffp(sa, kpg)*bb*dffd(nb, lij(ib, ib))
                            matr(kk) = matr(kk)+w*t1
                        end do
                    end do
!
! - TERME K:PG      KPG(nnop,nnog)
                    do rb = 1, nnog
                        kk = os+vg(rb)
                        t1 = -vffp(sa, kpg)*dboa*vffg(rb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end do
!
! - TERME K:PP = 0.D0      KPP(nnop,nnop)
!
                end do
            end if
        end if
    end do

! - Return code summary
    call codere(cod, npg, codret)
!
999 continue
end subroutine
