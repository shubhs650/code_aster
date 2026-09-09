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
! aslint: disable=W1306,W1504
!
subroutine nifipd(BEHInteg, &
                  ndim, nnod, nnog, nnop, npg, &
                  iw, vffd, vffg, vffp, idffd, &
                  vu, vg, vp, geomi, typmod, &
                  option, compor, lgpg, carcri, &
                  instm, instp, ddlm, ddld, &
                  sigm, vim, sigp, vip, &
                  lMatr, lVect, vect, matr, &
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
#include "asterfort/nmcomp.h"
#include "MeshTypes_type.h"
#include "asterfort/nmepsi.h"
#include "blas/ddot.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
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
    aster_logical, intent(in) :: lMatr, lVect
!
! --------------------------------------------------------------------------------------------------
!
!          CALCUL DES FORCES INTERNES POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES PETITES DEFORMATIONS
!          3D/D_PLAN/AXIS
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AU GONFLEMENT
! IN  NNO3    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  VFF1    : VALEUR  DES FONCTIONS DE FORME LIES AUX DEPLACEMENTS
! IN  VFF2    : VALEUR  DES FONCTIONS DE FORME LIES AU GONFLEMENT
! IN  VFF3    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
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
    aster_logical, parameter :: grand = ASTER_FALSE
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    aster_logical :: axi
    integer(kind=8) :: kpg, nddl
    integer(kind=8) :: ia, na, ra, sa, ib, nb, rb, sb, ja, jb
    integer(kind=8) :: os, kk
    integer(kind=8) :: vuiana, vgra, vpsa
    integer(kind=8) :: cod(27)
    real(kind=8) :: deplm(3*MT_NNOMAX3D), depld(3*MT_NNOMAX3D)
    real(kind=8) :: r, w, dff1(nnod, ndim)
    real(kind=8) :: presm(27), presd(27)
    real(kind=8) :: gonfm(27), gonfd(27)
    real(kind=8) :: gm, gd, pm, pd
    real(kind=8) :: fm(3, 3), epsm(6), deps(6)
    real(kind=8) :: sigma(6), sigmam(6), sigtr
    real(kind=8) :: dsidep(6, 6)
    real(kind=8) :: def(6, nnod, ndim), deftr(nnod, ndim), ddivu, divum
    real(kind=8) :: ddev(6, 6), devd(6, 6), dddev(6, 6)
    real(kind=8) :: iddid, devdi(6), iddev(6)
    real(kind=8) :: t1, t2
    blas_int :: b_incx, b_incy, b_n
    real(kind=8), parameter :: kr(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: idev(6, 6) = reshape((/2.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, &
                                                      -1.d0, 2.d0, -1.d0, 0.d0, 0.d0, 0.d0, &
                                                      -1.d0, -1.d0, 2.d0, 0.d0, 0.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 3.d0/), &
                                                    (/6, 6/))
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    nddl = nnod*ndim+nnop+nnog
    dsidep = 0.d0
    codret = 0
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
        epsm = 0.d0
        call dfdmip(ndim, nnod, axi, geomi, kpg, &
                    iw, vffd(1, kpg), idffd, r, w, &
                    dff1)
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dff1, deplm, fm, epsm)

! ----- Kinematic - Increment of strains
        deps = 0.d0
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dff1, depld, fm, deps)

! ----- Gonflement
        b_n = to_blas_int(nnog)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        gm = ddot(b_n, vffg(1, kpg), b_incx, gonfm, b_incy)
        b_n = to_blas_int(nnog)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        gd = ddot(b_n, vffg(1, kpg), b_incx, gonfd, b_incy)

! ----- Pressure
        b_n = to_blas_int(nnop)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pm = ddot(b_n, vffp(1, kpg), b_incx, presm, b_incy)
        b_n = to_blas_int(nnop)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pd = ddot(b_n, vffp(1, kpg), b_incx, presd, b_incy)

! - CALCUL DES ELEMENTS GEOMETRIQUES
        divum = epsm(1)+epsm(2)+epsm(3)
        ddivu = deps(1)+deps(2)+deps(3)

! - CALCUL DE LA MATRICE B EPS_ij=B_ijkl U_kl
! - DEF (XX,YY,ZZ,2/RAC(2)XY,2/RAC(2)XZ,2/RAC(2)YZ)
        if (ndim .eq. 2) then
            do na = 1, nnod
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = 0.d0
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                end do
            end do
!
! - TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                do na = 1, nnod
                    def(3, na, 1) = fm(3, 3)*vffd(na, kpg)/r
                end do
            end if
        else
            do na = 1, nnod
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = fm(ia, 3)*dff1(na, 3)
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                    def(5, na, ia) = (fm(ia, 1)*dff1(na, 3)+fm(ia, 3)*dff1(na, 1))/rac2
                    def(6, na, ia) = (fm(ia, 2)*dff1(na, 3)+fm(ia, 3)*dff1(na, 2))/rac2
                end do
            end do
        end if
!
! - CALCUL DE TRACE(B)
        do na = 1, nnod
            do ia = 1, ndim
                deftr(na, ia) = def(1, na, ia)+def(2, na, ia)+def(3, na, ia)
            end do
        end do
!
! - DEFORMATION POUR LA LOI DE COMPORTEMENT
        do ia = 1, 3
            epsm(ia) = epsm(ia)+(gm-divum)/3.d0
            deps(ia) = deps(ia)+(gd-ddivu)/3.d0
        end do
!
! - CONTRAINTE EN T- POUR LA LOI DE COMPORTEMENT
        do ia = 1, 3
            sigmam(ia) = sigm(ia, kpg)+sigm(2*ndim+1, kpg)
        end do
        do ia = 4, 2*ndim
            sigmam(ia) = sigm(ia, kpg)*rac2
        end do

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- Integrator
        sigma = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instm, instp, &
                    compor, carcri, multComp, &
                    6, epsm, deps, &
                    6, sigmam, &
                    vim(1, kpg), &
                    sigma, vip(1, kpg), &
                    36, dsidep, cod(kpg))
        if (cod(kpg) .eq. 1) then
            codret = 1
            ASSERT(lVect)
            goto 999
        end if

! - CALCUL DE LA FORCE INTERIEURE ET DES CONTRAINTES DE CAUCHY
        if (lVect) then
! - CONTRAINTES A L'EQUILIBRE
            sigtr = sigma(1)+sigma(2)+sigma(3)
            do ia = 1, 3
                sigma(ia) = sigma(ia)-sigtr/3.d0+(pm+pd)
            end do
!
! - VECTEUR FINT:U
            do na = 1, nnod
                do ia = 1, ndim
                    kk = vu(ia, na)
                    b_n = to_blas_int(2*ndim)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    t1 = ddot(b_n, sigma, b_incx, def(1, na, ia), b_incy)
                    vect(kk) = vect(kk)+w*t1
                end do
            end do
!
! - VECTEUR FINT:G
            t2 = (sigtr/3.d0-pm-pd)
            do ra = 1, nnog
                kk = vg(ra)
                t1 = vffg(ra, kpg)*t2
                vect(kk) = vect(kk)+w*t1
            end do
!
! - VECTEUR FINT:P
            t2 = (divum+ddivu-gm-gd)
            do sa = 1, nnop
                kk = vp(sa)
                t1 = vffp(sa, kpg)*t2
                vect(kk) = vect(kk)+w*t1
            end do
!
! - STOCKAGE DES CONTRAINTES
            do ia = 1, 3
                sigp(ia, kpg) = sigma(ia)
            end do
            do ia = 4, 2*ndim
                sigp(ia, kpg) = sigma(ia)/rac2
            end do
            sigp(2*ndim+1, kpg) = sigtr/3.d0-pm-pd
        end if
!
! - MATRICE TANGENTE
        if (lMatr) then
!
            devd = matmul(idev/3.d0, dsidep)
            ddev = matmul(dsidep, idev/3.d0)
            dddev = matmul(devd, idev/3.d0)
!
! - CALCUL DE D^DEV:ID ET ID:D^DEV ET ID:D:ID/9.D0
            iddid = 0.d0
            do ia = 1, 6
                devdi(ia) = devd(ia, 1)+devd(ia, 2)+devd(ia, 3)
                iddev(ia) = ddev(1, ia)+ddev(2, ia)+ddev(3, ia)
                do ja = 1, 3
                    iddid = iddid+kr(ia)*dsidep(ia, ja)
                end do
            end do
            iddid = iddid/9.d0
!
! - MATRICE SYMETRIQUE
! - TERME K:UX
            do na = 1, nnod
                do ia = 1, ndim
                    vuiana = vu(ia, na)
                    os = (vuiana-1)*vuiana/2
!
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
                    do nb = 1, nnod
                        do ib = 1, ndim
                            if (vu(ib, nb) .le. vuiana) then
                                kk = os+vu(ib, nb)
                                t1 = 0.d0
                                do ja = 1, 2*ndim
                                    do jb = 1, 2*ndim
                                        t1 = t1+def(ja, na, ia)*dddev(ja, jb)*def(jb, nb, ib)
                                    end do
                                end do
                                matr(kk) = matr(kk)+w*t1
                            end if
                        end do
                    end do
!
! - TERME K:UG      KUG(NDIM,NNO1,NNO2)
                    t1 = 0.d0
                    do ja = 1, 2*ndim
                        t1 = t1+def(ja, na, ia)*devdi(ja)
                    end do
                    t1 = t1/3.d0
!
                    do rb = 1, nnog
                        if (vg(rb) .lt. vuiana) then
                            kk = os+vg(rb)
                            matr(kk) = matr(kk)+w*t1*vffg(rb, kpg)
                        end if
                    end do
!
! - TERME K:UP      KUP(NDIM,NNO1,NNO3)
                    do sb = 1, nnop
                        if (vp(sb) .lt. vuiana) then
                            kk = os+vp(sb)
                            t1 = deftr(na, ia)*vffp(sb, kpg)
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
! - TERME K:GU      KGU(NDIM,NNO2,NNO1)
                do nb = 1, nnod
                    do ib = 1, ndim
                        if (vu(ib, nb) .lt. vgra) then
                            kk = os+vu(ib, nb)
                            t1 = 0.d0
                            do jb = 1, 2*ndim
                                t1 = t1+iddev(jb)*def(jb, nb, ib)
                            end do
                            matr(kk) = matr(kk)+w*t1*vffg(ra, kpg)/3.d0
                        end if
                    end do
                end do
!
! - TERME K:GG      KGG(NNO2,NNO2)
                do rb = 1, nnog
                    if (vg(rb) .le. vgra) then
                        kk = os+vg(rb)
                        t1 = vffg(ra, kpg)*iddid*vffg(rb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end if
                end do
!
! - TERME K:GP      KGP(NNO2,NNO3)
                do sb = 1, nnop
                    if (vp(sb) .lt. vgra) then
                        kk = os+vp(sb)
                        t1 = -vffg(ra, kpg)*vffp(sb, kpg)
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
! - TERME K:PU      KPU(NDIM,NNO3,NNO1)
                do nb = 1, nnod
                    do ib = 1, ndim
                        if (vu(ib, nb) .lt. vpsa) then
                            kk = os+vu(ib, nb)
                            t1 = vffp(sa, kpg)*deftr(nb, ib)
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
                end do
!
! - TERME K:PG      KPG(NNO3,NNO2)
                do rb = 1, nnog
                    if (vg(rb) .lt. vpsa) then
                        kk = os+vg(rb)
                        t1 = -vffp(sa, kpg)*vffg(rb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end if
                end do
!
! - TERME K:PP = 0.D0      KPP(NNO3,NNO3)
            end do
        end if
    end do

! - Return code summary
    call codere(cod, npg, codret)
!
999 continue
end subroutine
