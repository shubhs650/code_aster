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

subroutine nmspfm(BEHinteg, typmod, ndim, nno, nddl, nddlsym, &
                  nno_p, nno_s, nddl_p, nddl_s, npg, lgpg, &
                  wref, vff_s, vf_p, pgl, geom, nsigm, mate, matint, matpou, option, &
                  deplm, ddepl, sigm, sigp, fint, &
                  ktan, vim, vip, carcri, compor, &
                  tm, tp, coopg, matsym, lMatr, lVect, lSigm, lElas, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/jevech.h"
#include "asterfort/r8inir.h"
#include "asterfort/assert.h"
#include "asterfort/codere.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmspci.h"
#include "asterfort/ffpoutimo.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpplg.h"
#include "asterfort/utmess.h"
#include "asterfort/lonelesp.h"
#include "asterfort/lcspelas.h"
#include "MultiFiber_type.h"
#include "blas/ddot.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHinteg
    integer(kind=8) :: nddlsym
    integer(kind=8) :: ndim, nno, nddl, nno_p, nno_s, nddl_p, nddl_s, npg
    integer(kind=8) :: lgpg
    integer(kind=8) :: nsigm, mate, codret
    real(kind=8) :: geom(3*nno), wref(npg), vf_p(ndim, npg), vff_s(nno, npg)
    real(kind=8) :: pgl(3, 3), deplm(nddl), ddepl(nddl), tm, tp
    real(kind=8) :: fint(nddl), ktan(nddlsym), coopg(4, npg)
    real(kind=8) :: sigm(nsigm, npg), sigp(nsigm, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    character(len=8), intent(in) :: typmod(2)
    character(len=8), intent(in)   :: matpou, matint
    character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    aster_logical, intent(in) :: matsym, lMatr, lVect, lSigm, lElas
!
!-----------------------------------------------------------------------
!  CALCUL DE FULL_MECA, RIGI_MECA_TANG, RAPH_MECA
!      POUR LES ELEMENTS 3D_INTERF_POU (TE0366)
!-----------------------------------------------------------------------
! IN  NDIM    DIMENSION DU PROBLEME (=3)
! IN  NNO     NOMBRE DE NOEUDS TOTAL DE L'ELEMENT
! IN  NDDL    NOMBRE DE DEGRES DE LIBERTE EN DEPL TOTAL (3 PAR NOEUD)
! IN  NDDLSYM    NOMBRE DE DEGRES DE LIBERTE EN DEPL (MATR SYME)
! IN  NNO_P   NOMBRE DE NOEUDS DE LA PARTIE POUTRE
! IN  NNO_S   NOMBRE DE NOEUDS DE LA PARTIE SOL
! IN  NDDL_P  NOMBRE DE DEGRES DE LIBERTE EN DEPL DE LA POUTRE (6 PAR NOEUD)
! IN  NDDL_S  NOMBRE DE DEGRES DE LIBERTE EN DEPL DU SOL (3 PAR NOEUD)
! IN  NPG    NOMBRE DE POINTS DE GAUSS
! IN  LGPG   NOMBRE DE VARIABLES INTERNES
! IN  WREF   POIDS DE REFERENCE DES POINTS DE GAUSS
! IN  VF_P    COORDONNES DES POINTS DE GAUSS
! IN  VFF_S   VALEUR DES FONCTIONS DE FORME DE L'ELEMENT SOL
! IN  PGL     MATRICE DE PROJECTION DU REPERE GLOBAL --> LOCAL
! IN  GEOM    COORDONNEES DES NOEUDS
! IN  NSIGM   DIMENSION DES CONTRAINTES/FORCES D'INTERACTION
! IN  MATE    MATERIAU CODE
! IN  MATINT  NOM DU MATERIAU DE L'INTERFACE
! IN  MATPOU  NOM DU MATERIAU DE LA POUTRE
! IN  OPTION OPTION DE CALCUL
! IN  DEPLM  DEPLACEMENTS NODAUX AU DEBUT DU PAS DE TEMPS
! IN  DDEPL  INCREMENT DES DEPLACEMENTS NODAUX
! IN  SIGM    CONTR LOCALES AUX POINTS DE GAUSS - (FLX, FLY, FLZ)
! OUT SIGP   CONTR LOCALES AUX POINTS DE GAUSS + (FLX, FLY, FLZ)
! OUT FINT   FORCES NODALES
! OUT KTAN   MATRICE TANGENTE (STOCKEE EN TENANT COMPTE DE LA SYMETRIE)
! IN  VIM    VARIABLES INTERNES AU DEBUT DU PAS DE TEMPS
! OUT VIP    VARIABLES INTERNES A LA FIN DU PAS DE TEMPS
! IN  CRIT   VALEURS DE L'UTILISATEUR POUR LES CRITERES DE CONVERGENCE
! IN  COMPOR NOM DE LA LOI DE COMPORTEMENT
! IN  MATSYM INFORMATION SUR LA MATRICE TANGENTE : SYMETRIQUE OU PAS
! IN  COOPG  COORDONNEES GEOMETRIQUES DES PG + POIDS
! OUT CODRET CODE RETOUR DE L'INTEGRATION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    character(len=4), parameter :: fami = "RIGI"
    integer(kind=8) :: cod(npg)
    integer(kind=8) :: icage
    integer(kind=8) :: kpg, n, ni, mj, kk, p, q, lddl, iddl
    integer(kind=8) :: ddl_sp(3*nno_s+6*nno_p), ddl_s(3*nno_s)
    integer(kind=8) :: ddl_p(6*nno_p)
    integer(kind=8) :: no1, no2
    integer(kind=8) :: niloc, mjloc
    real(kind=8) :: xl
    real(kind=8) :: vff_p(18)
    real(kind=8) :: fintl(nddl), ktanl(nddlsym)
    real(kind=8) :: ur(3), dur(3), dsidep(3, 3), poids(npg)
    real(kind=8) :: b(3, nddl)
    real(kind=8) :: angmas(3)
    real(kind=8) :: sigmo(3), sigma(3)
    real(kind=8) :: sigmrot(nsigm, npg), sigprot(nsigm, npg)
    real(kind=8) :: deplmrot(nddl), ddeplrot(nddl)
!
    real(kind=8) :: hy1, hz1, d1, tsec
    real(kind=8) :: factg(3), factsig(3)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    ur = 0.d0
    dur = 0.d0
    cod = 0
    lddl = nddlsym
    ddl_s = [(iddl, iddl=1, 3*nno_s, 1)]
    ddl_p = (/7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6/)
    ddl_sp = [ddl_s, ddl_p+3*nno_s]

!
! - Get some caracteristics
!

! - cross-section
    call jevech('PCAGEPO', 'L', icage)
    tsec = zr(icage-1+4)
    if (tsec .eq. 1) then
        hy1 = zr(icage-1+1)
        hz1 = zr(icage-1+2)
        factg = (/2*hy1+2*hz1, hy1, hz1/)
    else if (tsec .eq. 2) then
        d1 = 2*zr(icage-1+3)
        factg = (/3.1416*d1, d1, d1/)
    else
        call utmess('F', 'ELEMENTS2_56')
    end if

! - pile length
    no1 = 23
    no2 = 25
    xl = lonelesp(geom, 27, no1, no2)

!
! - Initialize
!
    if (lVect) then
        call r8inir(nddl, 0.d0, fintl, 1)
    end if
    if (lMatr) then
        call r8inir(lddl, 0.d0, ktanl, 1)
    end if

    angmas = r8vide()

! - Calculate Jacobian
    poids = wref*xl/2

!
! - Rotate displacements and internal stresses from global to local coordinates
!
    call utpvgl(nno_s, 3, pgl, deplm, deplmrot)
    call utpvgl(nno_s, 3, pgl, ddepl, ddeplrot)
    call utpvgl(nno_p, 6, pgl, deplm(nddl_s+1), deplmrot(nddl_s+1))
    call utpvgl(nno_p, 6, pgl, ddepl(nddl_s+1), ddeplrot(nddl_s+1))
    if (lSigm) then
        call utpvgl(npg, nsigm, pgl, sigm, sigmrot)
    end if

!
! - Loop on Gauss points
!
    do kpg = 1, npg

! ----- Get matrix B giving relative displacement from node displacements
        call ffpoutimo(vf_p(1, kpg), xl, mate, matpou, vff_p)
        call nmspci(nno_p, nno_s, vff_p, vff_s(1, kpg), b)

! ----- Calculation of relative displacement
        ur = matmul(b, deplmrot(ddl_sp))
        if (lVect) then
            dur = matmul(b, ddeplrot(ddl_sp))
        end if

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHinteg)
        BEHinteg%behavESVA%behavESVAGeom%coorElga(kpg, 1:3) = coopg(1:3, kpg)
        BEHinteg%behavESVA%behavESVAOther%rotpg(1:3*3) = reshape(pgl, (/3*3/))

! ----- Compute behaviour
        sigmo = 0.d0
        do n = 1, 3
            sigmo(n) = sigmrot(n, kpg)/factg(n)
        end do
        sigma = 0.d0
        if (lElas) then
            call lcspelas(fami, kpg, ksp, ndim, &
                          mate, matint, carcri, tm, tp, 3, ur, &
                          dur, 3, sigmo, 1, vim(1, kpg), option, &
                          sigma, vip(1, kpg), 3*3, dsidep, cod(kpg), BEHinteg)
        else
            call nmcomp(BEHinteg, &
                        ndim, option, typmod, &
                        tm, tp, compor, carcri, '                ', &
                        3, ur, dur, 3, sigmo, &
                        vim(1, kpg), &
                        sigma, vip(1, kpg), 3*3, dsidep, cod(kpg))
        end if
        if (cod(kpg) .eq. 1) goto 900

! ----- Stresses
        if (lSigm) then
            do n = 1, 3
                sigprot(n, kpg) = factg(n)*sigma(n)
                factsig(n) = factg(n)*sigma(n)
            end do
        end if

! ----- Internal forces
        if (lVect) then
            ASSERT(lSigm)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            do ni = 1, nddl
                niloc = ddl_sp(ni)
                fintl(ni) = fintl(ni)+poids(kpg)*ddot(b_n, b(1, niloc), b_incx, factsig, b_incy)
            end do
        end if

! ----- Rigidity matrix
        if (lMatr) then
            if (matsym) then
                kk = 0
                do ni = 1, nddl
                    do mj = 1, ni
                        kk = kk+1
                        niloc = ddl_sp(ni)
                        mjloc = ddl_sp(mj)
                        do p = 1, 3
                            do q = 1, 3
                                ktanl(kk) = ktanl(kk)+ &
                                            poids(kpg)*factg(p)*b(p, niloc)*dsidep(p, q)*b(q, mjloc)
                            end do
                        end do
                    end do
                end do
            else
                kk = 0
                do ni = 1, nddl
                    do mj = 1, nddl
                        kk = kk+1
                        niloc = ddl_sp(ni)
                        mjloc = ddl_sp(mj)
                        do p = 1, 3
                            do q = 1, 3
                                ktanl(kk) = ktanl(kk)+ &
                                            poids(kpg)*factg(p)*b(p, niloc)*dsidep(p, q)*b(q, mjloc)
                            end do
                        end do
                    end do
                end do
            end if
        end if
    end do

900 continue

!
! - Local to global coordinate system
!
    if (lMatr) then
        if (matsym) then
            call utpslg(nno_s+2*nno_p, 3, pgl, ktanl, ktan)
        else
            call utpplg(nno_s+2*nno_p, 3, pgl, ktanl, ktan)
        end if
    end if
    if (lVect) then
        call utpvlg(nno_s, 3, pgl, fintl, fint)
        call utpvlg(nno_p, 6, pgl, fintl(nddl_s+1), fint(nddl_s+1))
    end if
    if (lSigm) then
        call utpvlg(npg, nsigm, pgl, sigprot, sigp)
    end if

!
! Return code
!
    call codere(cod, npg, codret)

end subroutine
