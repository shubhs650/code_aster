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
subroutine lcpllg(toler, itmax, mod, nbmat, mater, &
                  nr, nvi, deps, sigd, vind, &
                  seuil, cutLevel, sigf, vinf, devg, &
                  devgii, irtet)
!
    implicit none
#include "asterfort/calcpj.h"
#include "asterfort/codent.h"
#include "asterfort/codree.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcopil.h"
#include "asterfort/lglcov.h"
#include "asterfort/lgldom.h"
#include "asterfort/lglini.h"
#include "asterfort/lglite.h"
#include "asterfort/prjsom.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    integer(kind=8) :: itmax, nbmat, nr, nvi, cutLevel, irtet
    real(kind=8) :: toler, mater(nbmat, 2), deps(6), sigd(6)
    real(kind=8) :: vind(*), sigf(6), vinf(*), seuil, devg(6), devgii
    character(len=8) :: mod
! --- BUT : LOI DE COMPORTEMENT PLASTIQUE POUR LA MECANIQUE DES ROCHES -
! ------- : D'APRES LA LOI DE LAIGLE -----------------------------------
! ======================================================================
! IN  : TOLER  : VALEUR DE LA TOLERANCE DE CONVERGENCE -----------------
! --- :        : (RESI_INTE) --------------------------------------
! --- : ITMAX  : NOMBRE D'ITERATIONS MAXIMUM A CONVERGENCE -------------
! --- :        : (ITER_INTE_MAXI) --------------------------------------
! --- : MOD    : TYPE DE MODELISATION ----------------------------------
! --- : NBMAT  : NOMBRE DE PARAMETRES MATERIAU -------------------------
! --- : MATER  : PARAMETRES MATERIAU -----------------------------------
! --- : NR     : NOMBRE DE RELATIONS NON LINEAIRES ---------------------
! --- : NVI    : NOMBRE DE VARIABLES INTERNES --------------------------
! --- : DEPS   : ACCROISSEMENTS DE DEFORMATIONS A L'ITERATION COURANTE -
! --- : SIGD   : CONTRAINTES A L'INSTANT PRECEDENT ---------------------
! --- : VIND   : VARIABLES INTERNES A L'INSTANT PRECEDENT --------------
! --- : SEUIL  : VARIABLE SEUIL ELASTIQUE ------------------------------
! OUT : SIGF   : CONTRAINTES A L'INSTANT COURANT -----------------------
! --- : VINF   : VARIABLES INTERNES A L'INSTANT COURANT ----------------
! --- : DEVG   : DEVIATEUR DU TENSEUR G, DIRECTION D'ECOULEMENT --------
! --- : DEVGII : NORME DU DEVIATEUR DE G -------------------------------
! --- : IRTET  : CONTROLE DU REDECOUPAGE DU PAS DE TEMPS ---------------
! ======================================================================
    integer(kind=8) :: ii, ndt, ndi, iter, irteti, codret
    real(kind=8) :: sige(6), lgleps, gamp, se(6), siie, invare, yd(10)
    real(kind=8) :: gamps, invars, b, s(6), delta, dy(10), yf(10)
    real(kind=8) :: fiter, dkooh(6, 6), epsf(6), i1, traceg, trois
    real(kind=8) :: evp, evps
    character(len=10) :: ctol, citer
    character(len=24) :: valk(2)
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
! --- INITIALISATION DE PARAMETRE --------------------------------------
! ======================================================================
    parameter(trois=3.0d0)
    parameter(lgleps=1.0d-8)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
    call jemarq()
! ======================================================================
! --- INITIALISATION DES VARIABLES -------------------------------------
! ======================================================================
    irteti = 0
    delta = 0.0d0
    gamp = vind(1)
    evp = vind(2)
    sige(1:ndt) = sigf(1:ndt)
    call lcdevi(sige, se)
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    siie = ddot(b_n, se, b_incx, se, b_incy)
    siie = sqrt(siie)
    invare = trace(ndi, sige)
! ======================================================================
! --- INITIALISATION YD = (SIG, INVAR, GAMP, EVP, DELTA) ---------------
! ======================================================================
    yd(1:ndt) = se(1:ndt)
    yd(ndt+1) = invare
    yd(ndt+2) = gamp
    yd(ndt+3) = evp
    yd(ndt+4) = delta
! ======================================================================
! --- CALCUL A PRIORI DE LA PROJECTION AU SOMMET -----------------------
! ======================================================================
    call calcpj(nbmat, mater, gamp, evp, sigd, &
                sige, lgleps, invare, gamps, evps, &
                invars, b)
! ======================================================================
! --- FAUT-IL FAIRE UNE PROJECTION AU SOMMET DU DOMAINE ? --------------
! ======================================================================
    if (prjsom(nbmat, mater, invare, invars, b, siie, 'SUPERIEUR')) then
! ======================================================================
! --- LA PROJECTION AU SOMMET DU DOMAINE EST RETENUE -------------------
! ======================================================================
        do ii = 1, ndt
            sigf(ii) = 0.0d0
        end do
        do ii = 1, ndi
            sigf(ii) = invars/trois
        end do
        call lcopil('ISOTROPE', mod, mater(1, 1), dkooh)
        epsf(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigf(1:ndt))
        if (mod .eq. 'C_PLAN') then
            sigf(3) = 0.0d0
            epsf(3) = dkooh(3, 1)*sigf(1)+dkooh(3, 2)*sigf(2)+dkooh(3, 4)*sigf(4)
        end if
        vinf(1) = gamps
        vinf(2) = evps
        vinf(nvi) = 1.0d0
        irteti = 0
    else
! ======================================================================
! --- LA PROJECTION AU SOMMET DU DOMAINE N'EST PAS RETENUE -------------
! ======================================================================
! --- CALCUL INITIAL (ITERATION 0) -------------------------------------
! ======================================================================
        call lglini(yd, nbmat, mater, seuil, sigd, &
                    deps, devg, devgii, traceg, dy, &
                    codret)
        if (codret .ne. 0) goto 100
        iter = 0
1       continue
! ======================================================================
! --- ITERATION ITER ---------------------------------------------------
! ======================================================================
! --- INCREMENTATION DES VARIABLES -------------------------------------
! ======================================================================
        yf(1:nr-1) = yd(1:nr-1)+dy(1:nr-1)
! ======================================================================
! --- VERIFICATION DE LA COHERENCE DE GAMP -----------------------------
! ======================================================================
        if (yf(ndt+2) .lt. 0.0d0) then
! ======================================================================
! --- GAMP < 0 ---------------------------------------------------------
! --- PEUT-ON FAIRE UN DECOUPAGE DE L'INCREMENT DE DEPLACEMENT ? -------
! ======================================================================
            if (cutLevel .eq. 0 .or. cutLevel .eq. 1) then
                call codent(iter, 'G', citer)
                call codree(toler, 'E', ctol)
                valk(1) = citer
                valk(2) = ctol
                call utmess('I', 'ALGORITH2_57', nk=2, valk=valk)
                irteti = 3
                goto 100
            else
                call utmess('I', 'ALGELINE5_52')
!               CALL UTEXCM(23,'ALGELINE5_52',0,' ',1,VALI,1,VALR)
                codret = 2
            end if
        end if
! ======================================================================
        delta = delta+dy(nr)
        yf(nr) = delta
! ======================================================================
! --- CALCUL DE F A L'ITERATION ITER + 1 -------------------------------
! ======================================================================
        call lgldom(nbmat, mater, yf, fiter)
! ======================================================================
! --- A-T-ON CONVERGE ? ------------------------------------------------
! ======================================================================
        if (lglcov(fiter, toler)) then
! ======================================================================
! --- IL Y A CONVERGENCE -----------------------------------------------
! ======================================================================
! --- MISE A JOUR DES VARIABLES INTERNES -------------------------------
! ======================================================================
            s(1:ndt) = yf(1:ndt)
            i1 = yf(ndt+1)
            gamp = yf(ndt+2)
            evp = yf(ndt+3)
            do ii = 1, ndt
                sigf(ii) = s(ii)
            end do
            do ii = 1, ndi
                sigf(ii) = sigf(ii)+i1/trois
            end do
            call lcopil('ISOTROPE', mod, mater(1, 1), dkooh)
            epsf(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigf(1:ndt))
            if (mod .eq. 'C_PLAN') then
                sigf(3) = 0.0d0
                epsf(3) = dkooh(3, 1)*sigf(1)+dkooh(3, 2)*sigf(2)+dkooh(3, 4)*sigf(4)
            end if
            vinf(1) = gamp
            vinf(2) = evp
            vinf(nvi) = 1.0d0
            irteti = 0
        else
! ======================================================================
! --- IL N'Y A PAS CONVERGENCE -----------------------------------------
! ======================================================================
            if (iter .lt. itmax) then
! ======================================================================
! --- LE NOMBRE D'ITERATION MAXIMAL N'A PAS ETE ATTEINT ----------------
! ======================================================================
                iter = iter+1
! ======================================================================
! --- NOUVEAU CALCUL PLASTIQUE -----------------------------------------
! ======================================================================
                call lglite(yf, nbmat, mater, fiter, devg, &
                            devgii, traceg, dy, codret)
                irteti = 1
                if (codret .ne. 0) goto 100
            else
! ======================================================================
! --- ON NE CONVERGE VRAIMENT PAS ! ------------------------------------
! ======================================================================
! --- FAUT-IL PROJETER AU SOMMET DU DOMAINE ? --------------------------
! ======================================================================
                if (prjsom(nbmat, mater, invare, invars, b, siie, 'INFERIEUR')) then
! ======================================================================
! --- DECOUPAGE
! ======================================================================
                    if (cutLevel .eq. 0 .or. cutLevel .eq. 1) then
                        call codent(iter, 'G', citer)
                        call codree(toler, 'E', ctol)
                        valk(1) = citer
                        valk(2) = ctol
                        call utmess('I', 'ALGORITH2_57', nk=2, valk=valk)
                        irteti = 3
                        goto 100
                    else
                        call utmess('I', 'ALGELINE5_52')
!                     CALL UTEXCM(23,'ALGELINE5_52',0,' ',1,VALI,1,VALR)
                        codret = 2
                    end if
! ======================================================================
! --- ON PROJETE AU SOMMET DU DOMAINE ----------------------------------
! ======================================================================
! --- MISE A JOUR DES VARIABLES INTERNES -------------------------------
! ======================================================================
                    do ii = 1, ndt
                        sigf(ii) = s(ii)
                    end do
                    do ii = 1, ndi
                        sigf(ii) = invars/trois
                    end do
                    call lcopil('ISOTROPE', mod, mater(1, 1), dkooh)
                    epsf(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigf(1:ndt))
                    if (mod .eq. 'C_PLAN') then
                        sigf(3) = 0.0d0
                        epsf(3) = dkooh(3, 1)*sigf(1)+dkooh(3, 2)*sigf(2)+dkooh(3, 4)*sigf(4)
                    end if
                    vinf(1) = gamps
                    vinf(2) = evp
                    vinf(nvi) = 1.0d0
                    irteti = 0
                else
! ======================================================================
! --- IL N'Y A PAS CONVERGENCE -----------------------------------------
! --- PEUT-ON FAIRE UN DECOUPAGE DE L'INCREMENT DE DEPLACEMENT ? -------
! ======================================================================
                    if (cutLevel .eq. 0 .or. cutLevel .eq. 1) then
                        call codent(iter, 'G', citer)
                        call codree(toler, 'E', ctol)
                        valk(1) = citer
                        valk(2) = ctol
                        call utmess('I', 'ALGORITH2_57', nk=2, valk=valk)
                        irteti = 3
                        goto 100
                    else
                        call utmess('I', 'ALGELINE5_52')
!                     CALL UTEXCM(23,'ALGELINE5_52',0,' ',1,VALI,1,VALR)
                        codret = 2
                    end if
                end if
            end if
        end if
        if (irteti .eq. 1) goto 1
    end if
100 continue
    if (irteti .eq. 3) then
        irtet = 1
    else
        irtet = 0
    end if
    if (codret .eq. 2) irtet = 2
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
