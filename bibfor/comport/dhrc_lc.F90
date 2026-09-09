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
! aslint: disable=W1504
!
subroutine dhrc_lc(plateOrie, &
                   epsm, deps, vim, option, &
                   sig, vip, a0, c0, aa_t, &
                   ga_t, ab, gb, ac, gc, &
                   aa_c, ga_c, cstseu, crit, codret, &
                   dsidep, debug)
!
    use plate_type
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/dhrc_calc_a.h"
#include "asterfort/dhrc_calc_b.h"
#include "asterfort/dhrc_calc_c.h"
#include "asterfort/dhrc_calc_g.h"
#include "asterfort/dhrc_calc_n.h"
#include "asterfort/dhrc_jacob.h"
#include "asterfort/dhrc_mat_tan.h"
#include "asterfort/dhrc_sig.h"
#include "asterfort/dxefro.h"
#include "asterfort/jevech.h"
#include "asterfort/mgauss.h"
#include "asterfort/r8inir.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "blas/dgeev.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: epsm(6), deps(6), vim(*), crit(*), cstseu(6)
    real(kind=8), intent(in) :: a0(6, 6), c0(2, 2, 2)
    real(kind=8), intent(in) :: aa_t(6, 6, 2), ab(6, 2, 2), ac(2, 2, 2)
    real(kind=8), intent(in) :: ga_t(6, 6, 2), gb(6, 2, 2), gc(2, 2, 2)
    real(kind=8), intent(in) :: aa_c(6, 6, 2)
    real(kind=8), intent(in) :: ga_c(6, 6, 2)
    character(len=16), intent(in) :: option
    aster_logical, intent(in) :: debug
    integer(kind=8), intent(out) :: codret
    real(kind=8), intent(out) :: sig(8), vip(*)
    real(kind=8), intent(out) :: dsidep(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
!      LOI GLOBALE POUR LES PLAQUES/COQUES DKTG - DHRC
!
! --------------------------------------------------------------------------------------------------
!
!       A0     : RAIDEUR ELASTIQUE (D=0)
!       AA     : PARAMETRE ALPHA DE LA FONCTION D'ENDOMMAGEMENT
!       GA     : PARAMETRE GAMMA DE LA FONCTION D'ENDOMMAGEMENT
!       AB     : PARAMETRE ALPHA DE LA FONCTION D'ENDOMMAGEMENT
!       GB     : PARAMETRE GAMMA DE LA FONCTION D'ENDOMMAGEMENT
!       C0     : RAIDEUR ELASTIQUE (D=0)
!       AC     : PARAMETRE ALPHA DE LA FONCTION D'ENDOMMAGEMENT
!       GC     : PARAMETRE GAMMA DE LA FONCTION D'ENDOMMAGEMENT
!       VIM     : VARIABLES INTERNES EN T-
!       OPTION  : TOUTES
!       CRIT    : CRITERES DE CONVERGENCE LOCAUX
!              (1) = NB ITERATIONS MAXI A CONVERGENCE
!                    (ITER_INTE_MAXI == ITECREL)
!              (2) = TYPE DE JACOBIEN A T+DT
!                    (TYPE_MATR_COMP == MACOMP)
!                     0 = EN VITESSE     >SYMETRIQUE
!                     1 = EN INCREMENTAL >NON-SYMETRIQUE
!              (3) = VALEUR TOLERANCE DE CONVERGENCE
!                    (RESI_INTE == RESCREL)
!  !
! IN:            (5) = NOMBRE D'INCREMENTS POUR LE
!                    REDECOUPAGE LOCAL DU PAS DE TEMPS
!                    (ITER_INTE_PAS  == ITEDEC)
!                    -1,0,1 = PAS DE REDECOUPAGE
!                     N = NOMBRE DE PALIERS
!              (6) = TYPE D INTEGRATION LOCAL POUR LA LOI DE
!                    COMPORTEMENT (ALGO_INTE)
! OUT:
!       SIG     : CONTRAINTE
!               (NXX NYY NXY MXX MYY MXY)
!       VIP     : VARIABLES INTERNES EN T+
!       CODRET  : CODE RETOUR DE L'INTEGRATION DE LA LDC
!                 0 => PAS DE PROBLEME
!                 1 => ABSENCE DE CONVERGENCE
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: rigi, resi
    aster_logical :: lelas
    integer(kind=8) :: k, itmax, indi(6), nbact, l, i, iret, iter, iter2
    blas_int :: info
    real(kind=8) :: wr(6), wi(6), work(18), vl(1), vr(1)
    real(kind=8) :: eps(8)
    real(kind=8) :: vint(6)
    real(kind=8) :: a(6, 6), b(6, 2, 2), c(2, 2, 2), ates(6, 6)
    real(kind=8) :: ap1(6, 6), ap2(6, 6), as1(6, 6), as2(6, 6)
    real(kind=8) :: bp1(6, 2), bp2(6, 2), bs1(6, 2), bs2(6, 2)
    real(kind=8) :: cp1(2, 2), cp2(2, 2), cs1(2, 2), cs2(2, 2)
    real(kind=8) :: seuils(6), seuact(6), told
    real(kind=8) :: g1, g2
    real(kind=8) :: neta1(2), neta2(2)
    real(kind=8) :: jacob(6, 6), bocaj(6, 6), det
    real(kind=8) :: epsg(8), sigg(8), dsideg(6, 6)
    real(kind=8) :: dsidem(3, 3), dsidec(3, 3), dsidef(3, 3)
    real(kind=8) :: dsidmg(3, 3), dsidcg(3, 3), dsidfg(3, 3)
    real(kind=8) :: xab1(3, 3)
    blas_int :: b_lda, b_ldvl, b_ldvr, b_lwork, b_n
!
! --------------------------------------------------------------------------------------------------
!
    rigi = (option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL')
    resi = (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL')
    lelas = option .eq. 'RIGI_MECA       '

! - INITIALISATION
    if (lelas) then
        call r8inir(6, 0.0d0, epsm, 1)
        call r8inir(6, 0.0d0, deps, 1)
    end if
!
! --  CALCUL DES EPSILON INITIAUX
    if (resi) then
        do k = 1, 6
            eps(k) = epsm(k)+deps(k)
        end do
    else
        do k = 1, 6
            eps(k) = epsm(k)
        end do
    end if
!
    eps(7) = 0.0d0
    eps(8) = 0.0d0
    eps(3) = eps(3)*0.5d0
    eps(6) = eps(6)*0.5d0
    call r8inir(8, 0.0d0, epsg, 1)
    call dxefro(1, plateOrie%t2iu, eps, epsg)
!
    epsg(3) = epsg(3)*2.d0
    epsg(6) = epsg(6)*2.d0
!
    if (debug) then
        write (6, *) 'eps  :', eps
        write (6, *) 'epsg :', epsg
    end if
!
! ---------------------------------------------------------------------
!
! --  INITIALISATION DE D1, D2, EPSP1 ET EPSP2
! --  STOCKAGE DES VARIABLES INTERNES DANS UN VECTEUR VINT
! --  VINT=(D1,D2,EPSP1X,EPSP1Y,EPSP2X,EPSP2Y)
!
    if (lelas) then
        call r8inir(6, 0.0d0, vint, 1)
    else
        do k = 1, 6
            vint(k) = vim(k)
        end do
    end if
!
    if (debug) then
        write (6, *) 'vint :', vint
    end if
!
    do k = 1, 6
        indi(k) = 0
    end do
!
    if (resi) then
        iter = 0
!
        told = crit(3)
        itmax = nint(crit(1))
!
! --  CALCUL DES TENSEURS DE RAIDEUR A,B,C EN FONCTION DE
!     L'ENDOMMAGEMENT ET DE LEURS DERIVEES PAR RAPPORT A D1 ET D2
!
        call dhrc_calc_b(ab, gb, vint, b, bp1, &
                         bp2, bs1, bs2)
        call dhrc_calc_c(c0, ac, gc, vint, c, &
                         cp1, cp2, cs1, cs2)
        call dhrc_calc_a(a0, aa_t, ga_t, aa_c, ga_c, &
                         epsg, vint, a, ap1, ap2, &
                         as1, as2)
!
! ----------------------------------------------------------------------
! -------CALCUL DES FORCES THERMODYNAMIQUES -------
! ----------------------------------------------------------------------
!
        call dhrc_calc_n(epsg, vint, b, c, neta1, &
                         neta2)
        call dhrc_calc_g(epsg, vint, ap1, bp1, cp1, &
                         ap2, bp2, cp2, g1, g2)
!
! ----------------------------------------------------------------------
! -------CALCUL DES SEUILS-------
! ----------------------------------------------------------------------
!
!     SEUILS D'ENDOMMAGEMENT
!
        call r8inir(6, 0.0d0, seuils, 1)
!
        seuils(1) = g1/cstseu(1)-1.0d0
        seuils(2) = g2/cstseu(2)-1.0d0
!
!     SEUILS DE PLASTICITE
!
        seuils(3) = (neta1(1)/cstseu(3))**2.0d0-1.0d0
        seuils(4) = (neta1(2)/cstseu(4))**2.0d0-1.0d0
        seuils(5) = (neta2(1)/cstseu(5))**2.0d0-1.0d0
        seuils(6) = (neta2(2)/cstseu(6))**2.0d0-1.0d0
!
        if (debug) then
            write (6, *) 'seuils :', seuils
        end if
!
222     continue
!
! --  COMPTEUR D'ITERATIONS
        iter = iter+1
!
        if (iter .gt. itmax) then
            if (debug) then
                write (6, *) 'Non convergence iter1'
            end if
            codret = 1
            goto 999
        end if
!
! --  NOMBRE DE SEUILS ACTIVES
        nbact = 0
!
! --  CREATION DE L'INDICATRICE DES SEUILS ACTIVES
        do k = 1, 6
            if (seuils(k) .gt. told) then
                indi(k) = k
                nbact = nbact+1
            end if
        end do
!
        if (debug) then
            write (6, *) 'Iter 1 :', iter
            write (6, *) 'indi :', indi
        end if
!
! --  SI PAS DE SEUILS ATTEINTS, ON PASSE DIRECTEMENT AU CALCUL DES
!     CONTRAINTES
        if (nbact .eq. 0) then
            goto 555
        else
!
! --  SINON
!
! --  CREATION DU VECTEUR SEUILS ACTIVES
            call r8inir(6, 0.0d0, seuact, 1)
            nbact = 0
            do k = 1, 6
                if (k .eq. indi(k)) then
                    nbact = nbact+1
                    seuact(nbact) = seuils(k)
                end if
            end do
!
            if (debug) then
                write (6, *) 'seuact :', seuact
            end if
!
            iter2 = 0
!-----------------------------------------------------------------------
! --  BOUCLE DE RESOLUTION DE L'EVOLUTION DES VARIABLES --
! ----------------------------------------------------------------------
111         continue
!
! --  COMPTEUR D'ITERATIONS
            iter2 = iter2+1
!
            if (debug) then
                write (6, *) 'Iter 2 :', iter2
                write (6, *) 'seuils :', seuils
            end if
!
            if (iter2 .gt. itmax) then
                if (debug) then
                    write (6, *) 'Non convergence iter2'
                end if
                codret = 1
                goto 999
            end if
!
! --  CALCUL DE LA JACOBIENNE => JACOB(NBACT,NBACT)
!
            jacob(:, :) = 0.0d0
!
            call dhrc_jacob(epsg, vint, c, bp1, cp1, &
                            bp2, cp2, as1, bs1, cs1, &
                            as2, bs2, cs2, indi, neta1, &
                            neta2, cstseu, jacob)
!
!
! --  VERIFICATION DE LA CONVEXITE DE L'ENERGIE LIBRE
! --  POUR CE FAIRE ON REGARDE SI LA JACOBIENNE EST DEFINIE POSITIVE
!
! --  INVERSION DE LA JACOBIENNE => BOCAJ(NBACT,NBACT)
!
            bocaj(:, :) = 0.0d0
!
            do k = 1, nbact
                bocaj(k, k) = 1.0d0
            end do
!
            call mgauss('NFSP', jacob, bocaj, 6, nbact, &
                        6, det, iret)
!
! --  MISE A JOUR DES VARIABLES INTERNES
            l = 0
            do k = 1, 6
                if (k .eq. indi(k)) then
                    l = l+1
                    do i = 1, nbact
                        vint(k) = vint(k)-bocaj(l, i)*seuact(i)
                    end do
                end if
            end do
!
! --  VERIFICATION DE LA CROISSANCE DES D1 ET D2
            if (vint(1) .lt. vim(1)) then
                if (debug) then
                    write (6, *) 'd1m>d1p :'
                    write (6, *) 'd1m', vim(1)
                    write (6, *) 'd1p', vint(1)
                    write (6, *) 'ap1 :'
                    write (6, *) ap1(:, :)
                end if
                vint(1) = vim(1)
            end if
!
            if (vint(2) .lt. vim(2)) then
                if (debug) then
                    write (6, *) 'd2m>d2p'
                    write (6, *) 'd2m', vim(2)
                    write (6, *) 'd2p', vint(2)
                    write (6, *) 'ap2 :'
                    write (6, *) ap2(:, :)
                end if
                vint(2) = vim(2)
            end if
!
            if ((indi(1) .eq. 1) .or. (indi(2) .eq. 2)) then
! --  CALCUL DES TENSEURS DE RAIDEUR A,B,C EN FONCTION DE
!     L'ENDOMMAGEMENT ET DE LEURS DERIVEES PAR RAPPORT A D1 ET D2
!
                call dhrc_calc_b(ab, gb, vint, b, bp1, &
                                 bp2, bs1, bs2)
                call dhrc_calc_c(c0, ac, gc, vint, c, &
                                 cp1, cp2, cs1, cs2)
                call dhrc_calc_a(a0, aa_t, ga_t, aa_c, ga_c, &
                                 epsg, vint, a, ap1, ap2, &
                                 as1, as2)
            end if
!
! --  CALCUL DES SEUILS AVEC VARIABLES ACTUALISEES
!
! ----------------------------------------------------------------------
! ----CALCUL DES FORCES THERMODYNAMIQUES AVEC VARIABLES ACTUALISEES-----
! ----------------------------------------------------------------------
!
            call dhrc_calc_n(epsg, vint, b, c, neta1, &
                             neta2)
            call dhrc_calc_g(epsg, vint, ap1, bp1, cp1, &
                             ap2, bp2, cp2, g1, g2)
!
! ----------------------------------------------------------------------
! -------CALCUL DES SEUILS AVEC VARIABLES ACTUALISEES-------
! ----------------------------------------------------------------------
!
!     SEUILS D'ENDOMMAGEMENT
!
            seuils(1) = g1/cstseu(1)-1.0d0
            seuils(2) = g2/cstseu(2)-1.0d0
!
!     SEUILS DE PLASTICITE
!
            seuils(3) = (neta1(1)/cstseu(3))**2.0d0-1.0d0
            seuils(4) = (neta1(2)/cstseu(4))**2.0d0-1.0d0
            seuils(5) = (neta2(1)/cstseu(5))**2.0d0-1.0d0
            seuils(6) = (neta2(2)/cstseu(6))**2.0d0-1.0d0
!
            if (debug) then
                write (6, *) 'vint :', vint(1:6)
            end if
!
            l = 0
            do k = 1, 6
                if (k .eq. indi(k)) then
                    l = l+1
                    seuact(l) = seuils(k)
                end if
            end do
!
            do i = 1, nbact
                if (seuact(i) .gt. told) then
                    goto 111
                end if
            end do
!
! ----------------------------------------------------------------------
        end if
!
! ----------------------------------------------------------
! -- VERIFICATION QUE D'AUTRES SEUILS NE SONT PAS ACTIVES --
! ----------------------------------------------------------
        goto 222
!
555     continue
!
        do k = 1, 6
            vip(k) = vint(k)
        end do
! --  CALCUL DE LA DISSIPATION
        vip(7) = (vip(1)*cstseu(1)+vip(2)*cstseu(2))
        vip(8) = vim(8)+(abs(vip(3)-vim(3))*cstseu(3)+abs(vip(4)-vim(4))*cstseu(4)+abs(vip(5)-vi&
                 &m(5))*cstseu(5)+abs(vip(6)-vim(6))*cstseu(6))
        vip(10) = 1.d0-( &
                  a(1, 1)*a(2, 2)*a(3, 3))**(1.d0/3.d0)/(a0(1, 1)*a0(2, 2)*a0(3, 3))**(1.d0/3.d0 &
                                                                                       )
        vip(11) = 1.d0-( &
                  a(4, 4)*a(5, 5)*a(6, 6))**(1.d0/3.d0)/(a0(4, 4)*a0(5, 5)*a0(6, 6))**(1.d0/3.d0 &
                                                                                       )
!
    else
        do k = 1, 11
            if (lelas) then
                vip(k) = 0.0d0
            else
                vip(k) = vim(k)
            end if
        end do
    end if
!
    do k = 1, 2
        if (vip(k) .lt. vim(k)) then
            if (debug) then
                write (6, *) 'Evolution non positive de l endommagement'
                write (6, *) 'vim :', vim(1:9)
                write (6, *) 'vip :', vip(1:9)
            end if
            codret = 1
            goto 999
        end if
    end do
!
    call dhrc_calc_b(ab, gb, vip, b, bp1, &
                     bp2, bs1, bs2)
    call dhrc_calc_c(c0, ac, gc, vip, c, &
                     cp1, cp2, cs1, cs2)
    call dhrc_calc_a(a0, aa_t, ga_t, aa_c, ga_c, &
                     epsg, vip, a, ap1, ap2, &
                     as1, as2)
    call dhrc_calc_n(epsg, vip, b, c, neta1, &
                     neta2)
!
    if (resi) then
! --  CALCUL DES CONTRAINTES
        call dhrc_sig(epsg, vip, a, b, sigg)
    end if
!
! ----------------------------------------------------------------------
! --  CALCUL DE LA MATRICE TANGENTE
! ----------------------------------------------------------------------
!
    nbact = 0
    do k = 1, 6
        if (indi(k) .eq. k) then
            nbact = nbact+1
        end if
    end do
!
    if ((.not. rigi) .or. (nbact .eq. 0)) then
!
! --  CALCUL DE LA MATRICE ELASTIQUE
        do k = 1, 6
            do i = 1, 6
                dsideg(k, i) = a(k, i)
            end do
        end do
!
    else
!
! --  CALCUL DE LA MATRICE TANGENTE
!
! --  CALCUL DE LA JACOBIENNE => JACOB(NBACT,NBACT)
!
        jacob(:, :) = 0.0d0
!
        call dhrc_jacob(epsg, vip, c, bp1, cp1, &
                        bp2, cp2, as1, bs1, cs1, &
                        as2, bs2, cs2, indi, neta1, &
                        neta2, cstseu, jacob)
!
! --  INVERSION DE LA JACOBIENNE => BOCAJ(NBACT,NBACT)
!
        bocaj(:, :) = 0.0d0
!
        do k = 1, nbact
            bocaj(k, k) = 1.0d0
        end do
!
        call mgauss('NFSP', jacob, bocaj, 6, nbact, &
                    6, det, iret)
!
        call dhrc_mat_tan(a, ap1, ap2, b, bp1, &
                          bp2, bocaj, neta1, neta2, indi, &
                          cstseu, epsg, vip, dsideg)
    end if
!
    ates(:, :) = 0.0d0
!
    do k = 1, 6
        do i = 1, 6
            ates(k, i) = dsideg(k, i)
        end do
    end do
!
    if (debug) then
        write (6, *) 'dsideg :', dsideg
    end if
!
    b_ldvr = to_blas_int(1)
    b_ldvl = to_blas_int(1)
    b_lda = to_blas_int(6)
    b_n = to_blas_int(6)
    b_lwork = to_blas_int(18)
    call dgeev('N', 'N', b_n, ates, b_lda, &
               wr, wi, vl, b_ldvl, vr, &
               b_ldvr, work, b_lwork, info)
!
!     ECRITURE DES VALEURS PROPRES
!
    do k = 1, 6
        if (wr(k) .lt. 0.0d0) then
            call utmess('A', 'COMPOR4_71', si=k, sr=wr(k))
        end if
    end do
!
! ---  PASSAGE DES EFFORTS GENERALISES SIGG DU REPERE GLOBAL DE LA COQUE
! ---  AU REPERE INTRINSEQUE A L'ELEMENT => SIG
    if (resi) then
        call r8inir(8, 0.0d0, sig, 1)
        call dxefro(1, plateOrie%t2ui, sigg, sig)
    end if
!
    dsidmg(:, :) = 0.0d0
    dsidcg(:, :) = 0.0d0
    dsidfg(:, :) = 0.0d0
!
    do k = 1, 3
        do l = 1, 3
!     MEMBRANE
            dsidmg(k, l) = dsideg(k, l)
!     COUPLAGE MEMBRANE-FLEXION
            dsidcg(k, l) = dsideg(k, l+3)
!     FLEXION
            dsidfg(k, l) = dsideg(k+3, l+3)
        end do
    end do
!
    call utbtab('ZERO', 3, 3, dsidmg, plateOrie%t1ve, &
                xab1, dsidem)
    call utbtab('ZERO', 3, 3, dsidcg, plateOrie%t1ve, &
                xab1, dsidec)
    call utbtab('ZERO', 3, 3, dsidfg, plateOrie%t1ve, &
                xab1, dsidef)
!
    dsidep(:, :) = 0.0d0
!
    do k = 1, 3
        do l = 1, 3
            dsidep(k, l) = dsidem(k, l)
            dsidep(k, l+3) = dsidec(k, l)
            dsidep(l+3, k) = dsidec(k, l)
            dsidep(k+3, l+3) = dsidef(k, l)
        end do
    end do
!
    vip(9) = vip(7)+vip(8)
!
999 continue
!
!A DECOMMENTER POUR DEBUG DE L INTEGRATION DE LA LOI
!    if (codret .gt. 0) then
!        if (debug .eq. .false._1) then
!            write(6,*) 'Echec de l integration de la loi de comportement'
!            write(6,*) 'rigi  =', rigi
!            write(6,*) 'resi  =', resi
!            write(6,*) 'lelas =', lelas
!            write(6,*) 'alpha =', plateOrie%alpha
!            write(6,*) 'beta  =', plateOrie%beta
!            write(6,*) 't2ev2 =', plateOrie%t2ui
!            write(6,*) 't2ve2 =', plateOrie%t2iy
!            write(6,*) 'cosi  =', plateOrie%c
!            write(6,*) 'sinu  =', plateOrie%s
!            write(6,*) 'Donnees d entree'
!            write(6,*) 'epsm    = [',epsm   ,']'
!            write(6,*) 'deps    = [',deps   ,']'
!            write(6,*) 'epsg    = [',epsg   ,']'
!            write(6,*) 'vim     = [',vim(1:9)    ,']'
!            write(6,*) 'cstseu  = [',cstseu ,']'
!            write(6,*) 'crit    = [',crit(1:6)   ,']'
!            write(6,*) 'Matrice A'
!            write(6,*) 'a0  = [',a0 ,']'
!            write(6,*) 'aa_t  = [',aa_t ,']'
!            write(6,*) 'ga_t  = [',ga_t ,']'
!            write(6,*) 'aa_c  = [',aa_c ,']'
!            write(6,*) 'ga_c  = [',ga_c ,']'
!            write(6,*) 'Matrice B'
!            write(6,*) 'ab  = [',ab ,']'
!            write(6,*) 'gb  = [',gb ,']'
!            write(6,*) 'Matrice C'
!            write(6,*) 'c0  = [',c0 ,']'
!            write(6,*) 'ac  = [',ac ,']'
!            write(6,*) 'gc  = [',gc ,']'
!            write(6,*) 'Donnees de sortie'
!            write(6,*) 'vip  = [',vip(1:9) ,']'
!            write(6,*) 'sig  = [',sig ,']'
!            write(6,*) '-------------------------------------'
!            write(6,*) 'On rejout l integration en mode debug'
!        end if
!    end if
!
end subroutine
