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

module nmcine_line_gc_module

    implicit None

    real(kind=8), private :: DEUXMUP, TROISKP, SIGYP, RPRIM, PM, SIGEL(6), PRAGP

    private :: critere

contains

    subroutine nmcine_line_gc(fami, kpg, ksp, ndim, typmod, &
                              imate, relaComp, crit, epsm, deps, &
                              sigm, vim, option, sigp, vip, &
                              dsidep, iret)
        !
        ! ------------------------------------------------------------------------------------------
        !
        !     Écrouissage cinématique en contraintes planes
        !
        !     Très fortement inspiré de NMECMI
        !
        ! ------------------------------------------------------------------------------------------
        !
        ! IN
        !   fami        famille de point de gauss (rigi,mass,...)
        !   kpg, ksp    numéro du (sous)point de gauss
        !   ndim        dimension de l espace (3d=3,2d=2,1d=1)
        !   typmod      type de modélisation
        !   imate       adresse du matériau code
        !   crit        critères  locaux
        !                   crit(1) = nombre d'itérations maxi a convergence
        !                               iter_inte_maxi == itecrel
        !                   crit(2) = type de jacobien a t+dt (type_matr_comp == macomp)
        !                                   0 = en vitesse     > symétrique
        !                                   1 = en incrémental > non-symétrique
        !                   crit(3) = valeur de la tolérance de convergence (resi_inte == rescrel)
        !                   crit(5) = nombre d'incréments pour le redécoupage local du pas de temps
        !                             (iter_inte_pas == itedec)
        !                                   0 = pas de redécoupage
        !                                   n = nombre de paliers
        !   deps        incrément de déformation totale
        !   sigm        contrainte  à t-
        !   epsm        déformation à t-
        !   vim         variables internes a t-
        !                   attention "vim" variables internes a t- modifiées si redécoupage local
        !   option      option de calcul
        !                   'rigi_meca_tang'> dsidep(t)
        !                   'full_meca'     > dsidep(t+dt) , sig(t+dt)
        !                   'raph_meca'     > sig(t+dt)
        !
        ! OUT
        !   sigp        contrainte a t+
        !   vip         variables internes a t+
        !   dsidep      matrice de comportement tangent
        !   iret        code retour
        !                   0   Tout va bien
        !                   1   Redécoupage global ?
        !                   2   Redécoupage local  ?
        !
        !
        ! ------------------------------------------------------------------------------------------
        !
        !   Attention les tenseurs et matrices sont rangés dans l'ordre
        !       XX, YY, ZZ, SQRT(2)*XY, SQRT(2)*XZ, SQRT(2)*YZ
        !
        !   Codage en dur de "pragm" et "PRAGP" car pas d'écrouissage isotrope [R5.03.16 §5]
        !
        ! ------------------------------------------------------------------------------------------
        implicit none
        !
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/radial.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/zerofr.h"
        !
        integer(kind=8)         :: kpg, ksp, ndim, imate, iret
        real(kind=8)    :: crit(10)
        real(kind=8)    :: epsm(6), deps(6), sigm(6), sigp(6), dsidep(6, 6)
        !
        integer(kind=8), parameter  :: nbvari = 12
        real(kind=8)        :: vim(nbvari), vip(nbvari)
        !
        character(len=*)  :: fami
        character(len=8)  :: typmod(*)
        character(len=16) :: relaComp, option
        !
        ! ------------------------------------------------------------------------------------------
        !   Nom des index des variables internes pour VMIS_CINE_GC et VMIS_CINE_LINE
        integer(kind=8)  :: icels, icelu, iepsq, iplas, idiss, iwthe
        integer(kind=8)  :: ixxm, iyym, izzm, ixym, ixzm, iyzm
        integer(kind=8)  :: ivari
        !
        !   Variables globales déclarées dans le module en private
        !       DEUXMUP, TROISKP, SIGYP, RPRIM, PM, SIGEL(6), PRAGP
        ! ------------------------------------------------------------------------------------------
        !
        aster_logical   :: plasti, cinegc, cineli

        integer(kind=8) :: ndimsi, kk, ll, niter, ibid
        !
        real(kind=8)    :: sigedv(6), depsdv(6), sigmdv(6), sigpdv(6), sigdv(6), sigmp(6)
        real(kind=8)    :: depsth(6), xp(6), xm(6), epsicplan(6)
        real(kind=8)    :: depsmo, sigmmo, rp, hp, gp, g1, rpm, sgels, epelu, plast
        real(kind=8)    :: sieleq, sigeps, seuil, dp, cc, radi, prec, dx
        real(kind=8)    :: hsg, pp, precr, epsthe, sigmels, epsielu, epsmo
        real(kind=8)    :: em, num, deuxmum, troiskm, pragm, dsdem
        real(kind=8)    :: ep, nup, dsdep
        ! ------------------------------------------------------------------------------------------
        real(kind=8)            :: dp0, xap
        real(kind=8), parameter :: kron(6) = (/1.0, 1.0, 1.0, 0.0, 0.0, 0.0/)
        !
        ! ------------------------------------------------------------------------------------------
        real(kind=8)        :: valrm(3)
        character(len=16)   :: valkm(3)
        ! ------------------------------------------------------------------------------------------
        integer(kind=8)     :: icodre(4)
        real(kind=8)        :: valres(4)
        character(len=16)   :: nomres(4)
        ! ------------------------------------------------------------------------------------------
        !
        !   Est ce que l'on est bon
        cinegc = relaComp(1:12) .eq. 'VMIS_CINE_GC'
        cineli = relaComp(1:14) .eq. 'VMIS_CINE_LINE'
        !
        ! index des variables internes
        if (cinegc) then
            ! Pour VMIS_CINE_GC
            !       'CRITSIG', 'CRITEPS', 'EPSPEQ', 'INDIPLAS', 'DISSIP', 'DISSTHER',
            !       'XCINXX',  'XCINYY',  'XCINZZ', 'XCINXY', 'XCINXZ', 'XCINYZ',
            icels = 1; icelu = 2; iepsq = 3; iplas = 4; idiss = 5; iwthe = 6
            ixxm = 7; iyym = 8; izzm = 9; ixym = 10; ixzm = 11; iyzm = 12
            ivari = 12
        else if (cineli) then
            ! Pour VMIS_CINE_LINE
            !       'XCINXX', 'XCINYY', 'XCINZZ', 'XCINXY', 'XCINXZ', 'XCINYZ',
            !       'INDIPLAS', 'EPSPEQ'
            icels = 0; icelu = 0; iepsq = 0; iplas = 7; idiss = 0; iwthe = 0
            ixxm = 1; iyym = 2; izzm = 3; ixym = 4; ixzm = 5; iyzm = 6
            ivari = 7
        else
            ASSERT(ASTER_FALSE)
        end if
        ASSERT(typmod(1) .eq. 'C_PLAN')
        ASSERT(2*ndim .eq. 4)
        ! Initialisations
        ndimsi = 2*ndim
        iret = 0
        epsicplan = 0.d0
        !
        ! Mise au format des contraintes de rappel
        if (iepsq .ne. 0) then
            PM = vim(iepsq)
            pp = vim(iepsq)
        else
            PM = 0.d0
            pp = 0.d0
        end if
        plast = vim(iplas)
        ! Cinématique
        xm(1) = vim(ixxm)
        xm(2) = vim(iyym)
        xm(3) = vim(izzm)
        xm(4) = vim(ixym)*sqrt(2.0)
        dp = 0.0
        !
        sgels = 0.0
        epelu = 0.0
        ! Lecture des caractéristiques élastiques du matériau (temps- et +)
        nomres(1) = 'E'; nomres(2) = 'NU'
        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 2)
        em = valres(1)
        num = valres(2)
        deuxmum = em/(1.0+num)
        troiskm = em/(1.0-2.0*num)
        !
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 2)
        ep = valres(1)
        nup = valres(2)
        DEUXMUP = ep/(1.0+nup)
        TROISKP = ep/(1.0-2.0*nup)
        !
        call verift(fami, kpg, ksp, 'T', imate, epsth_=epsthe)
        !
        ! Lecture des caractéristiques d'écrouissage
        dsdem = 0.0
        if (cinegc) then
            nomres(1) = 'D_SIGM_EPSI'; nomres(2) = 'SY'
            nomres(3) = 'SIGM_LIM'; nomres(4) = 'EPSI_LIM'
            call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                        4, nomres, valres, icodre, 1)
            ! vérification que SIGM_LIM, EPSI_LIM sont présents
            if (icodre(3)+icodre(4) .ne. 0) then
                valkm(1) = 'VMIS_CINE_GC'
                valkm(2) = nomres(3)
                valkm(3) = nomres(4)
                call utmess('F', 'COMPOR1_76', nk=3, valk=valkm)
            end if
            !
            dsdem = valres(1)
            sgels = valres(3)
            epelu = valres(4)
        else if (cineli) then
            nomres(1) = 'D_SIGM_EPSI'; nomres(2) = 'SY'
            call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                        2, nomres, valres, icodre, 1)
            !
            dsdem = valres(1)
        end if
        if (dsdem .le. 0.0) then
            valrm(1) = dsdem
            valrm(2) = em
            call utmess('F', 'COMPOR1_53', nr=2, valr=valrm)
        end if
        if ((em-dsdem) .lt. r8miem()) then
            valrm(1) = dsdem
            valrm(2) = em
            call utmess('F', 'COMPOR1_54', nr=2, valr=valrm)
        end if
        !
        nomres(1) = 'D_SIGM_EPSI'; nomres(2) = 'SY'
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 1)
        dsdep = valres(1)
        SIGYP = valres(2)
        if (dsdep .le. 0.0) then
            valrm(1) = dsdep
            valrm(2) = ep
            call utmess('F', 'COMPOR1_53', nr=2, valr=valrm)
        end if
        if ((ep-dsdep) .lt. r8miem()) then
            valrm(1) = dsdep
            valrm(2) = ep
            call utmess('F', 'COMPOR1_54', nr=2, valr=valrm)
        end if
        !
        ! On code en dur ==> RPRIM = 0.0
        pragm = 2.0*dsdem/(1.0-dsdem/em)/3.0
        PRAGP = 2.0*dsdep/(1.0-dsdep/ep)/3.0
        ! Dans le cas d'un écrouissage cinématique ==> RPRIM = 0.0
        ! On garde RPRIM dans la suite pour être similaire à NMECMI
        ! RPRIM = dsdep*ep/(ep-dsdep)-1.50*PRAGP
        RPRIM = 0.0
        rpm = RPRIM*PM+SIGYP
        !
        ! Calcul de depsmo et depsdv
        epsicplan(1:ndimsi) = epsm(1:ndimsi)+deps(1:ndimsi)
        deps(3) = -nup/(1.0-nup)*(deps(1)+deps(2))+(1.0+nup)/(1.0-nup)*epsthe
        epsicplan(3) = -nup/(1.0-nup)*(epsicplan(1)+epsicplan(2))
        depsmo = 0.0
        epsmo = 0.0
        do kk = 1, 3
            depsth(kk) = deps(kk)-epsthe
            depsmo = depsmo+depsth(kk)
            epsmo = epsmo+epsicplan(kk)
        end do
        depsmo = depsmo/3.0
        epsmo = epsmo/3.0
        do kk = 4, ndimsi
            depsth(kk) = deps(kk)
        end do
        epsielu = 0.0
        do kk = 1, ndimsi
            depsdv(kk) = depsth(kk)-depsmo*kron(kk)
            epsielu = epsielu+(epsicplan(kk)-epsmo*kron(kk))**2
        end do
        epsielu = sqrt(2.0*epsielu/3.0)
        !
        ! Calcul de sigmp
        sigmmo = 0.0
        do kk = 1, 3
            sigmmo = sigmmo+sigm(kk)
        end do
        sigmmo = sigmmo/3.0
        do kk = 1, ndimsi
            sigmp(kk) = DEUXMUP/deuxmum*(sigm(kk)-sigmmo*kron(kk))+TROISKP*sigmmo*kron(kk)/troiskm
        end do
        ! Calcul de sigmmo, sigmdv, SIGEL, sieleq, seuil
        sigmmo = 0.0
        do kk = 1, 3
            sigmmo = sigmmo+sigmp(kk)
        end do
        sigmmo = sigmmo/3.0
        sieleq = 0.0
        do kk = 1, ndimsi
            sigmdv(kk) = sigmp(kk)-sigmmo*kron(kk)
            xm(kk) = xm(kk)*PRAGP/pragm
            SIGEL(kk) = sigmdv(kk)+DEUXMUP*depsdv(kk)-xm(kk)
            sieleq = sieleq+SIGEL(kk)**2
        end do
        sieleq = sqrt(1.50*sieleq)
        sigmels = 0.0
        seuil = sieleq-rpm
        hp = 1.0
        gp = 1.0
        !
        ! Calcul de sigp, sigpdv, vip, dp, rp
        if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
            ! Calcul de dp (et dx si C_PLAN) :
            if (seuil .le. 0.0) then
                plast = 0.
                dp = 0.0
                rp = rpm
            else
                plast = 1.0
                prec = abs(crit(3))
                niter = abs(nint(crit(1)))
                precr = prec*SIGYP
                ! Calcul de l'approximation : dp sans contrainte plane
                dp0 = sieleq-SIGYP-RPRIM*PM
                dp0 = dp0/(RPRIM+1.50*(DEUXMUP+PRAGP))
                xap = dp0
                call zerofr(0, 'DEKKER', critere, 0.d0, xap, precr, niter, dp, iret, ibid)
                if (iret .eq. 1) goto 999
                rp = SIGYP+RPRIM*(PM+dp)
            end if
            pp = PM+dp
            gp = 1.0+1.50*PRAGP*dp/rp
            hp = gp+1.50*DEUXMUP*dp/rp
            plasti = (plast .ge. 0.50)
            !
            ! Calcul de sigp
            if (plasti) then
                hsg = hp/gp
                dx = (hsg-1.0)*SIGEL(3)
                dx = dx/(DEUXMUP/1.5+TROISKP*hsg/3.0)
                depsmo = depsmo+dx/3.0
                depsdv(1) = depsdv(1)-dx/3.0
                depsdv(2) = depsdv(2)-dx/3.0
                depsdv(3) = depsdv(3)+dx*2.0/3.0
            end if
            do kk = 1, ndimsi
                sigedv(kk) = sigmdv(kk)+DEUXMUP*depsdv(kk)
                g1 = 1.50*PRAGP*dp/rp/hp
                xp(kk) = xm(kk)*(1.0-g1)+g1*sigedv(kk)
                sigpdv(kk) = sigedv(kk)*gp/hp+xm(kk)*1.50*DEUXMUP*dp/rp/hp
                sigp(kk) = sigpdv(kk)+(sigmmo+TROISKP*depsmo)*kron(kk)
                sigmels = sigmels+sigpdv(kk)**2
            end do
            sigmels = sqrt(1.50*sigmels)
        end if
        !
        ! Calcul de dsidep(6,6)
        if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
            plasti = (plast .ge. 0.50)
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                ! OPTION='RIGI_MECA_TANG' => SIGMA(T)
                do kk = 1, ndimsi
                    sigdv(kk) = sigmdv(kk)-xm(kk)
                end do
                rp = rpm
            else
                ! OPTION='FULL_MECA' => SIGMA(T+DT)
                do kk = 1, ndimsi
                    sigdv(kk) = sigpdv(kk)-xp(kk)
                end do
            end if
            !
            dsidep(:, :) = 0.0
            !
            ! Partie plastique
            sigeps = 0.0
            do kk = 1, ndimsi
                sigeps = sigeps+sigdv(kk)*depsdv(kk)
            end do
            if (plasti .and. sigeps .ge. 0.0) then
                cc = -(1.5*DEUXMUP/rp)**2/(1.5*(DEUXMUP+PRAGP)+RPRIM)*(1.0-dp*RPRIM/rp)/hp
                do kk = 1, ndimsi
                    do ll = 1, ndimsi
                        dsidep(kk, ll) = cc*sigdv(kk)*sigdv(ll)
                    end do
                end do
            end if
            !
            ! Partie élastique
            do kk = 1, 3
                do ll = 1, 3
                    dsidep(kk, ll) = dsidep(kk, ll)+TROISKP/3.0-DEUXMUP*gp/hp/3.0
                end do
            end do
            do kk = 1, ndimsi
                dsidep(kk, kk) = dsidep(kk, kk)+DEUXMUP*gp/hp
            end do
            !
            ! Correction pour les contraintes planes
            cykk: do kk = 1, ndimsi
                if (kk .eq. 3) cycle cykk
                cyll: do ll = 1, ndimsi
                    if (ll .eq. 3) cycle cyll
                    dsidep(kk, ll) = dsidep(kk, ll)-dsidep(kk, 3)*dsidep(3, ll)/dsidep(3, 3)
                end do cyll
            end do cykk
        end if
        !
        if (option(1:9) .ne. 'RIGI_MECA') then
            if (crit(10) .gt. 0.0) then
                call radial(ndimsi, sigm, sigp, vim(iplas), plast, 1, vim(ixxm), vip(ixxm), radi)
                if (radi .gt. crit(10)) then
                    iret = 2
                end if
            end if
        end if
        !
        ! Mise à jour des variables internes
        if ((option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA')) then
            vip(1:ivari) = vim(1:ivari)
            if (iepsq .ne. 0) vip(iepsq) = pp
            vip(iplas) = plast
            vip(ixxm) = xp(1)
            vip(iyym) = xp(2)
            vip(izzm) = xp(3)
            vip(ixym) = xp(4)/sqrt(2.0)
            ! Critère ELS et ELU
            if (cinegc) then
                vip(icels) = sigmels/sgels
                vip(icelu) = epsielu/epelu
            end if
        end if
        !
999     continue
        !
    end subroutine nmcine_line_gc

    ! ----------------------------------------------------------------------------------------------
    ! Variables globales déclarées dans le module en private
    !   DEUXMUP, TROISKP, SIGYP, RPRIM, PM, SIGEL(6), PRAGP
    !

    ! ----------------------------------------------------------------------------------------------
    ! Fonction dont on cherche le zéro pour la plasticité de Von Mises cinématique
    real(kind=8) function critere(xxx)
        implicit none
        real(kind=8) ::  xxx
        !
        !   Variables locales
        real(kind=8) :: rpp, gp, hp, hsg, demuc, dx, fp
        ! ------------------------------------------------------------------------------------------
        rpp = SIGYP+RPRIM*(PM+xxx)
        gp = 1.0+1.5*PRAGP*xxx/rpp
        hp = gp+1.5*DEUXMUP*xxx/rpp
        hsg = hp/gp
        demuc = DEUXMUP+PRAGP
        !
        dx = (hsg-1.0)*SIGEL(3)
        dx = dx/(DEUXMUP/1.5+TROISKP*hsg/3.0)
        !
        fp = (SIGEL(1)-DEUXMUP*dx/3.0)**2+(SIGEL(2)-DEUXMUP*dx/3.0)**2+ &
             (SIGEL(3)+2.0*DEUXMUP*dx/3.0)**2+SIGEL(4)**2
        !
        critere = 1.5*demuc*xxx-sqrt(1.50*fp)+rpp
        !
    end function critere

end module nmcine_line_gc_module
