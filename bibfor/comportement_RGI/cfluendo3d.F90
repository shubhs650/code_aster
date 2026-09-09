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
! aslint: disable=W0413
!
subroutine cfluendo3d(fami, kpg, ksp, ndim, imate, &
                      compor, carcri, instam, instap, epsm, &
                      deps, sigm, vim, option, sigp, &
                      vip, typmod, dsidep, codret)
!=====================================================================
!      Ficher de base de FLUAGE ET ENDOMMAGEMENT
!=====================================================================
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/fluendo3d.h"
#include "asterfort/getRgiPara.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
!
    character(len=*) :: fami
    integer(kind=8) :: kpg, ksp, ndim, imate
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigm(6), vim(*)
    character(len=16) :: option
    real(kind=8) :: sigp(6), vip(*)
    character(len=8) :: typmod(*)
    real(kind=8) :: dsidep(6, 6)
    integer(kind=8) :: codret

!
! DECLARATIONS LOCALES
    integer(kind=8) :: nvarflumax, nmatflumax
    parameter(nvarflumax=158, nmatflumax=165)
!
    integer(kind=8) :: nmatbe, nmatac, nmatflu, nvarflu, nmatbe2, nmatbe3
!   Nombre de paramètres relatifs au beton
    parameter(nmatbe=57, nmatbe2=59)
!   Nombre de paramtre relatif a l'acier (1+21*5)
    parameter(nmatac=106)
!
    integer(kind=8) :: nbelas3d
    parameter(nbelas3d=4)

    character(len=12) :: nomres(nmatflumax), nomres1(nbelas3d)
    real(kind=8) :: valres(nmatflumax), xmat(nbelas3d+nmatflumax), valres1(nbelas3d)
    integer(kind=8) :: nstrs, mfr, i, j
    integer(kind=8) :: retour(nmatflumax), retour1(nbelas3d), iteflumaxi
    real(kind=8) :: d(6, 6), e, nu, coef, coef1, coef2, coef3
    real(kind=8) :: zero, un, deux, rac2, var0(6), sig0(6)
!
!    real(kind=8) :: hydrm, hydrp, sechm, vgm, vgp
    real(kind=8) :: sechp, sref
    real(kind=8) :: alpham, alphap, teta13d, teta23d, sech, epstf3d(6)
!
!
    integer(kind=8) :: iret, ifour, ifour11
    real(kind=8) :: epsmc(6), depsc(6), epsc(6)
    real(kind=8) :: tm, tp, tref

!   variables de transfert de donnees ( a declarer suivant idvar4 et idvisc)
    integer(kind=8) :: nmat3dmax, nmat3d, nstrs3d, nvari3dmax, nvari3d, ierr1, mfr11
!

!   taille du pseudo vecteur des contraintes pour fluendo3d (tjrs 6
!   en raison de son utilisation dans fludes3d qui la suppose à 6)
    parameter(nstrs3d=6)
!
!   mettre a jour ici le nbre de parametres materiaux et variables
!   interne de l option du modele
    integer(kind=8) :: NMATAILX, NMATRAG, NVARRAG
    parameter(NMATRAG=0)
    parameter(NVARRAG=0)
    parameter(NMATAILX=0)
!
!
!   nombre totale de parametres et variables internes option comprise
    parameter(nmat3dmax=nbelas3d+nmatflumax+NMATRAG+NMATAILX)
    parameter(nvari3dmax=nvarflumax+NVARRAG)
!
!   nbre de parametres materiaux sans les tailles des elements
    real(kind=8) :: xmat3d(nmat3dmax), sig03d(6), sigf3d(6), depst3d(6)
    real(kind=8) :: var03d(nvari3dmax), varf3d(nvari3dmax), varf(nvari3dmax), sigf(6)
!   indicateur d isotropie initiale
    aster_logical :: iso1, local11, end3d, fl3d
!   temperatures debut et fin de pas , moyenne, pas de temps, volule rgi
    real(kind=8) :: dt3d
! - Indicateur du type de matrice
    aster_logical :: matrEndo
    real(kind=8)  :: matdech(6, 6)
!
    character(len=16) :: rela_name
!
!   --------------------------------------------------------------------
!   Nombre de paramètres matériau et de variables internes
    rela_name = compor(RELA_NAME)
    if (rela_name(1:12) .ne. 'RGI_BETON_BA') then
        nmatflu = nmatbe
        nvarflu = 114
        nmatbe3 = 0
    else
        nmatflu = nmatbe+2+nmatac
        nvarflu = 158
        nmatbe3 = nmatbe2
    end if
!   nombre totale de paramètres et variables internes (option comprise)
    nmat3d = nbelas3d+nmatflu+NMATRAG+NMATAILX
    nvari3d = nvarflu+NVARRAG

! - Flag for elastic matrix
    matrEndo = nint(carcri(TYPE_MATR_T)) .eq. 4
    if (.not. L_MATR(option)) matrEndo = ASTER_FALSE

    codret = 0
!
!
    sig03d = 0.d0
    sigf3d = 0.d0
    var03d = 0.d0
    varf3d = 0.d0
    varf = 0.d0
    depsc = 0.d0
    epsc = 0.d0
    epsmc = 0.d0
    xmat = 0.d0
    var0 = 0.d0
    sig0 = 0.d0
    sigf = 0.d0
    depst3d = 0.d0
    valres1 = 0.d0
    valres = 0.d0
    matdech = 0.d0
!
    iteflumaxi = int(carcri(ITER_INTE_MAXI))
!
!
! APPEL DE RCVARC POUR LE CALCUL DE LA TEMPERATURE
!
    call rcvarc('f', 'TEMP', '-', fami, kpg, &
                ksp, tm, iret)
    call rcvarc('f', 'TEMP', '+', fami, kpg, &
                ksp, tp, iret)
    call rcvarc('f', 'TEMP', 'REF', fami, kpg, &
                ksp, tref, iret)
!
! ------------------------------------------------
!!     RECUPERATION DE L HYDRATATION DEBUT DE PAS
!    call rcvarc(' ', 'HYDR', '-', fami, kpg,&
!                ksp, hydrm, codret)
!    if (codret .ne. 0) then
!        hydrm=0.d0
!    endif
!!
!! ------------------------------------------------
!!     RECUPERATION DE L HYDRATATION FIN DE PAS
!    call rcvarc(' ', 'HYDR', '+', fami, kpg,&
!                ksp, hydrp, codret)
!    if (codret .ne. 0) then
!        hydrp=0.d0
!    endif
!
! ------------------------------------------------
!     RECUPERATION DU SECHAGE
    call rcvarc(' ', 'SECH', '+', fami, kpg, &
                ksp, sechp, iret)
    if (iret .ne. 0) sechp = 0.d0
!    call rcvarc(' ', 'SECH', '-', fami, kpg,&
!                ksp, sechm, iret)
!    if (iret .ne. 0) sechm=0.d0
    call rcvarc(' ', 'SECH', 'REF', fami, kpg, &
                ksp, sref, iret)
    if (iret .ne. 0) sref = 0.d0
!
    sech = sechp
!
!   le séchage de référence doit être nul
!
    if (sref .ne. 0.d0) call utmess('F', 'COMPOR3_9', sk=rela_name)
!
! -----------------------------------------------
!     RECUPERATION DU VOLUME DE GEL DEBUT DE PAS
!    call rcvarc(' ', 'X1', '-', fami, kpg,&
!                ksp, vgm, codret)
!    if (codret .ne. 0) then
!        vgm=0.d0
!    endif
!
! ------------------------------------------------
!     RECUPERATION DU VOLUME DE GEL FIN DE PAS
!    call rcvarc(' ', 'X1', '+', fami, kpg,&
!                ksp, vgp, codret)
!    if (codret .ne. 0) then
!        vgp=0.d0
!    endif
!
! ------------------------------------------------
!
!   vérification que B_ENDOGE et K_DESSIC sont nuls
    nomres1(1) = 'B_ENDOGE'
    nomres1(2) = 'K_DESSIC'
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres1, valres1, retour1, 0)
    do i = 1, 2
        if (retour1(i) .eq. 0) then
            if (valres1(i) .ne. 0.d0) then
                call utmess('F', 'COMPOR3_40', nk=2, valk=[rela_name, nomres1(i)])
            end if
        end if
    end do
!
    nomres1(1) = 'E'
    nomres1(2) = 'NU'
    nomres1(3) = 'RHO'
    nomres1(4) = 'ALPHA'
!
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                4, nomres1, valres1, retour1, 2)
!
!
!        MODULES INSTANTANES ISOTROPES
    xmat(1:4) = valres1(1:4)
    alpham = valres1(4)
!
! --- EVALUATION PARAMETERES MATERIAU ELASTIQUES A T+
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                4, nomres1, valres1, retour1, 2)
!
    alphap = valres1(4)
!
! ------------------------------------------------------------------
! --  RETRAIT INCREMENT DE DEFORMATION DUE A LA DILATATION THERMIQUE
! ------------------------------------------------------------------
    if (ndim .eq. 2) then
        nstrs = 4
    else
        nstrs = 6
    end if
!
    if ((option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA')) then
        do i = 1, 3
            depsc(i) = deps(i)-(alphap*(tp-tref)-alpham*(tm-tref))
            epsmc(i) = epsm(i)-alpham*(tm-tref)
            epsc(i) = epsmc(i)+depsc(i)
        end do
        do i = 4, nstrs
            depsc(i) = deps(i)
            epsc(i) = epsm(i)+depsc(i)
        end do
    end if
!
! ------------------------------------------------------------------
! --  RECUPERATION PARAMETRES MATERIAU DU MODELE FLUENDO3D
! ------------------------------------------------------------------
!
    call getRgiPara(nomres, nmatbe2)
!
!
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', rela_name, 0, ' ', [0.d0], &
                nmatflu, nomres, valres, retour, 0, &
                nan='NON')
!
    do i = (nbelas3d+1), (nbelas3d+nmatflu)
        xmat(i) = valres(i-nbelas3d)
    end do
!
!-----VALEUR FIXEE PROVISOIREMENT POUR MFR
    mfr = 1
!-----------------------------------------
!
    if (typmod(1) (1:2) .eq. '3D') then
        ifour = 2
    else
        if (typmod(1) (1:6) .eq. 'D_PLAN') then
            ifour = -1
        else
            ifour = 0
        end if
    end if
    if (rela_name .eq. 'FLUA_PORO_BETON') then
        fl3d = .true.
        end3d = .false.
    else
        if (rela_name .eq. 'ENDO_PORO_BETON') then
            fl3d = .false.
            end3d = .true.
        else
            fl3d = .true.
            end3d = .true.
        end if
    end if
!
!    pas de temps
    dt3d = instap-instam
!
    do i = 1, nmat3d
        xmat3d(i) = xmat(i)
    end do
!
    if (L_CORR(option)) then
        do i = 1, nvari3d
            var03d(i) = vim(i)
            varf3d(i) = vip(i)
        end do
    else
        do i = 1, nvari3d
            var03d(i) = vim(i)
            varf3d(i) = 0.d0
        end do
    end if
!
!   initialisation des contraintes effectives si premier pas
    if (abs(var03d(64)-1.d0) .ge. r8prem()) then
!
        do i = 1, 6
!         initialisation de l etage elastique
            if (i .le. 3) then
!              on retire la pression intra poreuse au cas où
!              elle aurait été initialisée dans les vari
                var03d(18+i) = sig03d(i)-var03d(61)
            else
                var03d(18+i) = sig03d(i)
            end if
!         idem pour l etage de Kelvin
            var03d(49+i) = var03d(18+i)
        end do
!      on met l indicateur de 1er passage a 1
        varf3d(64) = 1.
    else
!     on reconduit l indicateur de 1er passage
        varf3d(64) = var03d(64)
    end if
!
!    autres parametres a renseigner
!    temperature moyenne sur le pas
    teta13d = tm
    teta23d = tp
!    initialisation indicateur d erreur
    ierr1 = 0
!    indicateur isostropie elastique et de resistance
    iso1 = .true.
!    numero de la formulation (33 pour poreux)
    mfr11 = mfr
!    type de formulation
    ifour11 = ifour
!
!
!    controle de regularisation en cas d endommagement
    local11 = .true.
    rac2 = sqrt(2.d0)
    do i = 1, 3
        depst3d(i) = depsc(i)
        epstf3d(i) = epsc(i)
    end do
    do i = 4, nstrs
        depst3d(i) = depsc(i)*rac2
        epstf3d(i) = epsc(i)*rac2
    end do
!
    do i = 1, 3
        sig03d(i) = sigm(i)
    end do
!
    do i = 4, nstrs
        sig03d(i) = sigm(i)/rac2
    end do
!
!
    call fluendo3d(xmat3d, sig03d, sigf3d, depst3d, nstrs3d, &
                   var03d, varf3d, nvari3d, nbelas3d, teta13d, &
                   teta23d, dt3d, epstf3d, ierr1, iso1, &
                   mfr11, end3d, fl3d, local11, ndim, &
                   nmatbe3, iteflumaxi, sech, matrEndo, matdech)
!
    do i = 1, 3
        sigp(i) = sigf3d(i)
    end do
    do i = 4, nstrs
        sigp(i) = sigf3d(i)*rac2
    end do
!
    if (L_CORR(option)) then
        do i = 1, nvari3d
            vip(i) = varf3d(i)
        end do
    end if
!
    if (L_MATR(option)) then
!
        if (matrEndo) then
            dsidep(1:nstrs, 1:nstrs) = matdech(1:nstrs, 1:nstrs)
        else
            zero = 0.d0
            un = 1.d0
            deux = 2.d0
!
            d(:, :) = zero
!
            e = xmat(1)
            nu = xmat(2)
!
            coef = un/((un+nu)*(un-deux*nu))
            coef1 = e*(un-nu)*coef
            coef2 = e*nu*coef
            coef3 = e/(un+nu)
!
            d(1, 1) = coef1
            d(1, 2) = coef2
            d(1, 3) = coef2
!
            d(2, 1) = coef2
            d(2, 2) = coef1
            d(2, 3) = coef2
!
            d(3, 1) = coef2
            d(3, 2) = coef2
            d(3, 3) = coef1
!
            d(4, 4) = 0.5d0*coef3
            d(5, 5) = 0.5d0*coef3
            d(6, 6) = 0.5d0*coef3
!
            do i = 1, nstrs
                do j = 1, nstrs
                    dsidep(i, j) = d(i, j)
                end do
            end do
!
        end if
!
    end if
!
end subroutine
