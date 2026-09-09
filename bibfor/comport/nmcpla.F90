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
subroutine nmcpla(BEHInteg, &
                  fami, kpg, ksp, ndim, typmod, jvMaterCode, &
                  comporPlas, comporCreep, carcri, &
                  timed, timef, neps, epsdt, depst, &
                  nsig, sigd, vind, option, &
                  sigf, vinf, ndsde, dsde, iret)
!
    use calcul_module, only: ca_ctempl_, ca_ctempr_, ca_ctempm_, ca_ctempp_
    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcopil.h"
#include "asterfort/nmcomp.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/Behaviour_type.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8) :: jvMaterCode, ndim, kpg, ksp, iret
    integer(kind=8) :: neps, nsig, ndsde
    character(len=16), intent(in) :: comporPlas(COMPOR_SIZE)
    character(len=16), intent(in) :: comporCreep(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: timed, timef, tempd, tempf, tref
    real(kind=8) :: epsdt(6), depst(6)
    real(kind=8) :: sigd(6), sigf(6)
    real(kind=8) :: vind(*), vinf(*)
    real(kind=8) :: dsde(ndsde)
    character(len=16) :: option
    character(len=*) :: fami
    character(len=8) :: typmod(2)
!
! --------------------------------------------------------------------------------------------------
!
! Comportment management
!
! KIT_DDI with creeping = GRANGER
!
! --------------------------------------------------------------------------------------------------
!
!       INTEGRATION DU COUPLAGE FLUAGE/FISSURATION, C'EST A DIRE LE
!       COUPLAGE D'UNE LOI DE COMPORTEMENT DE TYPE FLUAGE GRANGER
!       ET D'UNE LOI DE  COMPORTEMENT ELASTO PLASTIQUE
!               AVEC    . N VARIABLES INTERNES
!                       . UNE FONCTION SEUIL ELASTIQUE
!
!       INTEGRATION DES CONTRAINTES           = SIG(T+DT)
!       INTEGRATION DES VARIABLES INTERNES    = VIN(T+DT) (CUMUL DES
!              VARIABLES INTERNES DES DEUX LOIS)
!       ET CALCUL DU JACOBIEN ASSOCIE         = DS/DE(T+DT) OU DS/DE(T)
!
!       ================================================================
!       ARGUMENTS
!
!       IN      KPG,KSP  NUMERO DU (SOUS)POINT DE GAUSS
!               NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!               TYPMOD  TYPE DE MODELISATION
!               IMAT    ADRESSE DU MATERIAU CODE
!               COMP    COMPORTEMENT DE L ELEMENT
!                       COMP(1) = RELATION DE COMPORTEMENT
!                       COMP(2) = NB DE VARIABLES INTERNES
!                       COMP(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!               OPT     OPTION DE CALCUL A FAIRE
!                               'RIGI_MECA_TANG'> DSDE(T)
!                               'FULL_MECA'     > DSDE(T+DT) , SIG(T+DT)
!                               'RAPH_MECA'     > SIG(T+DT)
!               CRIT    CRITERES  LOCAUX
!                       CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                                 (ITER_INTE_MAXI == ITECREL)
!                       CRIT(2) = TYPE DE JACOBIEN A T+DT
!                                 (TYPE_MATR_COMP == MACOMP)
!                                 0 = EN VITESSE     > SYMETRIQUE
!                                 1 = EN INCREMENTAL > NON-SYMETRIQUE
!                       CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                                 (RESI_INTE == RESCREL)
!                       CRIT(5) = NOMBRE D'INCREMENTS POUR LE
!                                 REDECOUPAGE LOCAL DU PAS DE TEMPS
!                                 (ITER_INTE_PAS == ITEDEC)
!                                 0 = PAS DE REDECOUPAGE
!                                 N = NOMBRE DE PALIERS
!               TIMED   INSTANT T
!               TIMEF   INSTANT T+DT
!               EPSDT   DEFORMATION TOTALE A T
!               DEPST   INCREMENT DE DEFORMATION TOTALE
!               SIGD    CONTRAINTE A T
!               VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!       OUT     SIGF    CONTRAINTE A T+DT
!               VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSDE    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!               IRET    CODE RETOUR DE L'INTEGRATION INTEGRATION DU
!                       COUPLAGE FLUAGE/FISSURATION
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ABSENCE DE CONVERGENCE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8) :: ndt_local, ndi_local
    integer(kind=8) :: ndt, ndi
    integer(kind=8) :: nvi_flua, nvi_plas, idx_vi_plas
    integer(kind=8) :: i, retcom, ire2, iret1, iret2, iret3
    integer(kind=8) :: nume_plas, nume_flua
    integer(kind=8) :: k, iter, itemax, nvi, numlc
    integer(kind=8) :: cerr(5)
    character(len=8) :: elem_model, nomc(5)
    character(len=16) :: rela_flua, rela_plas
    real(kind=8) :: nu
    real(kind=8) :: espi_creep(6), epsfld(6), epsflf(6), depsfl(6)
    real(kind=8) :: deps(6), kooh(6, 6)
    real(kind=8) :: materd(5), materf(5), depst2(6), depsel(6)
    real(kind=8) :: epsicv, toler, ndsig, nsigf
    real(kind=8) :: dsigf(6), hydrd, hydrf, sechd, sechf, sref
    real(kind=8) :: epseld(6), epself(6), epsthe
    real(kind=8) :: sigf2(6)
    real(kind=8) :: tmpdmx, tmpfmx, epsth
    real(kind=8) :: alphad, alphaf, bendod, bendof, kdessd, kdessf
    aster_logical :: l_inte_forc, l_epsi_varc
!
! --------------------------------------------------------------------------------------------------
!
    common/tdim/ndt, ndi
!
! --------------------------------------------------------------------------------------------------
!
    l_inte_forc = option .eq. 'RAPH_MECA' .or. option .eq. 'FULL_MECA'
    rela_flua = comporCreep(RELA_NAME)
    rela_plas = comporPlas(RELA_NAME)
    elem_model = typmod(1)
    read (comporCreep(NVAR), '(I16)') nvi_flua
    read (comporPlas(NVAR), '(I16)') nvi_plas
    read (comporPlas(NUME), '(I16)') nume_plas
    read (comporCreep(NUME), '(I16)') nume_flua
    ASSERT(comporCreep(RELA_NAME) (1:13) .eq. 'BETON_GRANGER')
    nvi = BEHInteg%behavPara%nvi
    numlc = BEHInteg%behavPara%numlc
!
! - Get size for tensors
!
    ndt_local = 2*ndim
    ndi_local = 3
!
! - Get temperatures
!
    call rcvarc(' ', 'TEMP', '-', fami, kpg, ksp, tempd, iret1)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, ksp, tempf, iret2)
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, ksp, tref, iret3)
!
! - Get maximum temperature during load (for BETON_DOUBLE_DP)
!
    tmpdmx = tempd
    tmpfmx = tempf
    if (rela_plas .eq. 'BETON_DOUBLE_DP') then
        if (((iret1+iret2) .eq. 0) .and. (.not. isnan(vind(nvi_flua+3)))) then
            if (tmpdmx .lt. vind(nvi_flua+3)) then
                tmpdmx = vind(nvi_flua+3)
            end if
            if (tmpfmx .lt. vind(nvi_flua+3)) then
                tmpfmx = vind(nvi_flua+3)
            end if
        end if
    end if
!
! - Set temperature access
!
    ca_ctempl_ = 1
    ca_ctempm_ = tmpdmx
    ca_ctempp_ = tmpfmx
    ca_ctempr_ = tref
!
! - Get convergence criteria
!
    itemax = int(carcri(1))
    toler = carcri(3)
!
! = COUPLING = BEGIN
!
    iter = 1
    depst2(1:6) = depst(1:6)
!
20  continue
!
! - Solve creep law
!
    if (l_inte_forc) then
!
! ----- Solve creep law
!
        l_epsi_varc = ASTER_TRUE
        sigf2 = 0.d0
        BEHInteg%behavPara%nvi = nvi_flua
        BEHInteg%behavPara%numlc = nume_flua
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    timed, timef, &
                    comporCreep, carcri, multComp, &
                    neps, epsdt, depst2, &
                    nsig, sigd, &
                    vind, &
                    sigf2, vinf, &
                    ndsde, dsde, iret, &
                    l_epsi_varc_=l_epsi_varc)
        BEHInteg%behavPara%nvi = nvi
        BEHInteg%behavPara%numlc = numlc
!
! ----- Get material parameters
!
        nomc(1) = 'E       '
        nomc(2) = 'NU      '
        nomc(3) = 'ALPHA   '
        nomc(4) = 'B_ENDOGE'
        nomc(5) = 'K_DESSIC'
        call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                    ' ', 'ELAS', 1, 'TEMP', [tmpdmx], &
                    1, nomc(2), materd(2), cerr(1), 2)
        call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                    ' ', 'ELAS', 1, 'TEMP', [tmpfmx], &
                    1, nomc(2), materf(2), cerr(1), 2)
        materd(1) = 1.d0
        materf(1) = 1.d0
!
! ----- Creep strains - At beginning of step
!
        do k = 1, ndt_local
            espi_creep(k) = vind(8*ndt_local+k)
            do i = 1, 8
                espi_creep(k) = espi_creep(k)-vind((i-1)*ndt_local+k)
            end do
        end do
        call lcopil('ISOTROPE', elem_model, materd, kooh)
        epsfld(1:ndt) = matmul(kooh(1:ndt, 1:ndt), espi_creep(1:ndt))
!
! ----- Creep strains - At end of step
!
        do k = 1, ndt_local
            espi_creep(k) = vinf(8*ndt_local+k)
            do i = 1, 8
                espi_creep(k) = espi_creep(k)-vinf((i-1)*ndt_local+k)
            end do
        end do
        call lcopil('ISOTROPE', elem_model, materf, kooh)
        epsflf(1:ndt) = matmul(kooh(1:ndt, 1:ndt), espi_creep(1:ndt))
!
! ----- Creep strain increment
!
        do k = 1, ndt_local
            depsfl(k) = epsflf(k)-epsfld(k)
        end do
    end if
!
! - Unset temperature access
!
    ca_ctempl_ = 0
!
! - Total strains
!
    if (l_inte_forc) then
        do k = 1, ndt_local
            deps(k) = depst(k)-depsfl(k)
        end do
    else
        do k = 1, ndt_local
            deps(k) = depst(k)
        end do
    end if
!
! - Solve plasticity law
!
    l_epsi_varc = ASTER_FALSE
    idx_vi_plas = nvi_flua+1
    sigf = 0.d0
    BEHInteg%behavPara%nvi = nvi_plas
    BEHInteg%behavPara%numlc = nume_plas
    call nmcomp(BEHInteg, &
                ndim, option, typmod, &
                timed, timef, &
                comporPlas, carcri, multComp, &
                neps, epsdt, deps, &
                nsig, sigd, &
                vind(idx_vi_plas), &
                sigf, vinf(idx_vi_plas), &
                ndsde, dsde, &
                retcom, l_epsi_varc_=l_epsi_varc)
    BEHInteg%behavPara%nvi = nvi
    BEHInteg%behavPara%numlc = numlc
    if (retcom .eq. 1) then
        iret = 1
        goto 999
    end if
!
! - Change in common
!
    ndt = 2*ndim
    ndi = 3
!
! - Coupling algorithm
!
    if (l_inte_forc) then
!
! ----- Get material parameters
!
        call rcvalb(fami, kpg, ksp, '-', jvMaterCode, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    5, nomc(1), materd(1), cerr(1), 2)
        if (cerr(3) .ne. 0) materd(3) = 0.d0
        if (cerr(4) .ne. 0) materd(4) = 0.d0
        if (cerr(5) .ne. 0) materd(5) = 0.d0
        call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    5, nomc(1), materf(1), cerr(1), 2)
        if (cerr(3) .ne. 0) materf(3) = 0.d0
        if (cerr(4) .ne. 0) materf(4) = 0.d0
        if (cerr(5) .ne. 0) materf(5) = 0.d0
!
! ----- Elastic strain increment
!
        call lcopil('ISOTROPE', elem_model, materd, kooh)
        epseld(1:ndt) = matmul(kooh(1:ndt, 1:ndt), sigd(1:ndt))
        call lcopil('ISOTROPE', elem_model, materf, kooh)
        epself(1:ndt) = matmul(kooh(1:ndt, 1:ndt), sigf(1:ndt))
        do k = 1, ndt_local
            depsel(k) = epself(k)-epseld(k)
        end do
!
! --- CALCUL DE L'INCREMENT DE DEFORMATION ELASTIQUE
! --- + RETRAIT ENDOGENNE + RETRAIT DESSICCATION + RETRAIT THERMIQUE
!
        alphad = materd(3)
        alphaf = materf(3)
        bendod = materd(4)
        bendof = materf(4)
        kdessd = materd(5)
        kdessf = materf(5)
        call rcvarc(' ', 'HYDR', '-', fami, kpg, &
                    ksp, hydrd, ire2)
        if (ire2 .ne. 0) hydrd = 0.d0
        call rcvarc(' ', 'HYDR', '+', fami, kpg, &
                    ksp, hydrf, ire2)
        if (ire2 .ne. 0) hydrf = 0.d0
        call rcvarc(' ', 'SECH', '-', fami, kpg, &
                    ksp, sechd, ire2)
        if (ire2 .ne. 0) sechd = 0.d0
        call rcvarc(' ', 'SECH', '+', fami, kpg, &
                    ksp, sechf, ire2)
        if (ire2 .ne. 0) sechf = 0.d0
        call rcvarc(' ', 'SECH', 'REF', fami, kpg, &
                    ksp, sref, ire2)
        if (ire2 .ne. 0) sref = 0.d0
!
        if ((iret1+iret2+iret3) .ne. 0) then
            epsthe = 0.d0
        else
            epsthe = alphaf*(tempf-tref)-alphad*(tempd-tref)
        end if
        epsth = epsthe-bendof*hydrf+bendod*hydrd-kdessf*(sref-sechf)+kdessd*(sref-sechd)
        do k = 1, 3
            depsel(k) = depsel(k)+epsth
        end do
!
! ----- For plane stress
!
        if (elem_model(1:6) .eq. 'C_PLAN') then
            nu = materf(2)
            depsel(3) = -nu/(1.d0-nu)*(depsel(1)+depsel(2))+(1.d0+nu)/(1.d0-nu)*epsth
        end if
!
! ---    CALCUL DE L'INCREMENT DE DEFORMATION EN ENTREE DU CALCUL
! ---    DE FLUAGE POR L'ITERATION SUIVANTE
!
        do k = 1, ndt_local
            depst2(k) = depst2(k)+depsel(k)+depsfl(k)-depst2(k)
        end do
!
! ---    CRITERE DE CONVERGENCE - NORME DE SIGF2 - SIGF
!
        do k = 1, ndt_local
            dsigf(k) = sigf2(k)-sigf(k)
        end do
!
        ndsig = 0.d0
        nsigf = 0.d0
        do k = 1, ndt_local
            ndsig = ndsig+dsigf(k)*dsigf(k)
            nsigf = nsigf+sigf(k)*sigf(k)
        end do
!
        if (nsigf .gt. toler*toler) then
            epsicv = (ndsig/nsigf)**0.5d0
        else
            epsicv = ndsig**0.5d0
        end if
!
        if (epsicv .gt. toler) then
            if (iter .lt. itemax) then
                iter = iter+1
                goto 20
            else
                iret = 1
                goto 999
            end if
        end if
    end if
!
! = COUPLING = END
!
999 continue
end subroutine
