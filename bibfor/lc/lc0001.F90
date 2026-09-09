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
! aslint: disable=C1505,W0104,W1306,W1504

subroutine lc0001(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)

    use Behaviour_type
    use MaterialPara_type
    implicit none

#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/nmelas_elas.h"
#include "asterfort/nmelas_incr.h"
#include "asterfort/nmorth.h"
!
    type(Behaviour_Integ), intent(in):: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8)                 :: sigp(nsig)
    real(kind=8)                 :: vip(nvi)
    character(len=8), intent(in) :: typmod(2)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out):: codret
!
! --------------------------------------------------------------------------------------------------
!
!   RELATION='ELAS': COMPORTEMENT ELASTIQUE INCREMENTAL (ISOTROPE, ISOTROPE TRANSVERSE, ORTHOTROPE)
!
! --------------------------------------------------------------------------------------------------
!
!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!               TYPMOD  TYPE DE MODELISATION
!               IMATE    ADRESSE DU MATERIAU CODE
!               COMPOR    COMPORTEMENT DE L ELEMENT
!               INSTAM   INSTANT T
!               INSTAP   INSTANT T+DT
!               EPSM   DEFORMATION TOTALE A T
!               DEPS   INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!               OPTION     OPTION DE CALCUL A FAIRE
!               ANGMAS
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1)  :: poum
    aster_logical :: lMatr, lSigm, lVari
    integer(kind=8) :: ndimsi
    real(kind=8) :: sig(2*ndim), dsde(2*ndim, 2*ndim), vi(nvi), zero(2*ndim), eps(2*ndim)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = BEHInteg%materPara
    ASSERT(neps .eq. nsig)
    ASSERT(neps .ge. 2*ndim)

    ndimsi = 2*ndim
    codret = 0
    sig = 0
    vi = 0
    dsde = 0
    zero = 0
    eps = epsm(1:ndimsi)+deps(1:ndimsi)

    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)

    if (lVari) vip = 0

    if (materPara%elasID .eq. ELAS_ISOT) then
        if (compor(INCRELAS) .eq. 'COMP_INCR') then
            call nmelas_incr(BEHInteg, &
                             fami, kpg, ksp, typmod, &
                             jvMaterCode, deps(1:ndimsi), sigm(1:ndimsi), option, &
                             sig, vi, dsde)

        else if (compor(INCRELAS) .eq. 'COMP_ELAS') then
            call nmelas_elas(BEHInteg, &
                             fami, kpg, ksp, typmod, &
                             jvMaterCode, eps, option, &
                             sig, vi, dsde)
        else
            ASSERT(ASTER_FALSE)
        end if

! --------------------------------------------------------------------------------------------------
!  Elasticite isotrope transverse et orthotrope
! --------------------------------------------------------------------------------------------------

    else if (materPara%elasID .eq. ELAS_ORTH .or. &
             materPara%elasID .eq. ELAS_ISTR) then
        if (compor(INCRELAS) .eq. 'COMP_INCR') then
            poum = "T"
            call nmorth(materPara, ndim, poum, &
                        deps, sigm, option, &
                        sig, dsde)
            vi(1) = 0

        else if (compor(INCRELAS) .eq. 'COMP_ELAS') then
            poum = merge('-', '+', option(1:9) .eq. 'RIGI_MECA')
            call nmorth(materPara, ndim, poum, &
                        eps, zero, option, &
                        sig, dsde)
            vi(1) = 0

        else
            ASSERT(ASTER_FALSE)
        end if

        if (lSigm) sigp(1:ndimsi) = sig
        if (lVari) vip(1:nvi) = vi
        if (lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde

    else
        ASSERT(ASTER_FALSE)
    end if
    if (lSigm) sigp(1:ndimsi) = sig
    if (lVari) vip(1:nvi) = vi
    if (lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde

end subroutine
