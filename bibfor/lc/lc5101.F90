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

subroutine lc5101(BEHInteg, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, ndsde, &
                  dsidep, codret)

    use Behaviour_type
    use tenseur_dime_module, only: kron, identity, proten
    implicit none

#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/rcvalb.h"

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
!   RELATION='ELAS': COMPORTEMENT ELASTIQUE TOTAL (HYPERELASTIQUE) 1D
! --------------------------------------------------------------------------------------------------
!       in      fami    famille de point de gauss (rigi,mass,...)
!       in      kpg,ksp numero du (sous)point de gauss
!       in      ndim    dimension de l espace (3d=3,2d=2,1d=1)
!               typmod  type de modelisation
!               imate    adresse du materiau code
!               compor    comportement de l element
!               instam   instant t
!               instap   instant t+dt
!               epsm   deformation totale a t
!               deps   increment de deformation totale
!               sigm    contrainte a t
!               vim    variables internes a t-
!               option     option de calcul a faire
!               angmas
!       out     sigp    contrainte a t+dt
!               vip    variables internes a t+dt + indicateur etat t+dt
!               dsidep    matrice de comportement tangent a t+dt ou t
! --------------------------------------------------------------------------------------------------
    integer(kind=8)     :: iokel(1)
    integer(kind=8) :: ndimsi
    real(kind=8) :: valel(1), sig, dsde, vi(nvi), eps
    real(kind=8):: young
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
!  Data preparation
! --------------------------------------------------------------------------------------------------

    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 1)

    codret = 0
    sig = 0
    vi = 0
    dsde = 0

! --------------------------------------------------------------------------------------------------
!  Behaviour integration
! --------------------------------------------------------------------------------------------------

    ! Material parameters
    call rcvalb(fami, kpg, ksp, '+', jvMaterCode, ' ', 'ELAS', 0, ' ', [0.d0], &
                1, ['E'], valel, iokel, 2)
    young = valel(1)

    ! Strain
    eps = epsm(1)+deps(1)

    ! Stress and tangent matrix
    sig = young*eps
    dsde = young

! --------------------------------------------------------------------------------------------------
!  Results preparation
! --------------------------------------------------------------------------------------------------

    if (BEHInteg%behavPara%lSigm) sigp(1:ndimsi) = sig
    if (BEHInteg%behavPara%lVari) vip(1:nvi) = vi
    if (BEHInteg%behavPara%lMatr) dsidep(1:ndimsi, 1:ndimsi) = dsde

end subroutine
