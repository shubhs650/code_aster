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
! aslint: disable=W1306
!
subroutine ngpipe(BEHInteg, typmod, compor, ndim, npg, neps, nddl, b, ni2ldc, &
                  ddlm, ddld, ddl0, ddl1, lgpg, sigm_user, vim, &
                  dtau, etamin, etamax, copilo)

    use Behaviour_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/pi0000.h"
#include "blas/dgemv.h"
!
    type(Behaviour_Integ) :: BEHInteg
    character(len=8) :: typmod(2)
    character(len=16) ::compor(COMPOR_SIZE)
    integer(kind=8), intent(in):: ndim
    integer(kind=8), intent(in):: npg
    integer(kind=8), intent(in):: neps
    integer(kind=8), intent(in):: nddl
    real(kind=8), intent(in) :: b(neps, npg, nddl)
    real(kind=8), intent(in) :: ni2ldc(neps, npg)
    real(kind=8), intent(in) :: ddlm(nddl)
    real(kind=8), intent(in) :: ddld(nddl)
    real(kind=8), intent(in) :: ddl0(nddl)
    real(kind=8), intent(in) :: ddl1(nddl)
    integer(kind=8), intent(in) :: lgpg
    real(kind=8), intent(in) :: sigm_user(neps, npg)
    real(kind=8), intent(in) :: vim(lgpg, npg)
    real(kind=8), intent(in) :: dtau
    real(kind=8), intent(in) :: etamin
    real(kind=8), intent(in) :: etamax
    real(kind=8), intent(out) :: copilo(5, npg)
! --------------------------------------------------------------------------------------------------
!     BUT:  CALCUL  DES COEFFICIENTS DE PILOTAGE POUR PRED_ELAS
! --------------------------------------------------------------------------------------------------
! in  BEHInteg: descripteur de la loi de comportement
! in  typmod : type de modelisation
! in  compor : carte comportement
! in  ndim   : dimension de l'espace
! in  npg    : nombre de points de gauss
! in  neps   : nombre de composantes de deformations / contraintes
! in  nddl   : nombre de ddl dans l'element
! in  b      : matrice cinematique
! in  ni2ldc : conversion contrainte --> avec racine de deux
! in  lgpg   : "longueur" des variables internes pour 1 point de gauss
!             cette longueur est un majorant du nbre reel de var. int.
! in  ddlm   : ddl u,alpha,mu en t-
! in  sigm   : contraintes de cauchy en t- (inutile)
! in  vim    : variables internes en t-
! in  ddld   : increment de ddl u,alpha,mu a l'iteration newton courante
! in  ddl0   : correction de ddl u,alpha,mu pour forces fixes
! in  ddl1   : correction de ddl u,alpha,mu pour forces pilotees
! in  dtau   : incrément de pilotage
! out copilo : coefficients a0 et a1 pour chaque point de gauss
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: kpg, nepg
    real(kind=8) :: sigm_ldc(neps, npg)
    real(kind=8) :: epsm(neps, npg), epsd_pilo(neps, npg), epsd_cste(neps, npg)
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
! --------------------------------------------------------------------------------------------------

    ASSERT(compor(DEFO) .eq. 'PETIT')
    copilo = r8vide()
    nepg = neps*npg

    ! Deformations
    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddlm, b_incx, 0.d0, epsm, &
               b_incy)

    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddld+ddl0, b_incx, 0.d0, epsd_cste, &
               b_incy)

    b_lda = to_blas_int(nepg)
    b_m = to_blas_int(nepg)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('N', b_m, b_n, 1.d0, b, &
               b_lda, ddl1, b_incx, 0.d0, epsd_pilo, &
               b_incy)

! -- Contrainte au format loi de comportement
    sigm_ldc = sigm_user*ni2ldc

! - TRAITEMENT DE CHAQUE POINT DE GAUSS
    do kpg = 1, npg
        call pi0000(BEHInteg, compor, typmod, ndim, &
                    epsm(:, kpg), epsd_cste(:, kpg), epsd_pilo(:, kpg), &
                    sigm_ldc(:, kpg), vim(:, kpg), dtau, etamin, etamax, &
                    copilo(:, kpg))
    end do
!
end subroutine
