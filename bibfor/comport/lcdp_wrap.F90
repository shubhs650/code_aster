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
#include "asterfort/Behaviour_type.h"
!
subroutine lcdp_wrap(fami, kpg, ksp, ndim, imate, &
                     carcri, neps, epsm, &
                     deps, vim, option, sigp, &
                     vip, dsidep, codret)
!
    use lcdp_module, only: dp_material
    implicit none
!
#include "asterfort/lcdp_compute.h"
#include "asterfort/lcdp_material.h"
!
    integer(kind=8) :: imate, ndim, kpg, ksp, codret, neps
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: epsm(neps), deps(neps)
    real(kind=8) :: sigp(neps)
    real(kind=8) :: vim(*), vip(*)
    real(kind=8) :: dsidep(neps, neps)
    character(len=16) :: option
    character(len=*) :: fami
! ----------------------------------------------------------------------
    real(kind=8), parameter, dimension(6)::rac2 = [1.d0, 1.d0, 1.d0, &
                                                   sqrt(2.d0), sqrt(2.d0), sqrt(2.d0)]
! ----------------------------------------------------------------------
    type(dp_material) :: mat
    aster_logical :: elas, rigi, resi
    integer(kind=8) :: itemax, ndimsi, state
    real(kind=8) :: prec
    real(kind=8) :: ka
    real(kind=8) :: eps(2*ndim), ep(2*ndim), sig(2*ndim)
    real(kind=8) :: deps_sig(2*ndim, 2*ndim)
!
! INITIALISATION
!
    elas = option(11:14) .eq. 'ELAS'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    resi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH'
    itemax = nint(carcri(1))
    prec = carcri(3)
    ndimsi = 2*ndim
!
! EXTRACTION DES DONNEES CINEMATIQUES
    eps = epsm(1:ndimsi)
    if (resi) eps = eps+deps(1:ndimsi)
!
! LECTURE DU MATERIAU
!
    mat = lcdp_material(fami, kpg, ksp, imate, resi)
!
! LECTURE DES VARIABLES INTERNES
!
    ka = vim(1)
    state = nint(vim(3))
    ep = vim(4:3+ndimsi)*rac2(1:ndimsi)
!
! COMPORTEMENT
    codret = lcdp_compute( &
             resi, rigi, elas, itemax, prec, mat, eps, ep, ka, state, sig, deps_sig, vip)
    if (codret .ne. 0) goto 999
!
! MISE A JOUR DES CHAMPS
!
    if (resi) then
        vip(1) = ka
        vip(2) = sum(ep(1:3))
        vip(3) = state
        vip(4:3+ndimsi) = ep/rac2(1:ndimsi)
    end if
!
    if (resi) sigp(1:ndimsi) = sig
    if (rigi) dsidep(1:ndimsi, 1:ndimsi) = deps_sig
!
999 continue
end subroutine
