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
subroutine lcspelas(fami, kpg, ksp, ndim, &
                    mate, nomat, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, ndsde, dsidep, codret, BEHinteg_)

! aslint: disable=I1306
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/rcvalb.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    integer(kind=8) :: mate, ndim, neps, nsig, nvi, ndsde, kpg, ksp, codret
    real(kind=8) :: epsm(neps), deps(neps)
    real(kind=8) :: sigm(nsig), sigp(nsig), dsidep(nsig, neps)
    real(kind=8) :: vim(nvi), vip(nvi), instam, instap
    character(len=8) :: nomat
    real(kind=8) :: carcri(CARCRI_SIZE)
    character(len=16) :: option
    character(len=*) :: fami
    type(Behaviour_Integ), optional:: BEHinteg_
!
!-----------------------------------------------------------------------
!     LOI DE COMPORTEMENT  D'INTERFACE
!     POUR LES ELEMENTS D'INTERFACE BIPHASIQUE.
!
! IN : EPSM DÉPLACEMENT RELATIF INSTANT MOINS
! IN : DEPS INCR DE DÉPLACEMENT RELATIF
! IN : SIGM EFFORT RELATIF INSTANT MOINS
! IN : MATE, OPTION, VIM, COOROT,INSTAM, INSTAP
! OUT : SIGP , DSIDEP , VIP
!
! m : tout ce qui se rapporte à l'INSTANT précédent
!   (ie l'état convergé précédent)
! p : itération actuelle (ie l'état non convergé)
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nbpael
    parameter(nbpael=2)
    integer(kind=8) :: codel(nbpael)
    real(kind=8) :: kn, kt, valel(nbpael)
    real(kind=8) :: epsp(neps)
    real(kind=8) :: cel(nsig)
    character(len=16) :: nomel(nbpael)
    character(len=1) :: poum
    integer(kind=8) :: nb_para
    character(len=8) :: para_name(5)
    real(kind=8) :: para_vale(5)
    aster_logical :: resi, rigi
    blas_int :: b_incx, b_incy, b_n
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
! CALCUL DE CONTRAINTE (RESIDU)
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
! CALCUL DE LA MATRICE TANGEANTE (RIGIDITE)
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
!
! #####################################
! RECUPERATION DES PARAMETRES PHYSIQUES
! #####################################
    nomel(1) = 'K_N'
    nomel(2) = 'K_T'
!
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if

    nb_para = 0
    para_name = ' '
    para_vale = 0.d0
!   évolution des paramètres
    if (present(BEHinteg_)) then
        nb_para = nb_para+1
        para_name(nb_para) = 'X'
        para_vale(nb_para) = BEHinteg_%behavESVA%behavESVAGeom%coorElga(kpg, 1)
        nb_para = nb_para+1
        para_name(nb_para) = 'Y'
        para_vale(nb_para) = BEHinteg_%behavESVA%behavESVAGeom%coorElga(kpg, 2)
        nb_para = nb_para+1
        para_name(nb_para) = 'Z'
        para_vale(nb_para) = BEHinteg_%behavESVA%behavESVAGeom%coorElga(kpg, 3)
    end if
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                nomat, 'INTERF_POU_ELAS', nb_para, para_name, [para_vale], &
                nbpael, nomel, valel, codel, 2)
! DEFINITION DE PARAMETRES PHYSIQUE:
    kn = valel(1)
    kt = valel(2)
!
! #####################################
! INITIALISATION DE VARIABLES
! #####################################
!
    b_n = to_blas_int(neps)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, epsp, b_incy)
    if (resi) call daxpy(b_n, 1.d0, deps, b_incx, epsp, b_incy)
!
! #####################################
! CALCUL DE LA CONTRAINTE ET DE LA RIGIDITE ELASTIQUE
! #####################################
!
    sigp(:) = 0.d0
    dsidep(:, :) = 0.d0
    if (rigi) then
        dsidep = reshape((/kt, 0.d0, 0.d0, 0.d0, kn, 0.d0, 0.d0, 0.d0, kn/), (/3, 3/))
    end if
    if (resi) then
        cel = (/kt, kn, kn/)
        sigp = cel*epsp
    end if
!
end subroutine
