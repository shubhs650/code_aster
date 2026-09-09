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
#include "asterf_types.h"
!
subroutine pmimpr(prtLevel, &
                  timeCurr, iterNewt, &
                  loadType_, valeImpo_, &
                  epsi_, sigm_, nbVari_, vi_, resi_, &
                  ee_, eini_)
!
    implicit none
!
#include "asterfort/infniv.h"
!
    integer(kind=8), intent(in) :: prtLevel
    real(kind=8), intent(in) :: timeCurr
    integer(kind=8), intent(in) :: iterNewt
    integer(kind=8), optional, intent(in) :: loadType_(6)
    real(kind=8), optional, intent(in) :: valeImpo_(6)
    integer(kind=8), optional, intent(in) :: nbvari_
    real(kind=8), optional, intent(in) :: epsi_(6), sigm_(6), vi_(*), resi_(12)
    real(kind=8), optional, intent(in) :: ee_, eini_
!
! --------------------------------------------------------------------------------------------------
!
! SIMU_POINT_MAT
!
! Print
!
! --------------------------------------------------------------------------------------------------
!
! IN  IND    : 0 pour l'état initial
!                1 pour l'itaration courante
!                2 pour la convergence
!                3 pour l'erreur relative
!                4 pour l'erreur absolue
! IN  INST   : INSTANT ACTUEL
! IN  VALIMP : VALEUR DE LA CMP DE EPSI OU SIGM IMPOSEE
! IN  ITER   : NUMERO D'ITERATION
! IN  EPS    : DEFORMATIONS
! IN  SIG    : CONTRAINTES
! IN  VI     : VARIABLES INTERNES
! IN  NBVARI : Nombre de variables internes
! IN  R      : RESIDU ACTUEL
! IN  EE     : ERREUR
! IN  EINI   : ERREUR INITIALE
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: lDebug = ASTER_FALSE
    character(len=4), parameter :: epsiName(6) = (/'EPXX', 'EPYY', 'EPZZ', &
                                                   'EPXY', 'EPXZ', 'EPYZ'/)
    character(len=4), parameter :: sigmName(6) = (/'SIXX', 'SIYY', 'SIZZ', &
                                                   'SIXY', 'SIXZ', 'SIYZ'/)
    integer(kind=8) :: niv, ifm, i
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)

    if (niv .ge. 2) then
        if (prtLevel .eq. 0) then
            write (ifm, *) ' '
            write (ifm, *) ' ==============================================='
            write (ifm, *) 'INST', timeCurr
            do i = 1, 6
                if (valeImpo_(i) .ne. 0.d0) then
                    if (loadType_(i) .eq. 0) then
                        write (ifm, *) sigmName(i), ' IMPOSEE =', valeImpo_(i)
                    else if (loadType_(i) .eq. 1) then
                        write (ifm, *) epsiName(i), ' IMPOSEE =', valeImpo_(i)
                    end if
                end if
            end do
            if (lDebug) then
                write (ifm, *) ' ETAT INITIAL '
                write (ifm, '(1X,A4,6(1X,E12.5))') 'EPSM', epsi_
                write (ifm, '(1X,A4,6(1X,E12.5))') 'SIGM', sigm_
                write (ifm, '(1X,A4,6(1X,E12.5))') 'VIM', (vi_(i), i=1, min(6, nbVari_))
                if (nbVari_ .gt. 6) then
                    write (ifm, '(1X,A4,6(1X,E12.5))') '   ', (vi_(i), i=7, nbVari_)
                end if
                write (ifm, '(1X,A4,6(1X,E12.5))') 'RESI', resi_
            end if

        else if (prtLevel .eq. 1) then
            if (lDebug) then
                write (ifm, *) '  '
                write (ifm, *) ' ITERATION', iterNewt
                write (ifm, *) ' '
                write (ifm, '(1X,A4,6(1X,E12.5))') 'EPS', epsi_
                write (ifm, '(1X,A4,6(1X,E12.5))') 'SIG', sigm_
                write (ifm, '(1X,A4,6(1X,E12.5))') 'VAR', (vi_(i), i=1, min(6, nbVari_))
                if (nbVari_ .gt. 6) then
                    write (ifm, '(5X,6(1X,E12.5))') (vi_(i), i=7, nbVari_)
                end if
                write (ifm, '(1X,A4,6(1X,E12.5))') 'RESI', resi_
            end if

        else if (prtLevel .eq. 2) then
            if (lDebug) then
                write (ifm, *) '  '
                write (ifm, *) ' ==============================================='
                write (ifm, *) ' CONVERGENCE ITERATION ', iterNewt
                write (ifm, *) ' ==============================================='
                write (ifm, *) ' '
                write (ifm, '(1X,A4,6(1X,E12.5))') 'EPS', epsi_
                write (ifm, '(1X,A4,6(1X,E12.5))') 'SIG', sigm_
                write (ifm, '(1X,A4,6(1X,E12.5))') 'VAR', (vi_(i), i=1, min(6, nbVari_))
                if (nbVari_ .gt. 6) then
                    write (ifm, '(1X,A4,6(1X,E12.5))') '   ', (vi_(i), i=7, nbVari_)
                end if
                write (ifm, '(1X,A4,6(1X,E12.5))') 'RESI', resi_
                write (ifm, *) ' '
                write (ifm, *) ' ==============================================='
            end if

        else if (prtLevel .eq. 3) then
            write (ifm, *) ' -----------------------------------------------'
            write (ifm, '(1X,A4,E12.5,1X,A4,I5,1X,A12,E12.5,1X,A14,E12.5)') &
                'INST', timeCurr, 'ITER', iterNewt, 'ERR.RELATIVE', ee_, 'RESIDU INITIAL', &
                eini_

        else if (prtLevel .eq. 4) then
            write (ifm, *) ' -------------------------------------'
            write (ifm, '(1X,A4,E12.5,1X,A4,I5,1X,A7,E12.5)') 'INST', timeCurr, &
                'ITER', iterNewt, 'ERR_ABS', ee_

        end if
    end if
end subroutine
