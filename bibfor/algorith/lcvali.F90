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
subroutine lcvali(materPara, &
                  defoComp, epsm, deps, &
                  instam, instap, codret)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: defoComp
    real(kind=8), intent(in) :: deps(:), epsm(:)
    real(kind=8), intent(in) :: instam, instap
    integer(kind=8), intent(inout) :: codret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 4
    character(len=16), parameter :: propName(4) = (/'EPSI_MAXI', 'VEPS_MAXI', &
                                                    'TEMP_MINI', 'TEMP_MAXI'/)
    real(kind=8) :: propVale(4)
    integer(kind=8) :: propCode(4)
    integer(kind=8) :: iret1, iret2, iret3, iret
    real(kind=8) :: epsmax, vepsmax
    real(kind=8) :: dt, tmax, tmin, temp
!
! --------------------------------------------------------------------------------------------------
!
    iret1 = 0
    iret2 = 0
    iret3 = 0
    if (defoComp .ne. 'SIMO_MIEHE') then
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', materPara%jvMaterCode, &
                    ' ', 'VERI_BORNE', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, &
                    propCode, 0)

        if (propCode(1) .eq. 0) then
            epsmax = propVale(1)
            if (norm2(epsm+deps) .gt. epsmax) iret1 = 4
        end if

        if (propCode(2) .eq. 0) then
            vepsmax = propVale(2)
            dt = instap-instam
            if (norm2(deps/dt) .gt. vepsmax) iret2 = 4
        end if

        if (propCode(3) .eq. 0) then
            tmin = propVale(3)
            tmax = propVale(4)
            call rcvarc(' ', 'TEMP', '+', &
                        materPara%schemePara%fami, &
                        materPara%schemePara%kpg, &
                        materPara%schemePara%ksp, &
                        temp, iret)
            if (iret .eq. 0) then
                if ((temp .lt. tmin) .or. (temp .gt. tmax)) then
                    iret3 = 4
                end if
            end if
        end if
        codret = max(iret1, iret2)
        codret = max(codret, iret3)
!
    end if
end subroutine
