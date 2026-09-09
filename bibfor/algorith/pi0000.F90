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
subroutine pi0000(BEHInteg, compor, typmod, ndim, &
                  epsm, epsd_cste, epsd_pilo, &
                  sigm, vim, dtau, etamin, etamax, copilo)

    use Behaviour_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/pi0002.h"
#include "asterfort/pi0007.h"
#include "asterfort/pi0036.h"
#include "asterfort/pi0060.h"
#include "asterfort/pi0076.h"
#include "asterfort/pi0120.h"
#include "asterfort/pi5076.h"
#include "asterfort/pi6036.h"
#include "asterfort/pi6046.h"
#include "asterfort/pi6057.h"
#include "asterfort/pi7010.h"
#include "asterfort/pi7011.h"
#include "asterfort/pi7046.h"
#include "asterfort/pi9040.h"
#include "asterfort/pi9041.h"
#include "asterfort/pi9051.h"
#include "asterfort/pi9056.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=8), intent(in) :: typmod(2)
    integer(kind=8), intent(in):: ndim
    real(kind=8), intent(in) :: epsm(:)
    real(kind=8), intent(in) :: epsd_cste(:)
    real(kind=8), intent(in) :: epsd_pilo(:)
    real(kind=8), intent(in) :: sigm(:)
    real(kind=8), intent(in) :: vim(:)
    real(kind=8), intent(in) :: dtau
    real(kind=8), intent(in) :: etamin
    real(kind=8), intent(in) :: etamax
    real(kind=8), intent(out) :: copilo(:)
! --------------------------------------------------------------------------------------------------
!
!     PILOTAGE PRED_ELAS : BRANCHEMENT SELON COMPORTEMENT
!
! --------------------------------------------------------------------------------------------------
! in  neps    dimension des deformations
! in  dtau    increment de pilotage
! in  vim     variables internes en t-                       (pred_elas)
! in  sigm    contraintes en t- (si necessaire)              (pred_elas)
! in  epsm    champ de deformation en t-
! in  epsd_cste    increment fixe
! in  epsd_pilo    increment pilote
! in  etamin  borne inf du pilotage (si utile)               (pred_elas)
! in  etamax  borne sup du pilotage (si utile)               (pred_elas)
! out copilo  coefficient de pilotage : f := a0+a1*eta = dtau
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: numlc, nvi
    character(len=16):: relaComp
! --------------------------------------------------------------------------------------------------
!
    ! Behaviour index
    numlc = BEHInteg%behavPara%numlc+BEHInteg%behavPara%lawIndexOffset

    ! Number of internal variables
    nvi = BEHInteg%behavPara%nvi
    ASSERT(.not. BEHInteg%behavPara%lGdefLog)
    ASSERT(.not. BEHInteg%behavPara%lReguVisc)
    ASSERT(.not. BEHinteg%behavPara%lAnnealing)

    ! Constitutive law name
    relaComp = compor(RELA_NAME)

    select case (numlc)

    case (2)
        call pi0002(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo, relaComp)

    case (7)
        call pi0007(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (36)
        call pi0036(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (60)
        call pi0060(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (76)
        call pi0076(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (120)
        call pi0120(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (5076)
        call pi5076(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (6036)
        call pi6036(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (6046)
        call pi6046(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (6057)
        call pi6057(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (7010)
        call pi7010(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (7011)
        call pi7011(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (7046)
        call pi7046(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (9040)
        call pi9040(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (9041)
        call pi9041(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (9051)
        call pi9051(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case (9056)
        call pi9056(BEHInteg, typmod, ndim, epsm, epsd_cste, epsd_pilo, &
                    sigm, vim(1:nvi), dtau, etamin, etamax, copilo)

    case default
        call utmess('F', 'MECANONLINE_59')
    end select

end subroutine

! if (typmod(2) .eq. 'GRADVARI') then
!     if (typilo .eq. 'DEFORMATION') then
!         call pidegv(neps, dtau, epsm, epsd_cste, epsd_pilo, copilo)

!     else
!         if (etamin .eq. -r8gaem() .or. etamax .eq. r8gaem()) &
!             call utmess('F', 'MECANONLINE_60', sk=relaComp)

!         if (relaComp .eq. 'ENDO_SCALAIRE') then
!             call piesgv(neps, dtau, jvMaterCode, lcquma, vim, epsm, &
!                         epsd_cste, epsd_pilo, typmod, &
!                         lcquga, etamin, etamax, lcqubo, copilo)

!         else if (relaComp .eq. 'ENDO_FISS_EXP') then
!             call piesgv(neps, dtau, jvMaterCode, lcmfma, vim, epsm, &
!                         epsd_cste, epsd_pilo, typmod, &
!                         lcmfga, etamin, etamax, lcmfbo, copilo)

!         else if (relaComp .eq. 'ENDO_ISOT_BETON') then
!             call pieigv(neps, dtau, jvMaterCode, vim, epsm, &
!                         epsd_cste, epsd_pilo, typmod, &
!                         etamin, etamax, copilo)

!         else
!             call utmess('F', 'MECANONLINE_59')
!         end if
!     end if

! else if (typmod(2) .eq. 'INTERFAC') then

!     ASSERT(typilo .ne. 'DEFORMATION')
!     ndim = neps/2
!     su_cste = 0
!     su_pilo = 0
!     mu_cste = 0
!     mu_pilo = 0
!     su_cste(1:ndim) = epsm(1:ndim)+epsd_cste(1:ndim)
!     su_pilo(1:ndim) = epsd_pilo(1:ndim)
!     mu_cste(1:ndim) = epsm(ndim+1:2*ndim)+epsd_cste(ndim+1:2*ndim)
!     mu_pilo(1:ndim) = epsd_pilo(ndim+1:2*ndim)

!     if (relaComp .eq. 'CZM_TAC_MIX') then
!         call pipetc(jvMaterCode, su_cste, su_pilo, mu_cste, mu_pilo, &
!                     vim, dtau, copilo)
!     else if (relaComp .eq. 'CZM_OUV_MIX') then
!         call pipeou(jvMaterCode, su_cste, su_pilo, mu_cste, mu_pilo, &
!                     vim, dtau, copilo)
!     else if (relaComp .eq. 'CZM_EXP_MIX') then
!         call pipeex(jvMaterCode, &
!                     su_cste, su_pilo, mu_cste, mu_pilo, &
!                     vim, dtau, copilo)
!     else if (relaComp .eq. 'CZM_LAB_MIX') then
!         call pipeab(jvMaterCode, dtau, vim(:), &
!                     su_cste, su_pilo, mu_cste, mu_pilo, nsol, sol, sgn)

!         if (nsol .eq. 0) then
!             copilo(5) = 0.d0
!         else if (nsol .eq. 1) then
!             copilo(1) = dtau-sgn(1)*sol(1)
!             copilo(2) = sgn(1)
!         else if (nsol .eq. 2) then
!             copilo(1) = dtau-sgn(1)*sol(1)
!             copilo(2) = sgn(1)
!             copilo(3) = dtau-sgn(2)*sol(2)
!             copilo(4) = sgn(2)
!         end if

!     else
!         call utmess('F', 'MECANONLINE_59')
!     end if

! else
!     call utmess('F', 'MECANONLINE_61', sk=typmod(2))
! end if
