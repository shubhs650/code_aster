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
subroutine nmcerr(sddisc, iterGlobMaxi, iterGlobElas, pasMiniElas, resiGlobMaxi, &
                  resiGlobRela, newtKrylResi, ds_contact_)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmecrr.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: iterGlobMaxi, iterGlobElas
    real(kind=8), intent(in) :: pasMiniElas, newtKrylResi
    real(kind=8), intent(in) :: resiGlobMaxi, resiGlobRela
    type(NL_DS_Contact), optional, intent(in) :: ds_contact_
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Discretization management
!
! Management of parameters - Create
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! IN  iterGlobMaxi     : ITER_GLOB_MAXI
! IN  iterGlobElas     : ITER_GLOB_ELAS
! IN  pasMiniElas      : PAS_MINI_ELAS
! IN  newtKrylResi     : PRECISION INITIALE POUR NEWTON-KRYLOV
! IN  resiGlobMaxi     : RESI_GLOB_MAXI
! IN  resiGlobRela     : RESI_GLOB_RELA
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: resiType, nbIter, maxIter, minIter
    integer(kind=8) :: maxIterGeom, nbIterGeom
    integer(kind=8) :: iterSup, nbIterCont
    character(len=24) :: sddiscIfcvName
    real(kind=8), pointer :: sddiscIfcv(:) => null()
    character(len=24) :: sddiscIfreName
    real(kind=8), pointer :: sddiscIfre(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    maxIter = max(iterGlobMaxi, iterGlobElas)
    minIter = min(iterGlobMaxi, iterGlobElas)
    iterSup = 0

! - NOMBRE MAXIMUM D'ITERATIONS
    nbIter = maxIter

! - CREATION DU VECTEUR D'INFORMATIONS SUR LA CONVERGENCE
    sddiscIfcvName = sddisc(1:19)//'.IFCV'
    call wkvect(sddiscIfcvName, 'V V R', 10, vr=sddiscIfcv)

! - Select type of residual
    resiType = 0
    if (resiGlobRela .ne. r8vide()) then
        resiType = resiType+1
    end if
    if (resiGlobMaxi .ne. r8vide()) then
        resiType = resiType+2
    end if
    if (resiType .eq. 0) then
        resiType = 1
    end if

! - Set parameters
    call nmecrr(sddisc, 'MXITER', paraValeI_=maxIter)
    call nmecrr(sddisc, 'MNITER', paraValeI_=minIter)
    call nmecrr(sddisc, 'NBITER', paraValeI_=nbIter)
    call nmecrr(sddisc, 'PAS_MINI_ELAS', paraValeR_=pasMiniElas)
    call nmecrr(sddisc, 'ITERSUP', paraValeI_=iterSup)
    call nmecrr(sddisc, 'RESI_GLOB_RELA', paraValeR_=resiGlobRela)
    call nmecrr(sddisc, 'RESI_GLOB_MAXI', paraValeR_=resiGlobMaxi)
    call nmecrr(sddisc, 'TYPE_RESI', paraValeI_=resiType)
    call nmecrr(sddisc, 'INIT_NEWTON_KRYLOV', paraValeR_=newtKrylResi)
    call nmecrr(sddisc, 'ITER_NEWTON_KRYLOV', paraValeR_=newtKrylResi)

! - RECUPERATION NOMBRE DE REAC_GEOM EN CONTACT DISCRET
    nbIterCont = 0
    if (present(ds_contact_)) then
        if (ds_contact_%l_meca_cont) then
            maxIterGeom = cfdisi(ds_contact_%sdcont_defi, 'ITER_GEOM_MAXI')
            nbIterGeom = cfdisi(ds_contact_%sdcont_defi, 'NB_ITER_GEOM')
            nbIterCont = max(maxIterGeom, nbIterGeom)
        end if
    end if
    nbIter = nbIter+nbIterCont

! - CREATION DU VECTEUR DE STOCKAGE DES RESIDUS
    nbIter = nbIter+1
    sddiscIfreName = sddisc(1:19)//'.IFRE'
    call wkvect(sddiscIfreName, 'V V R', 3*nbIter, vr=sddiscIfre)
!
    call jedema()
end subroutine
