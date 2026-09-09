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

subroutine ccaccl(option, &
                  modelZ, materFieldZ, caraElemZ, ligrel, &
                  resultType, &
                  nbParaIn, lpain, lchin, lchout, &
                  codret)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cesvar.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: modelZ, materFieldZ, caraElemZ
    character(len=24), intent(in) :: ligrel
    character(len=16), intent(in) :: resultType
    integer(kind=8), intent(in) :: nbParaIn
    character(len=8), intent(in) :: lpain(100)
    character(len=24), intent(inout) :: lchin(100)
    character(len=24), intent(inout) :: lchout(1)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! Compute ELEM, ELNO and ELGA fields - Input and output fields for special cases
!
! --------------------------------------------------------------------------------------------------
!
!  ROUTINE DE GESTION DES GLUTES NECESSAIRES POUR CERTAINES OPTIONS
!  * POUTRE POUX, DCHA_*, RADI_*, ...
!  * APPEL A CESVAR POUR LES CHAMPS A SOUS-POINTS.
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION
!   MODELE  K8   NOM DU MODELE
!   RESUIN  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT IN
!   MATER   K8   NOM DU MATERIAU
!   CARAEL  K8   NOM DU CARAELE
!   LIGREL  K24  NOM DU LIGREL
!   NUMORD  I    NUMERO D'ORDRE COURANT
!   NORDM1  I    NUMERO D'ORDRE PRECEDENT
!   TYPESD  K16  TYPE DE LA STRUCTURE DE DONNEES RESULTAT
!   NBPAIN  I    NOMBRE DE PARAMETRES IN
!   LIPAIN  K8*  LISTE DES PARAMETRES IN
!   LICHOU  K24* LISTE DES CHAMPS OUT
!
! IN/OUT :
!   LICHIN  K24* LISTE MODIFIEE DES CHAMPS IN
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: canbva = '&&CCACCL.CANBVA'
    character(len=24), parameter :: chnlin = '&&CCACCL.PNONLIN'
    integer(kind=8) :: iret, paraIndx, inume, nbsp
    character(len=8) :: mesh, caraElemToApply, parain
    character(len=19) :: compor, comporToApply
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    codret = 0

! - Get mesh
    call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)

! - Special cases
    if (option .eq. 'EFGE_ELNO') then
! ----- EFGE_ELNO: create map for linear or non-linear case
        inume = 0
        if (resultType .eq. 'EVOL_NOLI') then
            inume = 1
        end if
        call mecact('V', chnlin, 'MAILLA', mesh, 'NEUT_I', &
                    ncmp=1, nomcmp='X1', si=inume)
        if (inume .eq. 0) then
            paraIndx = indik8(lpain, 'PCOMPOR', 1, nbParaIn)
            ASSERT(paraIndx .ge. 1)
            lchin(paraIndx) = materFieldZ(1:8)//'.COMPOR'
        end if
!
    else if (option .eq. 'VARI_ELNO') then
! ----- Behaviour map for dimensionning number of internal state variables
        paraIndx = indik8(lpain, 'PCOMPOR', 1, nbParaIn)
        compor = " "
        if (paraIndx .ne. 0) then
            compor = lchin(paraIndx) (1:19)
        end if
        call exisd('CARTE', compor, iret)
        if (iret .ne. 1) then
            call utmess('A', 'CALCULEL2_86')
            codret = 1
            goto 30
        end if
    end if

!     ---------------------------------------------------------------
!     -- AJOUT EVENTUEL DU CHAM_ELEM_S PERMETTANT LES SOUS-POINTS
!        ET LE BON NOMBRE DE VARIABLES INTERNES
!     ---------------------------------------------------------------
    if ((option .eq. 'EPEQ_ELGA') .or. (option .eq. 'EPEQ_ELNO') .or. &
        (option .eq. 'EPSI_ELGA') .or. (option .eq. 'EPSI_ELNO') .or. &
        (option .eq. 'SIEF_ELGA') .or. (option .eq. 'SIEF_ELNO') .or. &
        (option .eq. 'SIEQ_ELGA') .or. (option .eq. 'SIEQ_ELNO') .or. &
        (option .eq. 'SIGM_ELGA') .or. (option .eq. 'SIGM_ELNO') .or. &
        (option .eq. 'EPVC_ELGA') .or. (option .eq. 'EPVC_ELNO') .or. &
        (option .eq. 'EPME_ELGA') .or. (option .eq. 'EPME_ELNO') .or. &
        (option .eq. 'EPSP_ELGA') .or. (option .eq. 'EPSP_ELNO') .or. &
        (option .eq. 'VARI_ELNO') .or. (option .eq. 'DEPL_ELGA') .or. &
        (option .eq. 'TEMP_ELGA') .or. (option .eq. 'VARC_ELGA') .or. &
        (option .eq. 'VARC_ELNO')) then

! ----- Internal state variables
        if (option .eq. 'VARI_ELNO') then
            comporToApply = compor
        else
            comporToApply = ' '
        end if

! ----- For "sub-points"
        caraElemToApply = caraElemZ

! ----- POUR LES OPTIONS SUIVANTES, LE NOMBRE DE SOUS-POINTS
! ----- DU CHAMP "OUT" DEPEND D'UN CHAMP "IN" PARTICULIER :
        if (option .eq. 'EPEQ_ELGA') then
            parain = 'PDEFORR'
        else if (option .eq. 'EPEQ_ELNO') then
            parain = 'PDEFORR'
        else if (option .eq. 'EPSI_ELNO') then
            parain = 'PDEFOPG'
        else if (option .eq. 'SIEF_ELNO') then
            parain = 'PCONTRR'
        else if (option .eq. 'SIEQ_ELGA') then
            parain = 'PCONTRR'
        else if (option .eq. 'SIEQ_ELNO') then
            parain = 'PCONTRR'
        else if (option .eq. 'SIGM_ELNO') then
            parain = 'PCONTRR'
        else if (option .eq. 'SIGM_ELGA') then
            parain = 'PSIEFR'
        else if (option .eq. 'VARC_ELGA') then
            parain = 'PVARCPR'
        else if (option .eq. 'VARC_ELNO') then
            parain = 'PVARCGR'
        else
            parain = ' '
        end if
!
        if (parain .ne. ' ') then
            paraIndx = indik8(lpain, parain, 1, nbParaIn)
            ASSERT(paraIndx .ge. 1)
            call jeexin(lchin(paraIndx) (1:19)//'.CELD', iret)
            nbsp = 1
            if (iret .ne. 0) then
                call dismoi('MXNBSP', lchin(paraIndx), 'CHAM_ELEM', repi=nbsp)
            end if
            if (nbsp .le. 1) then
                caraElemToApply = ' '
            end if
        end if

! ----- Create object for output field
        if (caraElemToApply .ne. ' ' .or. comporToApply .ne. ' ') then
            call cesvar(caraElemToApply, comporToApply, ligrel, canbva)
            call copisd('CHAM_ELEM_S', 'V', canbva, lchout(1))
            call detrsd('CHAM_ELEM_S', canbva)
        end if
    end if
!
30  continue
    call jedema()
!
end subroutine
