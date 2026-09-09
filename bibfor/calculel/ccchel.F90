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
subroutine ccchel(option, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, listLoadZ, &
                  resultIn, resultOut, resultType, &
                  numeStore, numeStorePrev, &
                  ligrel, isTransient, postCompPoux, jvBase, &
                  fieldNameOut)
!
    use postComp_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ccaccl.h"
#include "asterfort/cclpci.h"
#include "asterfort/cclpco.h"
#include "asterfort/ccpara.h"
#include "asterfort/ccpoux.h"
#include "asterfort/detrsd.h"
#include "asterfort/meceuc.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ, listLoadZ
    character(len=8), intent(in) :: resultIn, resultOut
    character(len=16), intent(in) ::resultType
    integer(kind=8), intent(in) :: numeStore, numeStorePrev
    character(len=24), intent(in) :: ligrel
    aster_logical, intent(in) :: isTransient
    type(POST_COMP_POUX), intent(in) :: postCompPoux
    character(len=1), intent(in) :: jvBase
    character(len=24), intent(out) :: fieldNameOut
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! Compute ELEM, ELNO and ELGA fields
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION
!   MODELE  K8   NOM DU MODELE
!   RESUIN  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT IN
!   RESUOU  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT OUT
!   NUMORD  I    NUMERO D'ORDRE COURANT
!   NORDM1  I    NUMERO D'ORDRE PRECEDENT
!   MATER   K8   NOM DU MATERIAU
!   MATECO  K8   NOM DU MATERIAU CODE
!   CARAEL  K8   NOM DU CARAELE
!   TYPESD  K16  TYPE DE LA STRUCTURE DE DONNEES RESULTAT
!   LIGREL  K24  NOM DU LIGREL
!   EXIPOU  BOOL EXISTENCE OU NON DE POUTRES POUX
!   EXITIM  BOOL EXISTENCE OU NON DE L'INSTANT DANS LA SD RESULTAT
!   LISCHA  K19  NOM DE L'OBJET JEVEUX CONTENANT LES CHARGES
!   NBCHRE  I    NOMBRE DE CHARGES REPARTIES (POUTRES)
!   IOCCUR  I    NUMERO D'OCCURENCE OU SE TROUVE LE CHARGE REPARTIE
!   SUROPT  K24
!   jvBase  K1   BASE SUR LAQUELLE DOIT ETRE CREE LE CHAMP DE SORTIE
!
! OUT :
!   RESOUT  K24  NOM DU CHAMP OUT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, nbParaOut, nbParaIn
    character(len=8) :: lpain(100), lpaout(1)
    character(len=24) :: lchin(100), lchout(1)
!
! --------------------------------------------------------------------------------------------------
!
    fieldNameOut = ' '

! - Create generic input fields
    call ccpara(option, &
                modelZ, materFieldZ, caraElemZ, &
                resultIn, resultOut, &
                numeStore, numeStorePrev, isTransient)

! - Construct list of input fields
    call cclpci(option, &
                modelZ, materFieldZ, materCodeZ, caraElemZ, &
                resultIn, resultOut, &
                ligrel, numeStore, &
                nbParaIn, lpain, lchin)

! - Construct list of output fields
    call cclpco(option, &
                resultOut, numeStore, &
                nbParaOut, lpaout, lchout)
    ASSERT(nbParaOut .eq. 1)

! - Special for POUX beams
    if (postCompPoux%lPoux) then
        call ccpoux(postCompPoux, &
                    listLoadZ, modelZ, &
                    resultIn, resultType, numeStore, &
                    nbParaIn, lpain, lchin, &
                    iret)
        if (iret .ne. 0) then
            goto 999
        end if
    end if

! - Change input and output parameters for special cases
    call ccaccl(option, &
                modelZ, materFieldZ, caraElemZ, ligrel, &
                resultType, &
                nbParaIn, lpain, lchin, lchout, &
                iret)
    if (iret .ne. 0) then
        goto 999
    end if

! - Compute option with complex case
    call meceuc('C', option, caraElemZ, ligrel, &
                nbParaIn, lchin, lpain, nbParaOut, lchout, &
                lpaout, jvBase)
    fieldNameOut = lchout(1)

! - Clean
    call detrsd('CHAM_ELEM', '&&CALCOP.INT_0')
!
999 continue
!
end subroutine
