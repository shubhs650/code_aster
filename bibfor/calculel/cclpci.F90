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
subroutine cclpci(option, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, &
                  resultIn, resultOut, &
                  ligrel, numeStore, &
                  nbParaIn, lpain, lchin)
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: getOrieField
    implicit none
!
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ
    character(len=8), intent(in) :: resultIn, resultOut
    character(len=24), intent(in) :: ligrel
    integer(kind=8), intent(in) :: numeStore
    integer(kind=8), intent(out) :: nbParaIn
    character(len=8), intent(out) :: lpain(100)
    character(len=24), intent(out) :: lchin(100)
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! Compute ELEM, ELNO and ELGA fields - Set input fields
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION A CALCULER
!   MODELE  K8   NOM DU MODELE
!   RESUIN  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT IN
!   RESUOU  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT OUT
!   MATER   K8   NOM DU MATERIAU
!   CARAEL  K8   NOM DU CARAELE
!   LIGREL  K24  NOM DU LIGREL
!   NUMORD  I    NUMERO D'ORDRE COURANT
!
! OUT :
!   NBPAIN  I    NOMBRE DE PARAMETRES IN
!   LIPAIN  K8*  LISTE DES PARAMETRES IN
!   LICHIN  K8*  LISTE DES CHAMPS IN
!   CODRET  I    CODE RETOUR (0 SI OK, 1 SINON)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: cataInStep, cataInOption
    character(len=8) :: cataInParaName
    character(len=24) :: cataInType, cataInName
    integer(kind=8) :: cataNbIn, iCataIn
    integer(kind=8) :: optionNume, cataInOptionNume, ierd
    integer(kind=8) :: decal, nb_in_maxi
    character(len=8) :: mesh
    character(len=24) :: fieldIn
    integer(kind=8), pointer :: cataDescopt(:) => null()
    character(len=24), pointer :: cataLocalis(:) => null()
    character(len=8), pointer :: cataOptpara(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nbParaIn = 0
    lpain = " "
    lchin = " "

! - Get mesh
    call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)

! - Access to catalog of options
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), optionNume)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', optionNume), 'L', vi=cataDescopt)
    call jeveuo(jexnum('&CATA.OP.LOCALIS', optionNume), 'L', vk24=cataLocalis)
    call jeveuo(jexnum('&CATA.OP.OPTPARA', optionNume), 'L', vk8=cataOptpara)
    cataNbIn = cataDescopt(2)

! - Loop on input parameters for this option
    do iCataIn = 1, cataNbIn
! ----- Get properties of input parameter
        cataInName = cataLocalis(3*iCataIn-1)
        cataInStep = cataLocalis(3*iCataIn) (1:16)
        cataInType = cataLocalis(3*iCataIn-2)
        cataInParaName = cataOptpara(iCataIn) (1:8)

! ----- Set new input parameter
        nbParaIn = nbParaIn+1
        lpain(nbParaIn) = cataInParaName

! ----- Is dependency is an option ?
        cataInOption = cataInName(1:16)
        call jenonu(jexnom('&CATA.OP.NOMOPT', cataInOption), cataInOptionNume)

! ----- Get input field
        fieldIn = " "
        if ((cataInOptionNume .ne. 0) .or. cataInType .eq. 'RESU') then
            if (cataInStep .eq. 'NP1') then
                decal = 1
            else if (cataInStep(1:3) .eq. 'NM1') then
                decal = -1
            else
                decal = 0
            end if
            call rsexch(' ', resultIn, cataInOption, numeStore+decal, fieldIn, ierd)
            if (ierd .ne. 0) then
                call rsexch(' ', resultOut, cataInOption, numeStore+decal, fieldIn, ierd)
            end if
!
            if (ierd .ne. 0) then
                if ((option .eq. cataInOption)) then
!             CAS OU UN CHAMP DEPEND DE LUI MEME A L'INSTANT N-1
!             EXEMPLE : ENDO_ELGA
                    if (cataInStep .eq. 'NM1T') then
                        fieldIn = '&&CALCOP.INT_0'
                        call alchml(ligrel, cataInOption, lpain(nbParaIn), 'V', fieldIn, ierd, ' ')
                        if (ierd .gt. 0) then
                            call utmess('A', 'CALCCHAMP_19', sk=option)
                            cycle
                        end if
                    else
                        call rsexch(' ', resultOut, cataInOption, numeStore+decal, fieldIn, ierd)
                        call alchml(ligrel, cataInOption, lpain(nbParaIn), 'G', fieldIn, ierd, ' ')
                        if (ierd .gt. 0) then
                            call utmess('A', 'CALCCHAMP_19', sk=option)
                            cycle
                        end if
                        call rsnoch(resultOut, cataInOption, numeStore+decal)
                    end if
                else
                    fieldIn = ' '
                end if
            end if

        else if (cataInType .eq. 'MAIL') then
            fieldIn = mesh(1:8)//cataInName(1:16)

        else if (cataInType .eq. 'MODL') then
            fieldIn = modelZ(1:8)//cataInName(1:16)

        else if (cataInType .eq. 'CARA') then
            fieldIn = caraElemZ(1:8)//cataInName(1:16)

        else if (cataInType .eq. 'ORIE') then
            call getOrieField(caraElemZ, cataInName, fieldIn)

        else if (cataInType .eq. 'VOLA') then
            fieldIn = cataInName

        else if (cataInType .eq. 'CHMA') then
            fieldIn = materFieldZ(1:8)//cataInName(1:16)

        else if (cataInType .eq. 'MACO') then
            fieldIn = materCodeZ(1:8)//cataInName(1:16)

        end if

        lchin(nbParaIn) = fieldIn

        !WRITE (6, *) "Ajout du champ d'entrée <", nbParaIn, ">: ", &
        !    lpain(nbParaIn), "/", lchin(nbParaIn)

    end do

! - Add HHO fields
    nb_in_maxi = nbParaIn+3
    ASSERT(nb_in_maxi .le. 100)
    call hhoAddInputField(modelZ, nb_in_maxi, lchin, lpain, nbParaIn)
!
    call jedema()
!
end subroutine
