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
subroutine calcop(option, listOptJvZ, resultIn, resultOut, listStoreJv, &
                  nbStore, resultType, codret, jvBase_, tldist_)
!
    use postComp_type
    use postComp_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/ccchel.h"
#include "asterfort/ccchno.h"
#include "asterfort/ccliop.h"
#include "asterfort/cclodr.h"
#include "asterfort/cclord.h"
#include "asterfort/ccnett.h"
#include "asterfort/ccvepo.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/getelem.h"
#include "asterfort/getvid.h"
#include "asterfort/indk16.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/medom1.h"
#include "asterfort/pcptcc.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rs_get_model.h"
#include "asterfort/rsexc1.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/xthpos.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: listOptJvZ
    character(len=8), intent(in) :: resultIn, resultOut
    character(len=19), intent(in) :: listStoreJv
    integer(kind=8), intent(in) :: nbStore
    character(len=16), intent(in) :: resultType
    integer(kind=8), intent(out) ::  codret
    character(len=1), optional, intent(in) :: jvBase_
    aster_logical, optional, intent(in)  :: tldist_
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP
!
!  Main subroutine to compute option
!
! --------------------------------------------------------------------------------------------------
!
! In  option            : option to compute
! In  listOptionJv      : JEVEUX name object for list of options to compute
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
! In  nbStore           : number of storing indexes
! In  listStoreJv       : JEVEUX name object for list of storing indexes
! In  resultType        : type of results datastructure
! Out codret            : error code
!                           0 - Everything is OK
!                           1 - Something goes wrong
! In  jvBase            : jeveux base to save fields
! In  tlDist            : time distribution
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: jvBaseName = '&&CALCOP'
    character(len=24), parameter :: restCellJv = '&&OP0106.MES_MAILLES'
    aster_logical :: isTransient, isOptionFromUser, dbg_ob, dbgv_ob, lcpu, ltest, ldist
    aster_logical :: ligrelHasBeenChanged, lbid, lsdpar, l_pmesh
    aster_logical :: lRestCell
    mpi_int :: mpicou, mpibid
    integer(kind=8) :: nbOptEff, iOptEff, ibid, posopt, jvcham
    integer(kind=8) :: numeStoreMin, numeStoreMax, jlinst, iStoreToCompute, nbStoreToCompute
    integer(kind=8) :: rang, nbproc
    integer(kind=8) :: numeStore, iret, nbRestCell, codre2, nbOpt, ipas, nbpas, jldist
    integer(kind=8) :: numeStorePrev, nbLoad, p, k, numork
    integer(kind=8) :: ideb, ifin, irelat, ifm, niv, nbEqua, nbEquaNew
    real(kind=8), pointer :: noch(:) => null()
    real(kind=8), pointer :: prbid(:) => null()
    complex(kind=8), pointer :: nochc(:) => null()
    complex(kind=8), pointer :: pcbid(:) => null()
    integer(kind=8) :: nbStoreIn, nocc
    integer(kind=8), pointer :: listStoreIn(:) => null(), listStore(:) => null()
    integer(kind=8), pointer :: lacalc(:) => null()
    character(len=1) :: jvBase, kbid, ktyp
    character(len=5) :: numeOptStr
    character(len=8) :: model, modelNew, modelCmd, modelRefe
    character(len=8) :: caraElem, mesh
    character(len=16) :: optionEff
    character(len=19) :: listLoad, k19b, nochou, nochok, partsd
    character(len=24) :: materField, materCode
    character(len=24) :: fieldNameOut, ligrel, ligrelSave, vldist, vcham, vcnoch
    character(len=24), parameter :: k24b = " "
    character(len=24) :: listOptEffJv, listStoreOptJv
    character(len=24) :: chamno
    character(len=24) :: listOptJv
    character(len=16), pointer :: listOpt(:) => null(), listOptEff(:) => null()
    type(POST_COMP_REST) :: postCompRest
    type(POST_COMP_POUX) :: postCompPoux
    type(POST_COMP_RESU) :: postCompResu
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

! - Initializations
    codret = 1
    listOptJv = listOptJvZ
    listLoad = '&&CALCOP.LISCHA'
    nbLoad = 0

! - On reporte ici un post-traitement XFEM depuis OP0025
    if (option .eq. 'TEMP_ELGA') then
        call dismoi('NOM_MODELE', resultIn, 'RESULTAT', repk=model)
        call exixfe(model, iret)
        if (iret .ne. 0) then
            call xthpos(resultIn, resultOut)
            codret = 0
            goto 999
        end if
    end if

! - Access to list of storing index to compute options
    if (nbStore .lt. 1) then
        call utmess("F", "CALCCHAMP1_1")
    end if
    call jeveuo(listStoreJv, 'L', vi=listStore)

! - Access to list of options
    nbOpt = 0
    if (listOptJv .ne. ' ') then
        call jeveuo(listOptJv, 'E', vk16=listOpt)
        call jelira(listOptJv, 'LONMAX', nbOpt)
    end if

! - Update list of options with dependencies => effective list of options
    call ccliop(option, jvBaseName, listOptEffJv, nbOptEff)
    if (nbOptEff .eq. 0) goto 999
    call jeveuo(listOptEffJv, 'L', vk16=listOptEff)

! - Transient computation
    call jenonu(jexnom(resultIn//'           .NOVA', 'INST'), iret)
    isTransient = iret .ne. 0

! - Get minimum and maximum storage index in input result
    call rs_get_liststore(resultIn, nbStoreIn)
    ASSERT(nbStoreIn .gt. 0)
    AS_ALLOCATE(vi=listStoreIn, size=nbStoreIn)
    call rs_get_liststore(resultIn, nbStoreIn, listStoreIn)
    numeStoreMin = listStoreIn(1)
    numeStoreMax = listStoreIn(nbStoreIn)
    AS_DEALLOCATE(vi=listStoreIn)

! - Get model from command
    modelCmd = " "
    call getvid(' ', 'MODELE', scal=modelCmd, nbret=nocc)
    if (nocc .eq. 0) then
        modelCmd = " "
    end if

! - Get model from first storing index
    call rs_get_model(resultIn, numeStoreMin, modelRefe, codret)
    if (codret .lt. 0) then
        call utmess('F', 'CALCCHAMP1_44')
    end if

! - Get mesh
    call dismoi('NOM_MAILLA', modelRefe, 'MODELE', repk=mesh)
    l_pmesh = isParallelMesh(mesh)

! - Models are not consistents
    if (modelCmd .ne. " ") then
        if (modelRefe .ne. modelCmd) then
            call utmess('F', 'CALCCHAMP1_2')
        end if
    end if
    model = modelRefe

! - Prepare datastructure for management of resultat parameters
    call setResuPara(resultIn, resultOut, resultType, &
                     nbStore, listStore, listStoreJv, &
                     postCompResu)

! - Specific checks for beams
    call ccvepo(option, modelRefe, &
                postCompResu, &
                postCompPoux)

! - Create restricted list of cells for nodal option
    nbRestCell = 0
    if (option(6:9) .eq. 'NOEU') then
        call getelem(mesh, " ", 1, ' ', restCellJv, nbRestCell)
    end if
    lRestCell = nbRestCell .ne. 0

! - Management of multiple FED for partial computation of options
    call initCompRest(nbStore, postCompRest)

!     PREMIER PASSAGE POUR DETERMINER LES OPTIONS REELLEMENT A CALCULER
!     EN PRENANT EN COMPTE LA DEPENDANCE
!     PAR EXEMPLE SI SIGM_NOEU A BESOIN DE SIGM_ELNO QUI A BESOIN DE
!     SIGM_ELGA ET QUE SIGM_ELNO EST PRESENTE ALORS ON N'A PAS BESOIN
!     DE CALCULER SIGM_ELGA

! - Create list of storing index to compute all options
    do iOptEff = 1, nbOptEff
        optionEff = listOptEff(iOptEff)
        isOptionFromUser = (option .eq. optionEff)
        call cclord(iOptEff, nbStore, listStore, jvBaseName, isOptionFromUser, &
                    numeStoreMin, numeStoreMax, resultIn, resultOut)
    end do

! - Create object of flags for options to compute
    AS_ALLOCATE(vi=lacalc, size=nbOptEff)

! - First: everything has to been compute
    lacalc(1:nbOptEff) = 1

! - PUIS ON RETIRE LES OPTIONS DONT LE CALCUL N'EST PAS UTILE
    do iOptEff = nbOptEff-1, 1, -1
        optionEff = listOptEff(iOptEff)
        call cclodr(iOptEff, nbStore, listStore, jvBaseName, numeStoreMin, &
                    numeStoreMax, resultIn, resultOut, lacalc)
    end do

! - Loop on options to compute (dependencies)
    do iOptEff = 1, nbOptEff
! ----- Get current option
        optionEff = listOptEff(iOptEff)
        isOptionFromUser = (option .eq. optionEff)

! ----- This option does have to been computed ?
        if (lacalc(iOptEff) .eq. 0) cycle

! ----- Get information of storing indexes to compute
        call codent(iOptEff, 'D0', numeOptStr)
        listStoreOptJv = jvBaseName//'.OP'//numeOptStr
        call jeveuo(listStoreOptJv, 'L', jlinst)
        nbStoreToCompute = zi(jlinst)

! ----- Select JEVEUX base
        jvBase = 'G'
        if (isOptionFromUser) then
            if (present(jvBase_)) then
                jvBase = jvBase_
            end if
        else
            jvBase = 'V'
        end if

!       CE BLOC A ETE AJOUTE POUR LE CAS OU UNE OPTION1 A DECLENCHE
!       LE CALCUL D'UNE OPTION2 MAIS QUE CETTE OPTION2 EST ENSUITE
!       REDEMANDEE DANS LE MEME CALC_CHAMP PAR L'UTILISATEUR
        if (nbOpt .ne. 0) then
            posopt = indk16(listOpt, optionEff, 1, nbOpt)
            if (posopt .ne. 0) jvBase = 'G'
            if (.not. isOptionFromUser .and. posopt .ne. 0) then
                listOpt(posopt) = ' '
            end if
        end if
!
        if (isOptionFromUser .and. (nbStoreToCompute .eq. 0)) then
            call utmess('A', 'CALCCHAMP_1', sk=optionEff)
        end if

! SI PARALLELISME EN TEMPS: EVENTUELLE INITIALISATION CONTEXTE
        if (present(tldist_) .and. (optionEff(6:9) .eq. 'NOEU') &
            .and. (nbStoreToCompute .ge. 1)) then
            call pcptcc(101, ldist, dbg_ob, dbgv_ob, lcpu, ltest, rang, nbproc, mpicou, &
                        nbStoreToCompute, nbpas, vldist, vcham, k24b, ibid, k19b, &
                        k24b, k24b, lbid, &
                        ibid, ibid, ibid, ibid, ibid, &
                        k24b, ibid, ibid, kbid, k24b, prbid, pcbid)
            call jeveuo(vldist, 'L', jldist)
            if (ldist) call jeveuo(vcham, 'E', jvcham)
        else
! SI LA QUESTION NE SE POSE PAS (APPEL RECURSIF CCFNRN > CALCOP OU OPTION ELGA/ELNO)
            call pcptcc(102, ldist, dbg_ob, dbgv_ob, lcpu, ltest, rang, nbproc, mpicou, &
                        nbStoreToCompute, nbpas, vldist, vcham, k24b, ibid, k19b, &
                        k24b, k24b, lbid, &
                        ibid, ibid, ibid, ibid, ibid, &
                        k24b, ibid, ibid, kbid, k24b, prbid, pcbid)
            call jeveuo(vldist, 'L', jldist)
        end if

! ----- SI PARALLELISME EN TEMPS: ON DEBRANCHE L'EVENTUEL PARALLELISME EN ESPACE
        if (modelCmd .ne. " ") then
            if (modelRefe .ne. modelCmd) then
                call utmess('F', 'CALCCHAMP1_2')
            end if
        end if
        call pcptcc(2, ldist, dbg_ob, lbid, lbid, lbid, rang, ibid, mpibid, &
                    ibid, ibid, k24b, k24b, k24b, ibid, k19b, &
                    model, partsd, lsdpar, &
                    ibid, ibid, ibid, ibid, ibid, &
                    k24b, ibid, ibid, kbid, k24b, prbid, pcbid)
        ! if (nbproc .eq. 1 .and. niv > 1) then
        !     call utmess('I', 'PREPOST_25', sk=optionEff)
        ! else if (nbproc .gt. 1) then
        !     if (ldist) then
        !         ASSERT(.not. l_pmesh)
        !         call utmess('I', 'PREPOST_22', si=nbStore, sk=optionEff)
        !     elseif (.not. l_pmesh) then
        !         if (lsdpar) then
        !             call utmess('I', 'PREPOST_23', sk=optionEff)
        !         else
        !             call utmess('I', 'PREPOST_24', sk=optionEff)
        !         end if
        !     end if
        ! end if
!
        codre2 = 0
        ligrel = ' '
        ligrelSave = ' '
! SI PARALLELISME EN TEMPS: GESTION DE L'INDICE DE DECALAGE
        ipas = 1
        nbEqua = -9999
        do iStoreToCompute = 1, nbStoreToCompute
! --------- Current storing indexes
            numeStore = zi(jlinst+iStoreToCompute+2)
            numeStorePrev = numeStore-1

            if (((zi(jldist+iStoreToCompute-1) .eq. rang) .and. (ldist)) &
                .or. (.not. ldist)) then

! ------------- Prepare indexes for time distribution
                call pcptcc(4, ldist, dbg_ob, lbid, lbid, lbid, rang, nbproc, mpibid, &
                            ibid, nbpas, k24b, k24b, k24b, ibid, k19b, &
                            k24b, k24b, lbid, &
                            iStoreToCompute, ipas, ideb, ifin, irelat, &
                            k24b, ibid, ibid, kbid, k24b, prbid, pcbid)
                ligrelHasBeenChanged = .false.

                if (optionEff(6:9) .eq. 'NOEU') then
                    if ((ldist) .and. (ideb .ne. ifin)) then
! --------------------- SI PARALLELISME EN TEMPS ET NPAS NON ATTEINT: NBPROC CHAM_NOS SIMULTANES
                        p = 1
                        do k = ideb, ifin
                            numork = zi(jlinst+k+2)
! ------------------------- Get parameters
                            call medom1(modelNew, materField, materCode, caraElem, &
                                        listLoad, nbLoad, resultIn, numeStore)

! ------------------------- Get ligrel
                            ASSERT(option .ne. 'SIRO_ELEM')
                            call getLigrel(option, modelNew, resultOut, 'V', postCompRest, &
                                           ligrel)

! ------------------------- Some checks
                            call checkOption(optionEff, ligrel, resultType)

                            if (dbg_ob) then
                                write (ifm, *) '< ', rang, &
                                    'calcop> modele_avant/modele_apres=', model, modelNew
                            end if

                            if (model .ne. modelNew) then
                                call utmess('F', 'PREPOST_1')
                            end if
                            model = modelNew

! ------------------------- Get field from output result
                            call rsexc1(resultOut, optionEff, numork, nochok)
                            if (numork .eq. numeStore) then
                                nochou = nochok
                                ligrelHasBeenChanged = ligrelSave .ne. ligrel
                            end if
                            zk24(jvcham+p-1) = nochok

                            if (dbg_ob) then
                                write (ifm, *) '< ', rang, 'calcop> p/k/chamnk=', p, k, nochok
                            end if
                            p = p+1
                        end do

! --------------------- Compute nodal field
                        call ccchno(optionEff, numeStore, resultIn, resultOut, fieldNameOut, &
                                    lRestCell, nbRestCell, restCellJv, &
                                    mesh, model, caraElem, jvBase, &
                                    ligrel, ligrelHasBeenChanged, codre2, &
                                    nochou, &
                                    ideb_=ideb, ifin_=ifin, vcham_=vcham)
                    else
! --------------------- SINON, 1 SEUL A LA FOIS
! --------------------- SI PARALLELISME EN TEMPS et NPAS ATTEINT (RELIQUAT DE PAS DE TEMPS)
! --------------------- OU SI NON PARALLELISME EN TEMPS
! --------------------- Get parameters
                        call medom1(modelNew, materField, materCode, caraElem, &
                                    listLoad, nbLoad, resultIn, numeStore)

! --------------------- Get ligrel
                        ASSERT(option .ne. 'SIRO_ELEM')
                        call getLigrel(option, modelNew, resultOut, 'V', postCompRest, &
                                       ligrel)

! --------------------- Some checks
                        call checkOption(optionEff, ligrel, resultType)

! --------------------- SI PARALLELISME EN TEMPS: CAS DE FIGURE DU RELIQUAT DE PAS DE TEMPS
! --------------------- (AU CAS OU)
                        if (ldist) then
                            if (dbg_ob) then
                                write (ifm, *) '< ', rang, &
                                    'calcop> modele_avant/modele_apres=', model, modelNew
                            end if
                            if (model .ne. modelNew) then
                                call utmess('F', 'PREPOST_1')
                            end if
                        end if
                        model = modelNew
                        ligrelHasBeenChanged = ligrelSave .ne. ligrel

                        call rsexc1(resultOut, optionEff, numeStore, nochou)
                        call ccchno(optionEff, numeStore, resultIn, resultOut, fieldNameOut, &
                                    lRestCell, nbRestCell, restCellJv, &
                                    mesh, model, caraElem, jvBase, &
                                    ligrel, ligrelHasBeenChanged, codre2, &
                                    nochou)
                    end if
! ----------------- SI PARALLELISME EN TEMPS:  COM MPI CHAM_NOS.VALE DONT LES NOMS
! ----------------- SONT STOCKES DANS VCHAM
                    chamno = nochou(1:19)//'.VALE'
                    if (ldist) then
                        call dismoi('TYPE_SCA', chamno(1:19), 'CHAM_NO', repk=ktyp)
                        call jelira(chamno, 'LONMAX', nbEquaNew)
! --------------------- SI PARALLELISME EN TEMPS:
! --------------------- POUR L'INSTANT, ON SUPPOSE QUE TOUS LES CHAM_NOS SONT DE LONGUEUR IDENTIQUE
! --------------------- ON TESTE SI C'EST LE CAS SUR LES NBPROCS PAS DE TEMPS CONTIGUES ET
! --------------------- SUR LE PAS PRECEDENT
                        call pcptcc(6, ldist, dbg_ob, lbid, lbid, lbid, rang, ibid, mpibid, &
                                    ibid, ibid, k24b, k24b, k24b, ibid, k19b, &
                                    k24b, k24b, lbid, &
                                    ibid, ipas, ibid, ibid, ibid, &
                                    k24b, nbEquaNew, nbEqua, kbid, k24b, prbid, pcbid)
                        nbEqua = nbEquaNew

! --------------------- Access to nodal field
                        if (ktyp .eq. 'R') then
                            call jeveuo(chamno, 'L', vr=noch)
                        else if (ktyp .eq. 'C') then
                            call jeveuo(chamno, 'L', vc=nochc)
                        else
                            ASSERT(ASTER_FALSE)
                        end if
                        if (dbg_ob) then
                            write (ifm, *) '< ', rang, &
                                'calcop> chamno/ktyp/nbEqua=', chamno, ktyp, nbEqua
                        end if
                    end if

                    call pcptcc(7, ldist, dbg_ob, lbid, lbid, lbid, rang, nbproc, mpicou, &
                                ibid, ibid, k24b, vcham, k24b, ibid, k19b, &
                                k24b, k24b, lbid, &
                                ibid, ipas, ideb, ifin, irelat, &
                                k24b, ibid, nbEqua, ktyp, vcnoch, noch, nochc)

! ----------------- PARALLELISME EN TEMPS: TEST DE VERIFICATION
                    call pcptcc(8, ldist, lbid, dbgv_ob, lbid, lbid, ibid, ibid, mpibid, &
                                ibid, ibid, k24b, vcham, k24b, ibid, k19b, &
                                k24b, k24b, lbid, &
                                ibid, ibid, ideb, ifin, ibid, &
                                chamno, ibid, ibid, kbid, k24b, prbid, pcbid)
!
                else if (optionEff(6:7) .eq. 'EL') then

! ----------------- Get parameters
                    call medom1(model, materField, materCode, caraElem, &
                                listLoad, nbLoad, resultIn, numeStore)

! ----------------- Get ligrel
                    call getLigrel(option, model, resultOut, jvBase, postCompRest, &
                                   ligrel)

! ----------------- Some checks
                    call checkOption(optionEff, ligrel, resultType)

! ----------------- Compute elementary field
                    call ccchel(optionEff, &
                                model, materField, materCode, caraElem, listLoad, &
                                resultIn, resultOut, resultType, &
                                numeStore, numeStorePrev, &
                                ligrel, isTransient, postCompPoux, jvBase, &
                                fieldNameOut)
                    if (fieldNameOut .eq. ' ') goto 20

                else
                    ASSERT(ASTER_FALSE)
                end if

! ------------- Save field in result datastructure
                call exisd('CHAMP_GD', fieldNameOut, iret)
                if (jvBase .eq. 'G') then
                    if (iret .eq. 0) then
                        codret = 1
                        call utmess('A', 'CALCULEL2_89', sk=optionEff)
                    else
                        if ((ldist) .and. (ideb .ne. ifin)) then
                            do k = ideb, ifin
                                numork = zi(jlinst+k+2)
                                call rsnoch(resultOut, optionEff, numork)
                            end do
                        else
                            call rsnoch(resultOut, optionEff, numeStore)
                        end if
                    end if
                end if

! ------------- Cleaning
                if (postCompPoux%lPoux) then
                    call jedetc('V', '&&MECHPO', 1)
                end if
                call detrsd('CHAM_ELEM_S', fieldNameOut)
                ligrelSave = ligrel
                if (ldist) then
                    ipas = ipas+1
                end if
            end if
        end do
! SI PARALLELISME EN TEMPS: NETTOYAGE DU CONTEXTE
        call pcptcc(301, ldist, dbg_ob, lbid, lbid, lbid, rang, ibid, mpibid, &
                    ibid, ibid, vldist, vcham, k24b, ibid, k19b, &
                    model, partsd, lsdpar, &
                    ibid, ibid, ibid, ibid, ibid, &
                    k24b, ibid, ibid, kbid, vcnoch, prbid, pcbid)
20      continue
    end do
!
    codret = 0

! - Clean objects
    AS_DEALLOCATE(vi=lacalc)
    call ccnett(jvBaseName, nbOptEff)
    if (option(6:9) .eq. 'NOEU' .and. nbRestCell .ne. 0) call jedetr(restCellJv)
!
999 continue
!
    call deleteCompRest(postCompRest)
!
    call jedema()
!
end subroutine
