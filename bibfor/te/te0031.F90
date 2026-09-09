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
subroutine te0031(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    use plateMaterial_module, only: chckMultiLayer
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/cosiro.h"
#include "asterfort/dkqmas.h"
#include "asterfort/dkqrig.h"
#include "asterfort/dktmas.h"
#include "asterfort/dktnli.h"
#include "asterfort/dktrig.h"
#include "asterfort/dsqmas.h"
#include "asterfort/dsqrig.h"
#include "asterfort/dstmas.h"
#include "asterfort/dstrig.h"
#include "asterfort/dxbsig.h"
#include "asterfort/dxeffi.h"
#include "asterfort/dxiner.h"
#include "asterfort/dxroep.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nmtstm.h"
#include "asterfort/pmavec.h"
#include "asterfort/q4gmas.h"
#include "asterfort/q4grig.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/t3grig.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpslg2.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT, DST, Q4G
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          REFE_FORC_NODA, FORC_NODA
!          EPOT_ELEM, ECIN_ELEM
!          MASS_MECA, MASS_MECA_DIAG, MASS_MECA_EXPLI, M_GAMMA, MASS_INER
!
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: npge = 3, nbEfgeNd = 8
    integer(kind=8) :: codret, nno
    integer(kind=8) :: jvDispM, jvDispIncr
    integer(kind=8) :: i, j, k
    integer(kind=8) :: jvGeom, jvCompor, jvOmega, jvAcce
    integer(kind=8) :: jvCodret, jvVect, jvMatr, jvEner, jvMassIner
    integer(kind=8) :: nddl, nbTermSyme, iret, jvSief
    integer(kind=8) :: nbLayer, itab(7), nbsp
    integer(kind=8) :: n1, n2, ni, elasID
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), bsigmEner(24), effgt(32)
    real(kind=8) :: effref, momref
    real(kind=8) :: vecloc(24), ener(3), matp(24, 24), matv(300)
    real(kind=8) :: foref, moref
    character(len=16) :: defoComp
    aster_logical :: lElasAniso, lNonLine
!     ---> POUR DKT/DST MATELEM = 3 * 6 DDL = 171 TERMES STOCKAGE SYME
!     ---> POUR DKQ/DSQ MATELEM = 4 * 6 DDL = 300 TERMES STOCKAGE SYME
    real(kind=8) :: matrRigi(576), matrMass(576), matrTang(576)
    real(kind=8) :: rho, epais
!     --->   UML : DEPLACEMENT A L'INSTANT T- (REPERE LOCAL)
!     --->   DUL : INCREMENT DE DEPLACEMENT   (REPERE LOCAL)
    real(kind=8) :: uml(6, 4), dul(6, 4)
    aster_logical :: lVect, lMatr, lVari, lSigm, matsym, lComposite
    character(len=8), parameter :: typmod(2) = (/'C_PLAN  ', '        '/)
    integer(kind=8) :: jvInstmr, jvInstpr
    real(kind=8) :: instm, instp
    integer(kind=8) :: jvCarcri, jvMaterc
    character(len=16), pointer :: compor(:) => null()
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie, plateOrieSave
    type(RESI_REFE) :: refe
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, nno=nno)
    ASSERT(nno .eq. 3 .or. nno .eq. 4)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - Change frame for input stress
    if (option .eq. 'FORC_NODA') then
        call cosiro(plateCara, plateOrie, &
                    'PSIEFR', 'L', 'UI', 'G')
    elseif (option .ne. 'REFE_FORC_NODA') then
        call cosiro(plateCara, plateOrie, &
                    'PCONTMR', 'L', 'UI', 'G')
        call cosiro(plateCara, plateOrie, &
                    'PCONTRR', 'L', 'UI', 'G')
    end if
    plateOrieSave = plateOrie
!
    lNonLine = (option(1:9) .eq. 'FULL_MECA') .or. (option .eq. 'RAPH_MECA') .or. &
               (option(1:10) .eq. 'RIGI_MECA_')

! - Material parameters
    if (lNonLine .or. option(1:9) .eq. 'RIGI_MECA') then
! ----- Get material parameters
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), materPara)

! ----- No local coordinate system from user
        call initLCSNone(materPara)
    end if

! - Anisotropic case ?
    lElasAniso = ASTER_FALSE
    if (lNonLine .or. option(1:9) .eq. 'RIGI_MECA') then
        call jevech('PMATERC', 'L', jvMaterc)
        call get_elas_id(zi(jvMaterc), elasID)
        lElasAniso = elasID .eq. ELAS_SHELL .or. &
                     elasID .eq. ELAS_COMPOSITE .or. &
                     elasID .eq. ELAS_ORTH
    end if
    if (lNonLine) then
        if (lElasAniso) then
            call utmess('F', 'PLATE1_75')
        end if
    end if

! - Composite ?
    lComposite = ASTER_FALSE
    if (lNonLine .or. option(1:9) .eq. 'RIGI_MECA') then
        lComposite = materPara%elasID .eq. ELAS_COMPOSITE
    end if

! - Check consistency between DEFI_COQU_MULT/AFFE_CARA_ELEM
    if (lComposite) then
        ASSERT(lElasAniso)
        call chckMultiLayer(materPara, plateCara)
    end if

    if (option .eq. 'RIGI_MECA' .or. &
        option .eq. 'EPOT_ELEM') then
        if (nomte .eq. 'MEDKTR3') then
            call dktrig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrRigi, ener)
        else if (nomte .eq. 'MEDSTR3') then
            call dstrig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrRigi, ener)
        else if (nomte .eq. 'MEDKQU4') then
            call dkqrig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrRigi, ener)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqrig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrRigi, ener)
        else if (nomte .eq. 'MEQ4QU4') then
            call q4grig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrRigi, ener)
        else if (nomte .eq. 'MET3TR3') then
            call t3grig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrRigi, ener)
        end if
        if (option .eq. 'RIGI_MECA') then
            call jevech('PMATUUR', 'E', jvMatr)
            call utpslg(nno, 6, pgl, matrRigi, zr(jvMatr))
        else if (option .eq. 'EPOT_ELEM') then
            call jevech('PENERDR', 'E', jvEner)
            do i = 1, 3
                zr(jvEner-1+i) = ener(i)
            end do
        end if

    else if ((option .eq. 'MASS_MECA') .or. (option .eq. 'MASS_MECA_DIAG') .or. &
             (option .eq. 'MASS_MECA_EXPLI') .or. (option .eq. 'M_GAMMA') .or. &
             (option .eq. 'ECIN_ELEM')) then
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MET3TR3') then
            call dktmas(plateCara, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        else if (nomte .eq. 'MEDSTR3') then
            call dstmas(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        else if (nomte .eq. 'MEDKQU4') then
            call dkqmas(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqmas(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gmas(plateCara, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        end if
        if (option .eq. 'MASS_MECA') then
            call jevech('PMATUUR', 'E', jvMatr)
            call utpslg(nno, 6, pgl, matrMass, zr(jvMatr))
        else if (option .eq. 'ECIN_ELEM') then
            call jevech('PENERCR', 'E', jvEner)
            call jevech('POMEGA2', 'L', jvOmega)
            do i = 1, 3
                zr(jvEner-1+i) = zr(jvOmega)*ener(i)
            end do
        else if (option .eq. 'M_GAMMA') then
            call jevech('PACCELR', 'L', jvAcce)
            call jevech('PVECTUR', 'E', jvVect)
            nddl = 6*nno
            nbTermSyme = nddl*(nddl+1)/2
            call utpslg(nno, 6, pgl, matrMass, matv)
            call vecma(matv, nbTermSyme, matp, nddl)
            call pmavec('ZERO', nddl, matp, zr(jvAcce), zr(jvVect))
        else if (option .eq. 'MASS_MECA_DIAG' .or. option .eq. 'MASS_MECA_EXPLI') then
            call jevech('PMATUUR', 'E', jvMatr)
            nddl = 6*nno
            nbTermSyme = nddl*(nddl+1)/2
            do i = 1, nbTermSyme
                zr(jvMatr-1+i) = matrMass(i)
            end do
            if (option .eq. 'MASS_MECA_EXPLI') then
!               CORRECTION DES TERMES CORRESPONDANT AU DDL 6
!               NON PREVU PAR LA THEORIE DKT. ON RAJOUTE
!               UN TERME DIAGONAL NON ZERO EGAL A CELUI DU DDL 5.
!               CETTE CORRECTION A ETE INSPIRE PAR LA DEMARCHE DANS EUROPLEXUS
                do j = 1, nno
                    n1 = 6*(j-1)+5
                    n2 = 6*(j-1)+4
                    ni = 6*j
                    nbTermSyme = (ni+1)*ni/2
                    n1 = (n1+1)*n1/2
                    n2 = (n2+1)*n2/2
                    zr(jvMatr-1+nbTermSyme) = (zr(jvMatr-1+n1)+zr(jvMatr-1+n2))*0.5d0
                end do
            end if
        end if

    else if (option .eq. 'MASS_INER') then
        call jevech('PMASSINE', 'E', jvMassIner)
        call dxroep(plateCara, rho, epais)
        call dxiner(platecara, &
                    zr(jvGeom), rho, epais, zr(jvMassIner), &
                    zr(jvMassIner+1), zr(jvMassIner+4))

    else if (lNonLine) then
        call jevech('PDEPLMR', 'L', jvDispM)
        call jevech('PDEPLPR', 'L', jvDispIncr)
        call jevech('PINSTMR', 'L', jvInstmr)
        call jevech('PINSTPR', 'L', jvInstpr)
        instm = zr(jvInstmr)
        instp = zr(jvInstpr)

! ----- Get fields for non-linear behaviour
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', jvCarcri)
        defoComp = compor(DEFO)

! ----- Initialisation of behaviour datastructure
        call behaviourInit(BEHInteg)

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  instm, instp, &
                                  materPara, BEHInteg)

! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

! ----- Update configuration
        if (defoComp .eq. 'GROT_GDEP') then
            do i = 1, nno
                zr(jvGeom+3*(i-1)) = zr(jvGeom+3*(i-1))+ &
                                     zr(jvDispM+6*(i-1))+zr(jvDispIncr+6*(i-1))
                zr(jvGeom+3*(i-1)+1) = zr(jvGeom+3*(i-1)+1)+ &
                                       zr(jvDispM+6*(i-1)+1)+zr(jvDispIncr+6*(i-1)+1)
                zr(jvGeom+3*(i-1)+2) = zr(jvGeom+3*(i-1)+2)+ &
                                       zr(jvDispM+6*(i-1)+2)+zr(jvDispIncr+6*(i-1)+2)
            end do

! --------- Calculate the transformation: global coordinate system/intrinsic coordinate system
            call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! --------- Compute coordinate system for plate
            call compCoorSystPlate(pgl, plateCara, plateOrie)

! --------- Change coordinates of geometry
            call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

        end if

! ----- Change frame for displacements
        call utpvgl(nno, 6, pgl, zr(jvDispM), uml)
        call utpvgl(nno, 6, pgl, zr(jvDispIncr), dul)

! ----- Compute non-linear options
        call dktnli(plateCara, plateOrie, &
                    BEHInteg, option, typmod, &
                    instm, instp, &
                    xyzl, uml, dul, &
                    vecloc, matrTang, codret)

! ----- Output fields
        if (lMatr) then
            call nmtstm(zr(jvCarcri), jvMatr, matsym)
            if (matsym) then
                call utpslg(nno, 6, pgl, matrTang, zr(jvMatr))
            else
                call utpslg2(nno, 6, pgl, matrTang, zr(jvMatr))
            end if
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', jvVect)
            call utpvlg(nno, 6, pgl, vecloc, zr(jvVect))
        end if
        if (lSigm) then
            call jevech('PCODRET', 'E', jvCodret)
            zi(jvCodret) = codret
        end if

    else if (option .eq. 'FORC_NODA') then
        effgt = 0.d0
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, itab=itab)
        jvSief = itab(1)
        nbsp = itab(7)
        nbLayer = plateCara%nbLayer
        if (nbsp .ne. npge*nbLayer) then
            call utmess('F', 'PLATE1_4')
        end if
!
        call dxeffi(plateCara, plateOrie, &
                    option, nomte, zr(jvSief), nbEfgeNd, &
                    effgt)
!
        call tecach('NNO', 'PCOMPOR', 'L', iret, iad=jvCompor)
        if (jvCompor .ne. 0) then
            defoComp = zk16(jvCompor-1+DEFO)
            if (defoComp .eq. 'GROT_GDEP') then
                call jevech('PDEPLAR', 'L', jvDispM)

! ------------- Update configuration
                do i = 1, nno
                    zr(jvGeom+3*(i-1)) = zr(jvGeom+3*(i-1))+zr(jvDispM+6*(i-1))
                    zr(jvGeom+3*(i-1)+1) = zr(jvGeom+3*(i-1)+1)+zr(jvDispM+6*(i-1)+1)
                    zr(jvGeom+3*(i-1)+2) = zr(jvGeom+3*(i-1)+2)+zr(jvDispM+6*(i-1)+2)
                end do

! ------------- Calculate the transformation
                call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! ------------- Compute coordinate system for plate
                call compCoorSystPlate(pgl, plateCara, plateOrie)

! ------------- Change coordinates of geometry
                call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

            end if
        end if

! ----- CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
        call dxbsig(plateCara, plateOrie, &
                    nomte, option, &
                    xyzl, pgl, effgt, &
                    bsigmEner)

! ----- AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
        call jevech('PVECTUR', 'E', jvVect)
        k = 0
        do i = 1, nno
            do j = 1, 6
                k = k+1
                zr(jvVect+k-1) = bsigmEner(k)
            end do
        end do
!
    else if (option .eq. 'REFE_FORC_NODA') then
        call refe%Init(nomte)
        foref = refe%GetRef('EFFORT')
        moref = refe%GetRef('MOMENT')
        call refe%Check()
        do i = 1, nno
            do j = 1, 3
                effgt(nbEfgeNd*(i-1)+j) = foref
                effgt(nbEfgeNd*(i-1)+3+j) = moref
            end do
            effgt(nbEfgeNd*(i-1)+7) = 0.0d0
            effgt(nbEfgeNd*(i-1)+8) = 0.0d0
        end do

! ----- CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
        call dxbsig(plateCara, plateOrie, &
                    nomte, option, &
                    xyzl, pgl, effgt, &
                    bsigmEner)

! ----- AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
        call jevech('PVECTUR', 'E', jvVect)
        k = 0
        do i = 1, nno
            effref = (abs(bsigmEner(k+1))+abs(bsigmEner(k+2))+abs(bsigmEner(k+3)))/3.d0
            momref = (abs(bsigmEner(k+4))+abs(bsigmEner(k+5))+abs(bsigmEner(k+6)))/3.d0
            do j = 1, 6
                k = k+1
                if (j .lt. 4) then
                    zr(jvVect+k-1) = effref
                else
                    zr(jvVect+k-1) = momref
                end if
            end do
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (option .ne. 'REFE_FORC_NODA') then
        call cosiro(plateCara, plateOrieSave, &
                    'PCONTPR', 'E', 'IU', 'G')
    end if
!
end subroutine
