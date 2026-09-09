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
! aslint: disable=W1501
!
subroutine te0491(option, nomte)
!
    use Behaviour_module
    use BehaviourStrain_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/diago3.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/enelpg.h"
#include "asterfort/eps1mc.h"
#include "asterfort/epsvmc.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nmgeom.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte, option
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
!
! Options: ENEL_ELEM, ENTR_ELEM, ENER_TOTALE, INDIC_ENER, INDIC_SEUIL
!
! --------------------------------------------------------------------------------------------------
!
!  OPTION INDIC_ENER : CALCUL DE  L'INDICATEUR GLOBAL
!  =================   ENERGETIQUE DETERMINE PAR L'EXPRESSION SUIVANTE:
!
!            IE = (SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV)/V
!
!        OU  .OMEGA EST LA DENSITE D'ENERGIE TOTALE
!            (I.E. OMEGA = SOMME_0->T(SIGMA:D(EPS)/DT).DTAU
!            .PSI EST LA DENSITE D'ENERGIE ELASTIQUE 'TOTALE'
!            (I.E. ASSOCIEE A LA COURBE DE TRACTION SI ON
!                  CONSIDERAIT LE MATERIAU ELASTIQUE NON-LINEAIRE)
!            .V EST LE VOLUME DU GROUPE DE MAILLES TRAITE
!
!  OPTION INDIC_SEUIL : CALCUL DE  L'INDICATEUR GLOBAL
!  ==================   DETERMINE PAR L'EXPRESSION SUIVANTE :
!
!   IS = (SOMME_DOMAINE(1 - ((SIG-X):EPS_PLAST)/((SIG_Y+R)*P)).DV)/V
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!            .X         EST LE TENSEUR DE RAPPEL
!            .EPS_PLAST EST LE TENSEUR DES DEFORMATIONS PLASTIQUES
!            .SIG_Y     EST LA LIMITE D'ELASTICITE
!            .R         EST LA FONCTION D'ECROUISSAGE
!            .P         EST LA DEFORMATION PLASTIQUE CUMULEE
!            .V         EST LE VOLUME DU GROUPE DE MAILLES TRAITE
!
!  OPTION ENEL_ELEM : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  ================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!  EN HPP
!   ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES DE CAUCHY
!            .D         EST LE TENSEUR DE HOOKE
!
!  EN GRANDES DEFORMATIONS SIMO MIEHE POUR ELAS OU VMIS_ISOT
!   ENERLAS = ENERGIE ELASTIQUE SPECIFIQUE
!           = K/2(0.5(J^2-1)-lnJ)+0.5mu(tr(J^(-2/3)be)-3)
!           SI PRESENCE DE THERMIQUE, ON AJOUTE UNE CORRECTION
!           SPECIFIQUE PRESENTEE DANS LA DOC R
!  EN GRANDES DEFORMATIONS GDEF_LOG
!   ENERELAS = SOMME_VOLUME((T_T*(1/D)*T).DV)
!        OU  .T       EST LE TENSEUR DES CONTRAINTES DU FORMALISME
!            .D         EST LE TENSEUR DE HOOKE
!
!  OPTION ENTR_ELEM : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  =================   MODIFIEE DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!   ENELAS =  0.5*Lame*H(tr(EPS))*tr(EPS)**2+mu*SUM(H(Ei)*Ei**2)
!
!        OU  .EPS      EST LE TENSEUR DES DEFORMATIONS ELASTIQUES
!            .Ei       SONT (pour i=1..3) LES DEFORMATIONS PROPRES
!            .H        LA FONCTION D'HEAVISIDE
!
!
!  OPTION ENER_TOTALE : CALCUL DE L'ENERGIE DE DEFORMATION TOTALE
!   ENER_TOTALE =  ENELAS + EPLAS
!
!          AVEC : ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!                 ENELAS EST L'ENERGIE DE DEFORMATION ELASTIQUE
!
!           OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!               .D         EST LE TENSEUR DE HOOKE
!
!          ET   : EPLAS = SOMME_VOLUME((R(P))*D(P))
!                 EPLAS EST L'ENERGIE DE DEFORMATION PLASTIQUE
!
!           OU  .P         EST LA DEFORMATION PLASTIQUE CUMULEE
!           ET   R(P) EST CALCULE POUR LES COMPORTEMENTS SUIVANTS :
!                      .VMIS_ISOT_LINE
!                      .VMIS_ISOT_TRAC
!
!   REMARQUE : EN GRANDE DEFORMATION ON INTEGRE SUR LE VOLUME INITIAL
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: axi = ASTER_FALSE
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8), parameter :: mxcmel = 162, nbsgm = 6
    real(kind=8), parameter :: zero = 0.d0, undemi = 0.5d0, un = 1.d0
    real(kind=8), parameter :: deux = 2.d0, trois = 3.d0, untier = 1.d0/3.d0
    integer(kind=8), parameter :: nbProp = 2
    integer(kind=8) :: propCodret(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'D_SIGM_EPSI', 'SY         '/)
    real(kind=8) :: propVale(nbProp)
    real(kind=8), parameter :: nharm = 0.d0
    integer(kind=8) :: idconm, idene1, idene2, jvDisp, ideplm, jvDispPrev
    integer(kind=8) :: jvDBaseFunc, jvSigm, jvSigmPrev, jvVari, kpg, jvGeom, jvMaterc
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jprol, jvale, jvTime
    integer(kind=8) :: nbsig, nbsig2, nbval, nbvari, ndim, nno, npg, nbEpsi
    integer(kind=8) :: iret, iret1, i, jtab(7)
    real(kind=8) :: airep, c1, c2, deuxmu, dsde, e
    real(kind=8) :: lame, mu, vecp(3, 3), epm(3)
    real(kind=8) :: enerElas, eneldv, enelsp, enelto, eplaeq, enerPlas, epseq
    real(kind=8) :: enerElasTrac, enerElas1, enerElas2, trepstraction, welastr
    real(kind=8) :: omega, p, poids, psi, rp
    real(kind=8) :: rprim, sigeq, sigy, tempg, trepsm, trsig
    real(kind=8) :: volume, welas, wtotal
    real(kind=8) :: sigmEner(nbsgm), epsdv(nbsgm)
    real(kind=8) :: epsiElas(nbsgm), epsiPlas(nbsgm), x(nbsgm)
    real(kind=8) :: epsim(nbsgm), delta(nbsgm), sigmm(nbsgm)
    real(kind=8) :: epsi(nbsgm), epssm(mxcmel), epss(mxcmel)
    real(kind=8) :: time, integ, integ1
    real(kind=8) :: epsiMeca(mxcmel), integ2, nu, k, indigl, para_vale
    real(kind=8) :: f(3, 3), r, trav(81)
    character(len=8) :: para_type
    character(len=16) :: relaName, defoComp
    aster_logical :: largeStrain
    integer(kind=8) :: strainType
    aster_logical :: lStrainMeca
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    omega = zero
    psi = zero
    volume = zero
    indigl = zero
    enerElas = zero
    enerElasTrac = zero
    enerPlas = zero
    welas = zero
    welastr = zero
    wtotal = zero

! - CARACTERISTIQUES DU TYPE D'ELEMENT
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)

! - NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT :
    nbsig = nbsigm()
    nbEpsi = nbsig
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Get displacements
    call jevech('PDEPLR', 'L', jvDisp)

! - Get stresses
    call jevech('PCONTPR', 'L', jvSigm)

! - Get time
    time = r8vide()
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    end if

! - Get parameters from behaviour (linear and non-linear cases)
    call behaviourGetParameters(relaName, defoComp)
    if ((defoComp .eq. 'SIMO_MIEHE') .or. (defoComp .eq. 'GDEF_LOG')) then
        largeStrain = ASTER_TRUE
    else
        largeStrain = ASTER_FALSE
    end if
    if (option .eq. 'ENER_TOTALE') then
        if (largeStrain) then
            call utmess('F', 'ENERGY1_3', sk=defoComp)
        end if
    end if

! - ON TESTE LA RECUPERATION DU CHAMP DE CONTRAINTES DU PAS PRECEDENT
    if (option .eq. 'ENER_TOTALE') then
        if ((relaName(1:9) .ne. 'VMIS_ISOT') .and. (relaName(1:4) .ne. 'ELAS')) then
            call tecach('NNO', 'PCONTMR', 'L', iret, iad=idconm)
            if (idconm .ne. 0) then
                call jevech('PCONTMR', 'L', jvSigmPrev)
            end if
        end if
    end if

! - ON TESTE LA RECUPERATION DU CHAMP DE CONTRAINTES DU PAS PRECEDENT
    if (option .eq. 'ENER_TOTALE') then
        call tecach('NNO', 'PDEPLM', 'L', iret, iad=ideplm)
        if (ideplm .ne. 0) then
            call jevech('PDEPLM', 'L', jvDispPrev)
        end if
    end if

! - Internal state variables (only for non-linear)
    jvVari = 1
    nbvari = 0
    call tecach('ONO', 'PVARIPR', 'L', iret, nval=7, itab=jtab)
    if (iret .eq. 0) then
        jvVari = jtab(1)
        nbvari = max(jtab(6), 1)*jtab(7)
    end if

! - CALCUL DES DEFORMATIONS TOTALES DANS LE CAS DE ENERGIE TOTALE A L INSTANT COURANT ET
! - CELUI D AVANT
    if ((relaName(1:9) .ne. 'VMIS_ISOT') .and. (relaName(1:4) .ne. 'ELAS')) then
        if (option .eq. 'ENER_TOTALE') then
            call eps1mc(nno, ndim, nbsig, npg, &
                        jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                        zr(jvGeom), zr(jvDisp), nharm, &
                        epss)
            if (ideplm .ne. 0) then
                call eps1mc(nno, ndim, nbsig, npg, &
                            jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                            zr(jvGeom), zr(jvDispPrev), nharm, &
                            epssm)
            end if
        end if
    end if

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Compute mechanical strains
    strainType = STRAIN_TYPE_SMALL
    lStrainMeca = ASTER_TRUE
    call epsvmc(nno, ndim, nbEpsi, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), zr(jvDisp), &
                time, nharm, &
                strainType, lStrainMeca, materPara, &
                epsiMeca)

    if (option .eq. 'INDIC_ENER' .or. option .eq. 'ENEL_ELEM' .or. &
        option .eq. 'ENER_TOTALE' .or. option .eq. 'ENTR_ELEM') then
        do kpg = 1, npg
            omega = zero
            psi = zero

! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, materPara)

!---------- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT
            sigmEner = 0.d0
            do i = 1, nbsig
                sigmEner(i) = zr(jvSigm+(kpg-1)*nbsig+i-1)
            end do

! --------- CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT :
            call nmgeom(3, nno, axi, largeStrain, zr(jvGeom), kpg, &
                        jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                        zr(jvDisp), &
                        ASTER_TRUE, poids, trav, f, epsi, &
                        r)

! --------- CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
            call enelpg(materPara, time, &
                        relaName, defoComp, &
                        f, sigmEner, &
                        nbvari, zr(jvVari+(kpg-1)*nbvari), &
                        enerElas)

            if (option .eq. 'ENEL_ELEM') then
                welas = welas+enerElas*poids
                goto 10
            end if

! --------- Get elastic parameters
            if (materPara%elasID .ne. ELAS_ISOT) then
                call utmess("F", "ENERGY1_4", nk=2, valk=[option, materPara%elasKeyword])
            end if
            call get_elas_para(fami, zi(jvMaterc), '+', kpg, ksp, &
                               materPara%elasID, materPara%elasKeyword, &
                               e_=e, nu_=nu)
            deuxmu = e/(un+nu)
            k = untier*e/(un-deux*nu)
            lame = (e*nu)/((1+nu)*(1-2*nu))
            mu = e/(2*(1+nu))

! --- CALCUL DES DEFORMATIONS ELASTIQUES EN CONSIDERANT LE MATERIAU ISOTROPE :
! --- EPS_ELAS    = 1/D*SIGMA
! ---             = ((1+NU)/E)*SIGMA-(NU/E)*TRACE(SIGMA) :
            c1 = (un+nu)/e
            c2 = nu/e

            trsig = sigmEner(1)+sigmEner(2)+sigmEner(3)
            epsiElas(1) = c1*sigmEner(1)-c2*trsig
            epsiElas(2) = c1*sigmEner(2)-c2*trsig
            epsiElas(3) = c1*sigmEner(3)-c2*trsig
            epsiElas(4) = c1*sigmEner(4)
            epsiElas(5) = c1*sigmEner(5)
            epsiElas(6) = c1*sigmEner(6)

! --------- PARTIE SPHERIQUE DE L'ENERGIE DE DEFORMATION ELASTIQUE POSITIVE OU NULLE
            trepstraction = epsiElas(1)+epsiElas(2)+epsiElas(3)
            if (trepstraction .le. zero) then
                trepstraction = zero
            end if
            enerElas1 = undemi*lame*trepstraction*trepstraction

! --------- CALCUL DES DEFORMATIONS ELASTIQUES PRINCIPALES POSITIVES OU NULLES:
            call diago3(epsiElas, vecp, epm)
            if (epm(1) .le. zero) then
                epm(1) = zero
            end if
            if (epm(2) .le. zero) then
                epm(2) = zero
            end if
            if (epm(3) .le. zero) then
                epm(3) = zero
            end if

! --------- PARTIE DE L'ENERGIE DE DEFORMATION ELASTIQUE (DIRECTIONS PRINCIPALES) POSITIVE OU NULLE
            enerElas2 = mu*(epm(1)*epm(1)+epm(2)*epm(2)+epm(3)*epm(3))

! --------- ENERGIE DE DEFORMATION ELASTIQUE DE TRACTION
            enerElasTrac = enerElas1+enerElas2

            if (option .eq. 'ENTR_ELEM') then
                welastr = welastr+enerElasTrac*poids
                goto 10
            end if

! --------- Get elastic parameters
            if (materPara%elasID .ne. ELAS_ISOT) then
                call utmess("F", "ENERGY1_4", nk=2, valk=[option, materPara%elasKeyword])
            end if
            call get_elas_para(fami, zi(jvMaterc), '+', kpg, ksp, &
                               materPara%elasID, materPara%elasKeyword, &
                               e_=e, nu_=nu)

! --------- CALCUL DU TERME OMEGA REPRESENTANT L'ENERGIE TOTALE
! --------- OMEGA = SOMME_0->T(SIGMA:D(EPS)/DT).DTAU
            if (relaName .eq. 'VMIS_ISOT_LINE') then
! ------------- Parameters from linear traction
                call rcvalb(fami, kpg, ksp, '+', zi(jvMaterc), &
                            ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCodret, 2)
                dsde = propVale(1)
                sigy = propVale(2)

! ------------- RECUPERATION DE LA DEFORMATION PLASTIQUE CUMULEE
                p = zr(jvVari+(kpg-1)*nbvari+1-1)

! ------------- PENTE DE LA COURBE DE TRACTION DANS LE DIAGRAMME 'REDRESSE'
                rprim = e*dsde/(e-dsde)

! ------------- CONTRAINTE UNIAXIALE SUR LA COURBE DE TRACTION
                rp = sigy+rprim*p

! ------------- TRAVAIL PLASTIQUE 'EQUIVALENT'
                enerPlas = undemi*(sigy+rp)*p

            else if (relaName .eq. 'VMIS_ISOT_TRAC') then
! ------------- TEMPERATURE AU POINT D'INTEGRATION COURANT
                call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                            1, tempg, iret1)
                call rctype(zi(jvMaterc), 1, 'TEMP', [tempg], para_vale, &
                            para_type)
                if ((para_type(1:4) .eq. 'TEMP') .and. (iret1 .eq. 1)) then
                    call utmess('F', 'COMPOR5_5', sk=para_type)
                end if

! ------------- RECUPERATION DE LA COURBE DE TRACTION
                call rctrac(zi(jvMaterc), 1, 'SIGM', para_vale, jprol, &
                            jvale, nbval, e)

! ------------- RECUPERATION DE LA DEFORMATION PLASTIQUE CUMULEE
                p = zr(jvVari+(kpg-1)*nbvari+1-1)

! ------------- TRAVAIL PLASTIQUE 'EQUIVALENT'
                call rcfonc('V', 1, jprol, jvale, nbval, &
                            p=p, rp=rp, rprim=rprim, airerp=airep)
                enerPlas = airep

            end if

! --------- AFFECTATION A OMEGA QUI EST L'ENERGIE TOTALE
! --------- DE LA CONTRIBUTION AU POINT D'INTEGRATION DE L'ENERGIE
! --------- ELASTIQUE ET DE L'ENERGIE PLASTIQUE :
            omega = omega+enerElas+enerPlas

            if (option .eq. 'ENER_TOTALE') then
                if ((relaName(1:9) .ne. 'VMIS_ISOT') .and. (relaName(1:4) .ne. 'ELAS')) then
! ----------------- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION PRECEDENT
! ----------------- ON LE CALCULE SEULEMENT DANS LE CAS DE LOI DE COMPORTEMENT
! ----------------- NI VMIS_ISOT_ NI ELAS
                    if (idconm .ne. 0) then
                        do i = 1, nbsig
                            sigmm(i) = zr(jvSigmPrev+(kpg-1)*nbsig+i-1)
                        end do
                    end if

! ----------------- TENSEUR DES DEFORMATIONS AU POINT D'INTEGRATION COURANT
! ----------------- ON LE CALCULE SEULEMENT DANS LE CAS DE LOI DE COMPORTEMENT
! ----------------- NI VMIS_ISOT_ NI ELAS
                    do i = 1, nbsig
                        epsi(i) = epss(i+(kpg-1)*nbsig)
                    end do

! ----------------- TENSEUR DES DEFORMATIONS AU POINT D'INTEGRATION PRECEDENT
! ----------------- ON LE CALCULE SEULEMENT DANS LE CAS DE LOI DE COMPORTEMENT
! ----------------- NI VMIS_ISOT_ NI ELAS
                    if (ideplm .ne. 0) then
                        do i = 1, nbsig
                            epsim(i) = epssm(i+(kpg-1)*nbsig)
                        end do
                    end if

                    if ((idconm .ne. 0) .and. (ideplm .ne. 0)) then
                        do i = 1, nbsig
                            delta(i) = epsi(i)-epsim(i)
                        end do
! --------------------- CALCUL DES TERMES A SOMMER
                        integ1 = sigmm(1)*delta(1)+ &
                                 sigmm(2)*delta(2)+ &
                                 sigmm(3)*delta(3)+ &
                                 deux*sigmm(4)*delta(4)+ &
                                 deux*sigmm(5)*delta(5)+ &
                                 deux*sigmm(6)*delta(6)
                        integ2 = sigmEner(1)*delta(1)+ &
                                 sigmEner(2)*delta(2)+ &
                                 sigmEner(3)*delta(3)+ &
                                 deux*sigmEner(4)*delta(4)+ &
                                 deux*sigmEner(5)*delta(5)+ &
                                 deux*sigmEner(6)*delta(6)
                        integ = undemi*(integ1+integ2)*poids
                    else
! --------------------- CAS OU LE NUMERO D ORDRE EST UN
                        integ = sigmEner(1)*epsi(1)+ &
                                sigmEner(2)*epsi(2)+ &
                                sigmEner(3)*epsi(3)+ &
                                deux*sigmEner(4)*epsi(4)+ &
                                deux*sigmEner(5)*epsi(5)+ &
                                deux*sigmEner(6)*epsi(6)
                        integ = undemi*integ*poids
                    end if
                    wtotal = wtotal+integ
                else
! ----------------- DANS LE CAS VMIS_ISOT ON CALCULE L ENERGIE TOTALE
                    wtotal = wtotal+(enerElas+enerPlas)*poids
                end if
                goto 10
            end if

! --------- CALCUL DE LA DILATATION VOLUMIQUE AU POINT D'INTEGRATION COURANT
            trepsm = epsiMeca(1+(kpg-1)*nbsig)+epsiMeca(2+(kpg-1)*nbsig)+epsiMeca(3+(kpg-1)*nbsig)

! --------- DEVIATEUR DES DEFORMATIONS AU POINT D'INTEGRATION COURANT
            epsdv(1) = epsiMeca(1+(kpg-1)*nbsig)-untier*trepsm
            epsdv(2) = epsiMeca(2+(kpg-1)*nbsig)-untier*trepsm
            epsdv(3) = epsiMeca(3+(kpg-1)*nbsig)-untier*trepsm
            epsdv(4) = epsiMeca(4+(kpg-1)*nbsig)
            epsdv(5) = epsiMeca(5+(kpg-1)*nbsig)
            epsdv(6) = epsiMeca(6+(kpg-1)*nbsig)
!
! --- CALCUL DE LA DEFORMATION ELASTIQUE EQUIVALENTE AU
! --- POINT D'INTEGRATION COURANT :
!
            epseq = sqrt(trois/deux*(epsdv(1)*epsdv(1)+ &
                                     epsdv(2)*epsdv(2)+ &
                                     epsdv(3)*epsdv(3)+ &
                                     deux*(epsdv(4)*epsdv(4)+ &
                                           epsdv(5)*epsdv(5)+ &
                                           epsdv(6)*epsdv(6))))

! --- COEFFICIENTS DU MATERIAU (CONSTANTE ELASTIQUE DE CISAILLEMENT
! --- DE LAME : MU ET  MODULE ELASTIQUE DE DILATATION : K) :
            deuxmu = e/(un+nu)
            k = untier*e/(un-deux*nu)
!
! --- CALCUL DE LA CONTRAINTE ELASTIQUE EQUIVALENTE AU
! --- POINT D'INTEGRATION COURANT :
!
            sigeq = deuxmu*epseq
!
! --- PARTIE SPHERIQUE DE L'ENERGIE DE DEFORMATION ELASTIQUE :
!
            enelsp = undemi*k*trepsm*trepsm

            if (relaName .eq. 'VMIS_ISOT_LINE') then
! ------------- Parameters from linear traction
                call rcvalb(fami, kpg, ksp, '+', zi(jvMaterc), &
                            ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCodret, 2)
                dsde = propVale(1)
                sigy = propVale(2)

! ------------- PENTE DE LA COURBE DE TRACTION DANS LE DIAGRAMME 'REDRESSE'
                rprim = e*dsde/(e-dsde)

! ------------- DEFORMATION NON-LINEAIRE CUMULEE EQUIVALENTE :
                p = (sigeq-sigy)/(rprim+trois/deux*deuxmu)
                if (p .le. r8prem()) p = zero

! ------------- CONTRAINTE UNIAXIALE SUR LA COURBE DE TRACTION :
                rp = sigy+rprim*p

! ------------- TRAVAIL ELASTIQUE NON-LINEAIRE 'EQUIVALENT' :
                enerPlas = undemi*(sigy+rp)*p

! ------------- TRAITEMENT DU CAS DE L'ECROUISSAGE NON-LINEAIRE ISOTROPE :

            else if (relaName .eq. 'VMIS_ISOT_TRAC') then

! ------------- RECUPERATION DE LA COURBE DE TRACTION :

                call rctrac(zi(jvMaterc), 1, 'SIGM', tempg, jprol, &
                            jvale, nbval, e)

! ------------- CALCUL DE LA LIMITE ELASTIQUE SIGY :
                call rcfonc('S', 1, jprol, jvale, nbval, &
                            sigy=sigy)
!
! --- CALCUL DU TRAVAIL ELASTIQUE NON-LINEAIRE ET DE LA
! --- CONTRAINTE EQUIVALENTE :
!
                call rcfonc('E', 1, jprol, jvale, nbval, &
                            e=e, nu=nu, p=zero, rp=rp, rprim=rprim, &
                            airerp=airep, sieleq=sigeq, dp=p)
!
! --- TRAVAIL ELASTIQUE NON-LINEAIRE 'EQUIVALENT' :
!
                enerPlas = airep
            else
                call utmess('F', 'ENERGY1_6')
            end if
!
! --- PARTIE DEVIATORIQUE DE L'ENERGIE DE DEFORMATION ELASTIQUE
! --- TOTALE 'EQUIVALENTE' (I.E. ASSOCIEE A LA COURBE DE
! --- TRACTION SI ON CONSIDERAIT LE MATERIAU ELASTIQUE
! --- NON-LINEAIRE :
!
            if (p .le. r8prem()) then
                eneldv = epseq*epseq*deuxmu/trois
            else
                eneldv = rp*rp/deuxmu/trois
            end if
!
! --- ENERGIE DE DEFORMATION ELASTIQUE TOTALE AU POINT
! --- D'INTEGRATION COURANT :
!
            enelto = enelsp+eneldv+enerPlas
!
! --- AFFECTATION A PSI QUI EST L'ENERGIE ELASTIQUE TOTALE
! --- DE LA CONTRIBUTION AU POINT D'INTEGRATION DE CETTE ENERGIE :
!
            psi = psi+enelto
!
! --- VOLUME DE L'ELEMENT :
!
            volume = volume+poids
!
! --- INDICATEUR GLOBAL ENERGETIQUE (NON NORMALISE) :
!
            if (omega .ge. 1d04*r8prem()) then
                indigl = indigl+(un-psi/omega)*poids
            end if
!
10          continue
        end do

! ----   RECUPERATION ET AFFECTATION DES GRANDEURS EN SORTIE
! ----   AVEC RESPECTIVEMENT LA VALEUR DE L'INDICATEUR GLOBAL SUR
! ----   L'ELEMENT ET LE VOLUME DE L'ELEMENT POUR L'OPTION
! ----   INDIC_ENER
! ----   AFFECTATION DE L'ENERGIE DE DEFORMATION ELASTIQUE
! ----   ET DE L'ENERGIE DE DEFORMATION TOTALE RESPECTIVEMENT
! ----   POUR LES OPTIONS ENEL_ELEM ET ENER_TOTALE :
        if (option .eq. 'INDIC_ENER') then
            call jevech('PENERD1', 'E', idene1)
            zr(idene1) = indigl
            call jevech('PENERD2', 'E', idene2)
            zr(idene2) = volume
        else if (option .eq. 'ENEL_ELEM') then
            call jevech('PENERD1', 'E', idene1)
            zr(idene1) = welas
        else if (option .eq. 'ENTR_ELEM') then
            call jevech('PENTRD1', 'E', idene1)
            zr(idene1) = welastr
        else if (option .eq. 'ENER_TOTALE') then
            call jevech('PENERD1', 'E', idene1)
            zr(idene1) = wtotal
        end if

    else if (option .eq. 'INDIC_SEUIL') then
        do kpg = 1, npg
! --------- Get elastic parameters
            if (materPara%elasID .ne. ELAS_ISOT) then
                call utmess("F", "ENERGY1_4", nk=2, valk=[option, materPara%elasKeyword])
            end if
            call get_elas_para(fami, zi(jvMaterc), '+', kpg, ksp, &
                               materPara%elasID, materPara%elasKeyword, &
                               e_=e, nu_=nu)
            c1 = (un+nu)/e
            c2 = nu/e

! --------- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT
            sigmEner = 0.d0
            do i = 1, nbsig
                sigmEner(i) = zr(jvSigm+(kpg-1)*nbsig+i-1)
            end do

! --- CALCUL DES DEFORMATIONS ELASTIQUES AU POINT
! --- D'INTEGRATION COURANT EN CONSIDERANT LE MATERIAU ISOTROPE :
! --- EPS_ELAS    = 1/D*SIGMA
! ---             = ((1+NU)/E)*SIGMA-(NU/E)*TRACE(SIGMA) :
            trsig = sigmEner(1)+sigmEner(2)+sigmEner(3)
            epsiElas(1) = c1*sigmEner(1)-c2*trsig
            epsiElas(2) = c1*sigmEner(2)-c2*trsig
            epsiElas(3) = c1*sigmEner(3)-c2*trsig
            epsiElas(4) = c1*sigmEner(4)
            epsiElas(5) = c1*sigmEner(5)
            epsiElas(6) = c1*sigmEner(6)

! --- CALCUL DES DEFORMATIONS PLASTIQUES AU POINT
! --- D'INTEGRATION COURANT
! --- EPS_PLAST = EPS_TOT - EPS_ELAS - EPSTH
! --- EPS_PLAST = EPS_HORS_THERMIQUE - EPS_ELAS
            epsiPlas(1) = epsiMeca(1+(kpg-1)*nbsig)-epsiElas(1)
            epsiPlas(2) = epsiMeca(2+(kpg-1)*nbsig)-epsiElas(2)
            epsiPlas(3) = epsiMeca(3+(kpg-1)*nbsig)-epsiElas(3)
            epsiPlas(4) = epsiMeca(4+(kpg-1)*nbsig)-epsiElas(4)
            epsiPlas(5) = epsiMeca(5+(kpg-1)*nbsig)-epsiElas(5)
            epsiPlas(6) = epsiMeca(6+(kpg-1)*nbsig)-epsiElas(6)

            if (relaName .eq. 'VMIS_CINE_LINE') then
                nbsig2 = 7
                ASSERT(jvVari .ne. 1)
                do i = 1, nbsig
                    x(i) = zr(jvVari+(kpg-1)*nbsig2+i-1)
                end do
                sigmEner = 0.d0
                do i = 1, nbsig
                    sigmEner(i) = sigmEner(i)-x(i)
                end do
            end if

! --------- CALCUL DU TRAVAIL PLASTIQUE AU POINT D'INTEGRATION COURANT
            enerPlas = sigmEner(1)*epsiPlas(1)+ &
                       sigmEner(2)*epsiPlas(2)+ &
                       sigmEner(3)*epsiPlas(3)+ &
                       deux*(sigmEner(4)*epsiPlas(4)+ &
                             sigmEner(5)*epsiPlas(5)+ &
                             sigmEner(6)*epsiPlas(6))

! --------- CALCUL DU TRAVAIL PLASTIQUE EQUIVALENT
! --------- TRAITEMENT DU CAS DE L'ECROUISSAGE LINEAIRE
            if (relaName .eq. 'VMIS_ISOT_LINE' .or. relaName .eq. 'VMIS_CINE_LINE') then
! ------------- Parameters from linear traction
                call rcvalb(fami, kpg, ksp, '+', zi(jvMaterc), &
                            ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCodret, 2)
                dsde = propVale(1)
                sigy = propVale(2)

! ------------- CALCUL DE LA DEFORMATION PLASTIQUE EQUIVALENTE
                epseq = sqrt(trois/deux*( &
                             epsiPlas(1)*epsiPlas(1)+ &
                             epsiPlas(2)*epsiPlas(2)+ &
                             epsiPlas(3)*epsiPlas(3)+ &
                             deux*(epsiPlas(4)*epsiPlas(4)+ &
                                   epsiPlas(5)*epsiPlas(5)+ &
                                   epsiPlas(6)*epsiPlas(6))))

! ------------- DEFORMATION PLASTIQUE CUMULEE (TEMPORAIRE POUR L'ECROUISSAGE CINEMATIQUE)
                if (relaName .eq. 'VMIS_CINE_LINE') then
                    p = epseq
                else if (relaName .eq. 'VMIS_ISOT_LINE') then
                    p = zr(jvVari+(kpg-1)*nbvari+1-1)
                end if

! ------------- PENTE DE LA COURBE DE TRACTION DANS LE DIAGRAMME 'REDRESSE' :
                rprim = e*dsde/(e-dsde)

! ------------- CONTRAINTE UNIAXIALE SUR LA COURBE DE TRACTION :
                rp = sigy+rprim*p

! ------------- TRAVAIL PLASTIQUE 'EQUIVALENT' :
                eplaeq = rp*p

            else if (relaName .eq. 'VMIS_ISOT_TRAC') then
! ------------- RECUPERATION DE LA COURBE DE TRACTION :
                call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                            ksp, tempg, iret1)
                call rctrac(zi(jvMaterc), 1, 'SIGM', tempg, jprol, &
                            jvale, nbval, e)

! ------------- RECUPERATION DE LA DEFORMATION PLASTIQUE CUMULEE
                p = zr(jvVari+(kpg-1)*nbvari+1-1)

! ------------- TRAVAIL PLASTIQUE 'EQUIVALENT'
                call rcfonc('V', 1, jprol, jvale, nbval, &
                            p=p, rp=rp, rprim=rprim, airerp=airep)
                eplaeq = rp*p

            else
                call utmess('F', 'ENERGY1_5')
            end if

! --------- CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT
            call dfdm3d(nno, kpg, jvGaussWeight, jvDBaseFunc, zr(jvGeom), &
                        poids)

! --------- VOLUME DE L'ELEMENT
            volume = volume+poids

! --------- INDICATEUR GLOBAL ENERGETIQUE (NON NORMALISE)
            if (eplaeq .ge. 1d04*r8prem()) then
                indigl = indigl+(un-enerPlas/eplaeq)*poids
            end if
        end do

! ----- RECUPERATION ET AFFECTATION DES GRANDEURS EN SORTIE
! ----- AVEC RESPECTIVEMENT LA VALEUR DE L'INDICATEUR GLOBAL SUR
! ----- L'ELEMENT ET LE VOLUME DE L'ELEMENT
        call jevech('PENERD1', 'E', idene1)
        zr(idene1) = indigl
        call jevech('PENERD2', 'E', idene2)
        zr(idene2) = volume
    end if
!
end subroutine
