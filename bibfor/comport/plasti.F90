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
! aslint: disable=W1504
!
subroutine plasti(BEHInteg, &
                  option, typmod, &
                  fami, kpg, ksp, jvMaterCode, &
                  compor, carcri, instam, instap, &
                  epsdt, depst, &
                  sigm, &
                  nvi, vim, &
                  sigp, vip, &
                  dsidep, codret, &
                  multComp_)
!
    use Behaviour_type
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/get_varc.h"
#include "asterfort/lccnvx.h"
#include "asterfort/lcdedi.h"
#include "asterfort/lcdehy.h"
#include "asterfort/lcelas.h"
#include "asterfort/lcelpl.h"
#include "asterfort/lcmate.h"
#include "asterfort/lcotan.h"
#include "asterfort/lcplas.h"
#include "asterfort/lcpopl.h"
#include "asterfort/lcsmelas.h"
#include "blas/dcopy.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam, instap
    real(kind=8), intent(in) :: epsdt(9), depst(9)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: sigm(6), vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(6), vip(nvi)
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(out) :: dsidep(6, *)
    integer(kind=8), intent(out) :: codret
    character(len=16), optional, intent(in) :: multComp_
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - The PLASTI environment (prefer MFront !)
!
! Main subroutine
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHInteg         : parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES CONTRAINTES           = SIGF(T+DT)
!     CALCUL DES VARIABLES INTERNES    = VINF(T+DT)
!     CALCUL DU JACOBIEN ASSOCIE       = DS/DE(T+DT) OU DS/DE(T)
!     CONVENTION :
!                 SUFFIXE D : DEBUT DU PAS DE TEMPS
!                 SUFFIXE F : FIN DU PAS DE TEMPS
!     ==================================================================
!     ARGUMENTS
!
!     IN FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!        KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!        TYPMOD  TYPE DE MODELISATION
!        IMAT    ADRESSE DU MATERIAU CODE
!        COMP    COMPORTEMENT DE L ELEMENT
!                COMP(1) = RELATION DE COMPORTEMENT (ROUSSELIER.)
!                COMP(2) = NB DE VARIABLES INTERNES
!                COMP(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!        CRIT    CRITERES  LOCAUX
!                CRIT(1) = NOMBRE D ITERATIONS MAXI (ITER_INTE_MAXI)
!                CRIT(3) = TOLERANCE DE CONVERGENCE(RESI_INTE)
!                CRIT(4) = THETA
!                CRIT(5) = ITER_INTE_PAS (UTILISE PAR REDECE EN AMONT)
!                CRIT(6) = ALGO_INTE(NEWTON, NEWTON_PERT, NEWTON_RELI)
!        TIMED   INSTANT T
!        TIMEF   INSTANT T+DT
!        CES PARAMETRES DE TEMPERATURE NE SONT PAS PRIS EN COMPTE EN
!        MECANIQUE PURE (ON UTILISE LES VARIABLES DE COMMANDES)
!
!        EPSDT   DEFORMATION TOTALE A T
!        DEPST   INCREMENT DE DEFORMATION TOTALE
!        SIGD    CONTRAINTE A T
!        VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!        OPT     OPTION DE CALCUL
!                        'RIGI_MECA_TANG'> DSDE(T)
!                        'FULL_MECA'     > DSDE(T+DT), SIGF, VINF
!                        'RAPH_MECA'     > SIGF, VINF
!        ANGMAS  ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
!                +  0 SI NAUTIQUIES OU 2 SI EULER
!                + LES 3 ANGLES D'EULER
!     OUT
!        SIGF    CONTRAINTE A T+DT
!        VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!        DSDE    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!        ICOMP   COMPTEUR POUR LE REDECOUPAGE DU PAS DE TEMPS
!        NVI     NB DE VARIABLES INTERNES
!        IRTETI  CODE RETOUR =0 OK, =1 => REDECOUPAGE DU PAS DE TEMPS
!     ------------------------------------------------------------------
!     INFO    MATERD        (*,1) = CARACTERISTIQUES ELASTIQUES A T
!                           (*,2) = CARACTERISTIQUES PLASTIQUES A T
!             MATERF        (*,1) = CARACTERISTIQUES ELASTIQUES A T+DT
!                           (*,2) = CARACTERISTIQUES PLASTIQUES A T+DT
!             MATCST          'OUI' SI MATERIAU CST ENTRE T ET T+DT
!                             'NON' SINON
!             NDT             NB DE COMPOSANTE TOTALES DES TENSEURS
!                                     = 6  3D
!                                     = 4  AXIS  C_PLAN  D_PLAN
!             NDI             NB DE COMPOSANTE DIRECTES DES TENSEURS
!             NR              NB EQUATION SYSTEME INTEGRE A RESOUDRE
!     ------------------------------------------------------------------
!     ATTENTION
!     SI OPT = 'RIGI_MECA_TANG' NE PAS TOUCHER AUX VARIABLES SIGF,VINF
!     QUI N ONT PAS DE PLACE MEMOIRE ALLOUEE
!
!     SIG EPS DEPS  ONT DEJA LEURS COMPOSANTES DE CISAILLEMENT
!     MULTIPLIES PAR RACINE DE 2 > PRISE EN COMPTE DES DOUBLES
!     PRODUITS TENSORIELS ET CONSERVATION DE LA SYMETRIE
!
!     ------------------------------------------------------------------
!
!     NMAT = NOMBRE MAXI DE PARMETRES MATERIAU
!     POUR LE MONOCRISTAL, DIMENSIONS MAX
!     NSG=NOMBRE DE SYSTEMES DE GLISSEMENT MAXIMUM
!     NFS=NOMBRE DE FAMILLES DE SYSTEMES DE GLISSEMENT MAXIMUM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: rungeKutta = 0
    integer(kind=8), parameter :: nmat = 90, nsg = 30, nfs = 5, nrm = nfs*nsg+6
    real(kind=8), allocatable :: drdy(:)
    character(len=3) :: matcst
    character(len=7) :: etatd, etatf
    character(len=8) :: typmod1, typma
    character(len=16) :: relaComp, defoComp, multComp
    character(len=24) :: cpmono(5*nmat+1)
    aster_logical :: l_temp
    integer(kind=8) :: ndt, ndi, nr, itmax, irtet, cutLevel
    integer(kind=8) :: nbcomm(nmat, 3), numhsr(1), irr, decirr, nbsyst, decal, gdef
    real(kind=8) :: toler, epsi, materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: epsd(9), deps(9)
    real(kind=8) :: seuil, theta, dt, devg(6), devgii
    real(kind=8) :: vp(3), vecp(3, 3), pgl(3, 3)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
    real(kind=8) :: tempd, tempf, tref
!     POUR BETON_BURGER - ATTENTION DIMENSION MAXI POUR CE MODELE
    parameter(epsi=1.d-15)
    aster_logical :: resi, rigi
    blas_int :: b_incx, b_incy, b_n
    common/tdim/ndt, ndi
    common/polycr/irr, decirr, nbsyst, decal, gdef
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = BEHInteg%materPara
    cutLevel = BEHInteg%behavPara%cutLevel
    codret = 0
    itmax = int(carcri(1))
    toler = carcri(3)
    theta = carcri(4)
    relaComp = compor(RELA_NAME)
    defoComp = compor(DEFO)
    multComp = ' '
    if (present(multComp_)) then
        multComp = multComp_
    end if
    typmod1 = typmod(1)
    dt = instap-instam
    resi = option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA'
    rigi = option(1:9) .eq. 'RIGI_MECA' .or. option(1:9) .eq. 'FULL_MECA'
    gdef = 0
    if (defoComp .eq. 'SIMO_MIEHE') then
        gdef = 1
    end if
    numhsr(1) = 1
!
    typma = 'VITESSE '

! - Allocate big objects
    allocate (drdy(nrm*nrm))

! - Get temperatures
    call get_varc(fami, kpg, ksp, 'T', tempd, &
                  tempf, tref, l_temp)

! - Glute pour LKR
    if (.not. l_temp .and. relaComp .eq. 'LKR') then
        tempd = 0.d0
        tempf = 0.d0
        tref = 0.d0
    end if

! - RECUPERATION COEF MATERIAU A T ET/OU T+DT
    call lcmate(materPara, &
                carcri, relaComp, typmod1, &
                nmat, tempd, tempf, &
                tref, rungeKutta, typma, hsr, materd, materf, matcst, &
                nbcomm, cpmono, pgl, itmax, &
                toler, ndt, ndi, nr, &
                nvi, vim, nfs, nsg, toutms, &
                1, numhsr, sigm, multComp)

    if (gdef .eq. 1) then
!        GDEF_MONO : PAS DE DEFORM. THERMIQUE
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, depst, b_incx, deps, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsdt, b_incx, epsd, b_incy)
    else
! --     RETRAIT INCREMENT DE DEFORMATION DUE A LA DILATATION THERMIQUE
        call lcdedi(fami, kpg, ksp, nmat, materd, &
                    materf, tempd, tempf, tref, depst, &
                    epsdt, deps, epsd)
! --     RETRAIT ENDOGENNE ET RETRAIT DE DESSICCATION
        call lcdehy(fami, kpg, ksp, nmat, materd, &
                    materf, deps, epsd)
    end if
!
! --    SEUIL A T > ETAT ELASTIQUE OU PLASTIQUE A T
    if (abs(vim(nvi)) .le. epsi) then
        etatd = 'ELASTIC'
    else
        etatd = 'PLASTIC'
    end if
!
!
    if (option(1:10) .eq. 'RIGI_MECA_' .and. gdef .eq. 1 .and. relaComp .eq. 'MONOCRISTAL') then
        call lcsmelas(epsd, deps, dsidep, nmat=nmat, materd_=materd)
        codret = 0
        goto 999
    end if
!
!     ----------------------------------------------------------------
!     OPTIONS 'FULL_MECA' ET 'RAPH_MECA' = CALCUL DE SIG(T+DT)
!     ----------------------------------------------------------------
!
    if (resi) then
!
        if (gdef .eq. 1) then
!           GDEF_MONO : PAS DE SEUIL CAR C'EST PLUS COMPLIQUE
            seuil = 1.d0
        else
! --        INTEGRATION ELASTIQUE SUR DT
            call lcelas(relaComp, typmod1, &
                        nmat, materd, materf, matcst, &
                        deps, sigm, vim, &
                        sigp, theta)
!
! --        PREDICTION ETAT ELASTIQUE A T+DT : F(SIG(T+DT),VIN(T)) = 0 ?
            seuil = 1.d0
            call lccnvx(fami, kpg, ksp, relaComp, &
                        jvMaterCode, nmat, materf, sigm, sigp, &
                        deps, vim, vip, nbcomm, cpmono, &
                        pgl, nvi, vp, vecp, hsr, &
                        nfs, nsg, toutms, instam, instap, &
                        seuil)
!
        end if
!
        if (seuil .ge. 0.d0) then
! --        PREDICTION INCORRECTE > INTEGRATION ELASTO-PLASTIQUE SUR DT
            etatf = 'PLASTIC'
!
            call lcplas(BEHInteg, &
                        fami, kpg, ksp, relaComp, &
                        toler, itmax, typmod1, jvMaterCode, nmat, &
                        materd, materf, nr, nvi, instam, &
                        instap, deps, epsd, sigm, vim, &
                        sigp, vip, compor, nbcomm, cpmono, &
                        pgl, nfs, nsg, toutms, hsr, &
                        irtet, theta, vp, &
                        seuil, devg, devgii, drdy, carcri)
!
!
            if (irtet .eq. 1) then
                goto 1
            else if (irtet .eq. 2) then
                goto 2
            end if
        else
! --        PREDICTION CORRECTE > INTEGRATION ELASTIQUE FAITE
            etatf = 'ELASTIC'
! ---       MISE A JOUR DE VINF EN FONCTION DE LA LOI
!           ET POST-TRAITEMENTS POUR DES LOIS PARTICULIERES
            call lcelpl(relaComp, nmat, materf, deps, nvi, &
                        vim, vip)
        end if
!
!        POST-TRAITEMENTS PARTICULIERS
        call lcpopl(relaComp, nmat, materd, materf, &
                    typmod1, sigp, vim, &
                    vip)
!
    end if
!
!     ----------------------------------------------------------------
!     OPTIONS 'FULL_MECA' ET 'RIGI_MECA_TANG' = CALCUL DE DSDE
!     ----------------------------------------------------------------
!     EVALUATION DU JACOBIEN DSDE A (T+DT) POUR 'FULL_MECA'
!     ET CALCUL ELASTIQUE    ET   A (T)    POUR 'RIGI_MECA_TANG'
!     ----------------------------------------------------------------
!
    if (rigi) then
        call lcotan(option, etatd, etatf, fami, &
                    kpg, ksp, relaComp, typmod1, jvMaterCode, &
                    nmat, materd, materf, epsd, deps, &
                    sigm, sigp, nvi, vim, vip, &
                    drdy, vp, vecp, theta, dt, &
                    devg, devgii, instam, instap, compor, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, nr, itmax, toler, &
                    typma, dsidep, irtet)
        if (irtet .ne. 0) goto 1
!
    end if
!
!       ----------------------------------------------------------------
!
    codret = 0
    goto 999
1   continue
    codret = 1
    goto 999
!
2   continue
    codret = 2
    goto 999
!
999 continue

! - De-Allocate big objects
    deallocate (drdy)
!
end subroutine
