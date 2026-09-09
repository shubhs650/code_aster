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
! aslint: disable=W1504,W0413
!
subroutine pminit(tablName, tablNbPara, tablType, &
                  tablParaName, tablParaType, tablVale, &
                  pgl, lRota, &
                  epsiPrev, sigmPrev, &
                  nbVari, vim, vip, &
                  loadEpsiType, loadType, loadFunc, coefImpo, &
                  coefAdim, typeMatrPred, lMatrElas, matrElas, lPrintMatr, option, &
                  variName, nbVariTabl, &
                  sddisc, ds_conv, ds_algopara, sderro, materPara)
!
    use MaterialPara_module
    use MaterialPara_type
    use NonLin_Datastructure_type
    implicit none
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/diinst.h"
#include "asterfort/dmat3d.h"
#include "asterfort/fointe.h"
#include "asterfort/fozero.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/nmcrga.h"
#include "asterfort/nmcrli.h"
#include "asterfort/nmcrsu.h"
#include "asterfort/nmdocn.h"
#include "asterfort/nmdomt.h"
#include "asterfort/nonlinDSAlgoParaCreate.h"
#include "asterfort/nonlinDSConvergenceCreate.h"
#include "asterfort/nonlinDSConvergenceInit.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcinp.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbVari
    character(len=8), intent(out) :: tablName
    integer(kind=8), intent(out) :: tablNbPara, tablType
    character(len=16), allocatable, intent(out) :: tablParaName(:)
    character(len=8), allocatable, intent(out) :: tablParaType(:)
    real(kind=8), allocatable, intent(out) :: tablVale(:)
    real(kind=8), intent(out) :: pgl(3, 3)
    aster_logical, intent(out) :: lRota
    real(kind=8), intent(out) :: epsiPrev(9), sigmPrev(6)
    real(kind=8), intent(out) :: vim(nbVari), vip(nbVari)
    integer(kind=8), intent(out) :: loadEpsiType, loadType(9)
    character(len=8), intent(out) :: loadFunc(9)
    real(kind=8), intent(out) :: coefImpo(6, 12), coefAdim
    integer(kind=8), intent(out) :: typeMatrPred
    aster_logical, intent(out) :: lMatrElas
    real(kind=8), intent(out) :: matrElas(6, 6)
    aster_logical, intent(out) :: lPrintMatr
    character(len=16), intent(out) :: option
    character(len=8), intent(out) :: variName(nbVari)
    integer(kind=8), intent(out) :: nbVariTabl
    character(len=19), intent(out) :: sddisc
    type(NL_DS_Conv), intent(inout) :: ds_conv
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    character(len=24), intent(out) :: sderro
    type(Material_Para), intent(inout) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
! SIMU_POINT_MAT
!
! Initializations
!
! --------------------------------------------------------------------------------------------------
!
! IN   IMATE  : ADRESSE MATERIAU CODE
! IN   NBVARI : NOMBRE DE VARIABLES INTERNES
! IN   NDIM   : 3
! IO  ds_conv          : datastructure for convergence management
! IO  ds_algopara      : datastructure for algorithm parameters
! OUT  TYPMOD : 3D
! OUT  TABLE  : TABLE RESULTAT
! OUT  NBPAR  : NOMBRE DE PARAMETRES DE LA TABLE RESULTAT
! OUT  NOMPAR : NOMS DES PARAMETRES DE LA TABLE RESULTAT
! OUT  ANGL_NAUT : ANGLES DU MOT-CLE MASSIF
! OUT  PGL    : MATRICE DE ROTATION AUTOUR DE Z
! OUT  IROTA  : =1 SI ROTATION AUTOUR DE Z
! OUT  EPSM   : DEFORMATIONS INITIALES
! OUT  SIGM   : CONTRAINTES INITIALES
! OUT  VIM    : VARIABLES INTERNES INITIALES
! OUT  VIP    : VARIABLES INTERNES NULLES
! OUT  DEFIMP : =1 SI LES 6 CMP DE EPSI DONT DONNEES
! OUT  COEF   : COEF POUR ADIMENSIONNALISER LE PB
! OUT  INDIMP : TABLEAU D'INDICES =1 SI EPS(I) DONNE
! OUT  FONIMP : FONCTIONS IMPOSEES POUR EPSI OU SIGM
! OUT  CIMPO  : = 1 POUR LA CMP DE EPSI OU SIGM IMPOSEE
! OUT  KEL    : OPERATEUR D'ELASTICITE
! OUT  SDDISC : SD DISCRETISATION
! OUT  PRED   : TYPE DE PREDICTION = 1 SI TANGENTE
! OUT  MATREL : MATRICE TANGENTE = 1 SI ELASTIQUE
! OUT  OPTION : FULL_MECA OU RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: tablNbParaMaxi = 9999
    integer(kind=8), parameter :: ndim = 3
    real(kind=8), parameter  :: rac2 = sqrt(2.d0)
    complex(kind=8), parameter :: c16Dummy = (0.d0, 0.d0)
    character(len=4), parameter :: epsiName(6) = (/'EPXX', 'EPYY', 'EPZZ', &
                                                   'EPXY', 'EPXZ', 'EPYZ'/)
    character(len=4), parameter :: sigmName(6) = (/'SIXX', 'SIYY', 'SIZZ', &
                                                   'SIXY', 'SIXZ', 'SIYZ'/)
    character(len=4), parameter :: gradName(9) = (/'F11', 'F12', 'F13', &
                                                   'F21', 'F22', 'F23', &
                                                   'F31', 'F32', 'F33'/)
    real(kind=8), parameter :: id(9) = (/1.d0, 0.d0, 0.d0, &
                                         0.d0, 1.d0, 0.d0, &
                                         0.d0, 0.d0, 1.d0/)
    character(len=8), parameter :: f0 = '&&CPM_F0'
    integer(kind=8) :: n1, i, j, k, nbocc
    integer(kind=8) :: nbEpsiLoad, nbSigmLoad, nbGradLoad, nbCmpEpsi, nbVariInit
    integer(kind=8) :: iligne, icolon, nbcol, numeInst, ier
    character(len=8) :: k8b, gradValue(9), vk8(2), answer
    character(len=8) :: epsiVale(6), sigmVale(6), valef
    character(len=16) :: tablFormat, cmdName, dsType
    character(len=19) :: listInst
    real(kind=8) :: timePrev, vale, timeInit
    real(kind=8) :: sigmInitVale, epsiInitVale
    real(kind=8) :: ang1(1), dsidep(36)
    real(kind=8) :: sigmInit(6), epsiInit(6), gradInitImpo(9)
    aster_logical :: lSetLinearRela, lGrad
!
! --------------------------------------------------------------------------------------------------
!
    tablName = " "
    tablNbPara = 0
    sddisc = '&&OP0033.SDDISC'
    sderro = '&&OP0033.ERRE.'

! - Create zero function
    call fozero(f0)

! - Create convergence management datastructure
    call nonlinDSConvergenceCreate(ds_conv)

! - Create algorithm parameters datastructure
    call nonlinDSAlgoParaCreate(ds_algopara)

! - Read parameters for convergence
    call nmdocn(ds_conv)

! - Read parameters for algorithm management
    call nmdomt(ds_algopara)

! - Create datastructure for events in algorithm
    call nmcrga(sderro)
!
! - Initializations for convergence management
    call nonlinDSConvergenceInit(ds_conv, sderro)

! - Get output tablName and its parameters
    call getres(tablName, dsType, cmdName)
    ASSERT(dsType .eq. 'TABLE_SDASTER')
    ASSERT(cmdName .eq. 'CALC_POINT_MAT')
    tablType = 0
    call getvtx(' ', 'FORMAT_TABLE', scal=tablFormat, nbret=n1)
    if (n1 .ne. 0) then
        if (tablFormat .eq. 'CMP_LIGNE') then
            tablType = 1
        end if
    end if

! - Count number of variables in table
    nbVariTabl = nbVari
    call getvis(' ', 'NB_VARI_TABLE', scal=k, nbret=n1)
    if (n1 .gt. 0) then
        nbVariTabl = k
    end if
    nbVariTabl = min(nbVariTabl, nbVari)

! - Number of components for strain
    nbCmpEpsi = 6
    lGrad = ASTER_FALSE
    call getvid(' ', gradName(1), scal=gradValue(1), nbret=n1)
    if (n1 .ne. 0) then
        nbCmpEpsi = 9
        lGrad = ASTER_TRUE
    end if

! - Print tangent operator ?
    lPrintMatr = ASTER_FALSE
    call getvtx(' ', 'OPER_TANGENT', scal=answer, nbret=n1)
    if (n1 .ne. 0) then
        if (answer .eq. 'OUI') then
            lPrintMatr = ASTER_TRUE
            ASSERT(.not. lGrad)
        end if
    end if

! - Size of table
    nbcol = 1+nbCmpEpsi+6+2+nbVariTabl+1+36
    if (nbcol .gt. tablNbParaMaxi .and. tablType .eq. 0) then
        call utmess('F', 'COMPOR1_68', si=nbcol)
    end if

! - Number of parameters in table
    if (tablType .eq. 0) then
        tablNbPara = 1+nbCmpEpsi+6+2+nbVariTabl+1
        if (lPrintMatr) then
            tablNbPara = tablNbPara+36
        end if
    else
        tablNbPara = 4
    end if

! - Allocate objects
    allocate (tablParaType(tablNbPara))
    allocate (tablParaName(tablNbPara))
    allocate (tablVale(tablNbPara))
    tablParaType(1:tablNbPara) = " "
    tablParaName(1:tablNbPara) = " "
    tablVale(1:tablNbPara) = 0.d0

! - Names of parameters in table
    tablParaName(1) = 'INST'
    if (tablType .eq. 0) then
        if (lGrad) then
            do i = 1, nbCmpEpsi
                tablParaName(1+i) = gradName(i)
            end do
        else
            do i = 1, nbCmpEpsi
                tablParaName(1+i) = epsiName(i)
            end do
        end if
        do i = 1, 6
            tablParaName(1+nbCmpEpsi+i) = sigmName(i)
        end do
        tablParaName(1+nbCmpEpsi+6+1) = 'TRACE'
        tablParaName(1+nbCmpEpsi+6+2) = 'VMIS'
        do i = 1, nbVariTabl
            tablParaName(1+nbCmpEpsi+6+2+i) (1:1) = 'V'
            call codent(i, 'G', tablParaName(1+nbCmpEpsi+6+2+i) (2:16))
        end do
        if (lPrintMatr) then
            do i = 1, 6
                do j = 1, 6
                    k = 1+nbCmpEpsi+6+2+nbVari+6*(i-1)+j
                    write (tablParaName(k), '(A,I1,I1)') 'K', i, j
                end do
            end do
        end if
        tablParaName(tablNbPara) = 'NB_ITER'
        tablParaType(1:tablNbPara) = 'R'
    else
        tablParaName(2) = 'GRANDEUR'
        tablParaName(3) = 'CMP'
        tablParaName(4) = 'VALEUR'
        tablParaType(1) = 'R'
        tablParaType(2) = 'K8'
        tablParaType(3) = 'K8'
        tablParaType(4) = 'R'
    end if

! - Create table and set list of parameters
    call tbcrsd(tablName, 'G')
    call tbajpa(tablName, tablNbPara, tablParaName, tablParaType)

! - Get local coordinate system
    call getUserLCSCommand(ndim, materPara%lcsPara)

    if (nbCmpEpsi .eq. 6) then
        epsiPrev = 0.d0
    else
        epsiPrev = id
    end if

! - ANGLE DE ROTATION
    lRota = ASTER_FALSE
    pgl = 0.d0
    call getvr8(' ', 'ANGLE', scal=ang1(1), nbret=n1)
    if ((n1 .ne. 0) .and. (ang1(1) .ne. 0.d0)) then
        lRota = ASTER_TRUE
        ang1(1) = ang1(1)*r8dgrd()
        pgl(1, 1) = cos(ang1(1))
        pgl(2, 2) = cos(ang1(1))
        pgl(1, 2) = sin(ang1(1))
        pgl(2, 1) = -sin(ang1(1))
        pgl(3, 3) = 1.d0
    end if

! - Get initial state - Stresses
    sigmPrev = 0.d0
    call getfac('SIGM_INIT', nbocc)
    if (nbocc .gt. 0) then
        do i = 1, 6
            call getvr8('SIGM_INIT', sigmName(i), iocc=1, scal=sigmInitVale, nbret=n1)
            if (n1 .ne. 0) then
                sigmPrev(i) = sigmInitVale
            end if
        end do
        sigmPrev(4:6) = sigmPrev(4:6)*rac2
    end if

! - Get initial state - Strains
    epsiPrev = 0.d0
    if (lGrad) then
        ASSERT(nbCmpEpsi .eq. 9)
        epsiPrev = id
    else
        ASSERT(nbCmpEpsi .eq. 6)
        epsiPrev = 0.d0
    end if
    call getfac('EPSI_INIT', nbocc)
    if (nbocc .gt. 0) then
        do i = 1, 6
            call getvr8('EPSI_INIT', epsiName(i), iocc=1, scal=epsiInitVale, nbret=n1)
            if (n1 .ne. 0) then
                epsiPrev(i) = epsiInitVale
            end if
        end do
        epsiPrev(4:6) = epsiPrev(4:6)*rac2
    end if

! - Get initial state - Internal state variables
    vim = 0.d0
    vip = 0.d0
    call getfac('VARI_INIT', nbocc)
    if (nbocc .gt. 0) then
        call getvr8('VARI_INIT', 'VALE', iocc=1, nbval=nbVari, vect=vim, nbret=nbVariInit)
        if (nbVariInit .ne. nbVari) then
            call utmess('F', 'COMPOR1_72', ni=2, vali=[nbVariInit, nbVari])
        end if
    end if

! - Get loads (strain or stress ?)
    coefImpo = 0.d0
    loadType = 0
    loadFunc = f0
    nbEpsiLoad = 0
    do i = 1, 6
        call getvid(' ', epsiName(i), scal=epsiVale(i), nbret=n1)
        if (n1 .ne. 0) then
            coefImpo(i, 6+i) = 1.d0
            loadFunc(i) = epsiVale(i)
            nbEpsiLoad = nbEpsiLoad+1
            loadType(i) = 1
        end if
    end do
    nbSigmLoad = 0
    do i = 1, 6
        call getvid(' ', sigmName(i), scal=sigmVale(i), nbret=n1)
        if (n1 .ne. 0) then
            coefImpo(i, i) = 1.d0
            loadFunc(i) = sigmVale(i)
            nbSigmLoad = nbSigmLoad+1
            loadType(i) = 0
        end if
    end do
    nbGradLoad = 0
    do i = 1, 9
        call getvid(' ', gradName(i), scal=gradValue(i), nbret=n1)
        if (n1 .ne. 0) then
            loadFunc(i) = gradValue(i)
            nbGradLoad = nbGradLoad+1
            loadType(i) = 2
        end if
    end do

! - Set type of strain load
    loadEpsiType = 0
    if (nbEpsiLoad .eq. 6) loadEpsiType = 1
    if (nbGradLoad .eq. 9) loadEpsiType = 2

! - Linear relations
    lSetLinearRela = ASTER_FALSE
    call getfac('MATR_C1', nbocc)
    if (nbocc .ne. 0) then
        lSetLinearRela = ASTER_TRUE
        do i = 1, nbocc
            call getvis('MATR_C1', 'NUME_LIGNE', iocc=i, scal=iligne, nbret=n1)
            call getvis('MATR_C1', 'NUME_COLONNE', iocc=i, scal=icolon, nbret=n1)
            call getvr8('MATR_C1', 'VALE', iocc=i, scal=vale, nbret=n1)
            coefImpo(iligne, icolon) = vale
        end do
    end if

    call getfac('MATR_C2', nbocc)
    if (nbocc .ne. 0) then
        lSetLinearRela = ASTER_TRUE
        do i = 1, nbocc
            call getvis('MATR_C2', 'NUME_LIGNE', iocc=i, scal=iligne, nbret=n1)
            call getvis('MATR_C2', 'NUME_COLONNE', iocc=i, scal=icolon, nbret=n1)
            call getvr8('MATR_C2', 'VALE', iocc=i, scal=vale, nbret=n1)
            coefImpo(iligne, icolon+6) = vale
        end do
    end if
    call getfac('VECT_IMPO', nbocc)
    if (nbocc .ne. 0) then
        do i = 1, nbocc
            call getvis('VECT_IMPO', 'NUME_LIGNE', iocc=i, scal=iligne, nbret=n1)
            call getvid('VECT_IMPO', 'VALE', iocc=i, scal=valef, nbret=n1)
            loadFunc(iligne) = valef
        end do
    end if
    if (lSetLinearRela) then
        do i = 1, 6
            k = 0
            do j = 1, 12
                if (coefImpo(i, j) .ne. 0.d0) then
                    k = 1
                end if
            end do
            if (k .eq. 0) then
                coefImpo(i, i) = 1.d0
            end if
        end do
        loadEpsiType = -1
    end if

!  RECUPERER LES VALEURS INITIALES DE F "GRAD_IMPOSE"
    gradInitImpo = 0.d0
    timeInit = 0.d0
    if (loadEpsiType == 2) then
        ASSERT(nbGradLoad .eq. 9)
        do i = 1, 9
            call fointe('F', loadFunc(i), 1, ['INST'], [timeInit], gradInitImpo(i), ier)
        end do
    end if

! - Set initial state in table
    if (tablType .eq. 0) then
        tablVale(1) = timeInit
        if (loadEpsiType == 2) then
            ASSERT(nbGradLoad .eq. 9)
            tablVale(2:2+nbCmpEpsi-1) = gradInitImpo(1:nbCmpEpsi)
        else
            epsiInit(1:6) = epsiPrev(1:6)
            epsiInit(4:6) = epsiInit(4:6)*1.d0/rac2
            tablVale(2:2+nbCmpEpsi-1) = epsiInit(1:nbCmpEpsi)
        end if
        sigmInit(1:6) = sigmPrev(1:6)
        sigmInit(4:6) = sigmInit(4:6)*1.d0/rac2
        tablVale(2+nbCmpEpsi:2+nbCmpEpsi+6-1) = sigmInit(1:6)
        tablVale(2+nbCmpEpsi+6) = 0.d0
        tablVale(2+nbCmpEpsi+7) = 0.d0
        tablVale(2+nbCmpEpsi+8:2+nbCmpEpsi+8+nbVariTabl-1) = vim(1:nbVariTabl)

        if (lPrintMatr) then
            dsidep = 0.d0
            tablVale(2+nbCmpEpsi+8+nbVariTabl:2+nbCmpEpsi+8+nbVariTabl+36-1) = dsidep(1:36)
        end if
        tablVale(tablNbPara) = 0.d0
        call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, [c16Dummy], k8b, 0)
    else

        vk8(1) = 'EPSI'
        do i = 1, nbCmpEpsi
            tablVale(2) = epsiPrev(i)
            vk8(2) = epsiName(i)
            call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, [c16Dummy], vk8, 0)
        end do
        vk8(1) = 'SIGM'
        do i = 1, nbCmpEpsi
            tablVale(2) = sigmPrev(i)
            vk8(2) = sigmName(i)
            call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, [c16Dummy], vk8, 0)
        end do
        vk8(1) = 'VARI'
        do i = 1, nbVariTabl
            tablVale(2) = vim(i)
            vk8(2) (1:1) = 'V'
            call codent(i, 'G', vk8(2) (2:8))
            variName(i) = vk8(2)
            call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, [c16Dummy], vk8, 0)
        end do
    end if

! - Create time discretization datastructure
    call getvid('INCREMENT', 'LIST_INST', iocc=1, scal=listInst, nbret=n1)
    call nmcrli(listInst, sddisc)

! - First time step
    numeInst = 0
    timePrev = diinst(sddisc, numeInst)

! - Get parameters for Newton algorithm
    lMatrElas = ASTER_FALSE
    option = 'FULL_MECA'
    if (ds_algopara%matrix_corr .eq. 'ELASTIQUE') then
        lMatrElas = ASTER_TRUE
        typeMatrPred = 0
        option = 'RAPH_MECA'
    end if
    typeMatrPred = 1
    if (ds_algopara%matrix_pred .eq. 'ELASTIQUE') then
        typeMatrPred = 0
    else if (ds_algopara%matrix_pred .eq. 'EXTRAPOLE') then
        typeMatrPred = -1
    end if
! - Automatic management of time stepping
    call nmcrsu(sddisc, listInst, ds_conv)

! - CALCUL DES VARIABLES DE COMMANDE
    call vrcinp(2, timePrev, timePrev)

! - Compute elastic matrix
    matrElas = 0.d0
    call dmat3d(materPara, '+', timePrev, matrElas)

!     DMAT ECRIT MU POUR LES TERMES DE CISAILLEMENT
    coefAdim = max(matrElas(1, 1), matrElas(2, 2), matrElas(3, 3))
    do j = 4, 6
        matrElas(j, j) = matrElas(j, j)*2.d0
        coefAdim = max(coefAdim, matrElas(j, j))
    end do
    if (lSetLinearRela) then
        coefAdim = 1.d0
    end if
!
end subroutine
