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
subroutine op0033()
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/detrsd.h"
#include "asterfort/diinst.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveut.h"
#include "asterfort/lcdetf.h"
#include "asterfort/matinv.h"
#include "asterfort/mgauss.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmcrcv.h"
#include "asterfort/nmfinp.h"
#include "asterfort/pmactn.h"
#include "asterfort/pmadat.h"
#include "asterfort/pmconv.h"
#include "asterfort/pmdocc.h"
#include "asterfort/pmdocr.h"
#include "asterfort/pmdrdy.h"
#include "asterfort/pmimpr.h"
#include "asterfort/pminit.h"
#include "asterfort/pmmaco.h"
#include "asterfort/pmsta1.h"
#include "asterfort/pmstab.h"
#include "asterfort/pmvtgt.h"
#include "asterfort/tnsvec.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcinp.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! CALC_POINT_MAT
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: typmod(2) = (/"3D", "  "/)
    integer(kind=8), parameter :: ndim = 3, ksp = 1, kpg = 1
    character(len=4), parameter :: fami = "PMAT"
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    character(len=8), allocatable :: tablParaType(:)
    character(len=16), allocatable :: tablParaName(:)
    real(kind=8), allocatable :: tablVale(:)
    integer(kind=8) :: tablNbPara, tablType
    integer(kind=8) :: iret, nbVari, i, ier
    integer(kind=8) :: jvMaterCode, iterNewt, ncmp
    integer(kind=8) :: loadEpsiType, liccvg(5)
    integer(kind=8) :: loadType(9), numeInst, newtLoopAction, action, itgt
    integer(kind=8) :: nbVariTabl, typeMatrPred, typeMatrTang, nbMaterDefi
    character(len=4) :: cargau
    character(len=8) :: tablName, loadFunc(9), materDefi(30)
    character(len=16) :: comporList(COMPOR_SIZE)
    character(len=16) :: option, optionNoNewton, optionPred
    character(len=16) :: relaComp, defoComp, multComp
    character(len=19) :: codi
    real(kind=8) :: timePrev, timeCurr, r8b, carcriList(CARCRI_SIZE), fem(9)
    real(kind=8) :: epsiIncr(9), sigmPrev(6), sigmCurr(6), epsiPrev(9)
    real(kind=8) :: valeImpo(9), resi(12), resiInit(12), dy(12), ddy(12), y(12)
    real(kind=8) :: dsidep(6, 9), drdy(12, 12), matrElas(6, 6), coefImpo(6, 12), ym(12)
    real(kind=8) :: work(10), sdeps(6), sigmCurrPred(6), smatr(36), resiEval(12)
    real(kind=8) :: matper(36), varia(2*36), epsilo, pgl(3, 3), matrRotaLoad(3, 3)
    real(kind=8) :: vectRotaLoad(3, 3), coefAdim, jm, jp, jd, coefextra
    aster_logical :: lastTimeStep, lIterNewtMaxi, conver, lPrintMatr, lMatrElas, lRota, lLoadGrad
    aster_logical :: lNewtonAlgo
    integer(kind=8) :: jvVim, jvVip, jvVimCopy, jvVipCopy, jvVariName
    type(NL_DS_Conv) :: ds_conv
    type(NL_DS_AlgoPara) :: ds_algopara
    type(Behaviour_Integ) :: BEHinteg
    blas_int :: b_incx, b_incy, b_n
    character(len=19) :: sddisc
    character(len=19), parameter :: sdcrit = '&&OP0033.SDCRIT'
    character(len=24) :: sderro
    character(len=19), parameter :: variNameJv = '&&OP0033.NOMVI'
    character(len=19), parameter :: vimJvName = '&&OP0033.VIM', vipJvName = '&&OP0033.VIP'
    character(len=19), parameter :: vipCopyJvName = '&&OP0033.SVIP'
    character(len=19), parameter :: vimCopyJvName = '&&OP0033.VIM2'
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call infmaj()
    call jemarq()

! - Initializations
    work = 0.d0
    dsidep = 0.d0
    iterNewt = 0
    action = 1
    lastTimeStep = ASTER_FALSE
    lIterNewtMaxi = ASTER_FALSE
    liccvg = 0

! - Prepare CALCUL parameters for external state variables
    call vrcinp(1, 0.d0, 0.d0)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Get material parameters
    call getvid(' ', 'MATER', nbval=6, vect=materDefi, nbret=nbMaterDefi)
!
! - Get list of parameters for constitutive law
    call pmdocc(comporList, nbVari, relaComp, defoComp, multComp)
!
! - Get list of parameters for integration of constitutive law
    call pmdocr(carcriList)

! - Create working vectors
    call wkvect(vimJvName, 'V V R', nbVari, jvVim)
    call wkvect(vipJvName, 'V V R', nbVari, jvVip)
    call wkvect(vipCopyJvName, 'V V R', nbVari, jvVipCopy)
    call wkvect(vimCopyJvName, 'V V R', nbVari, jvVimCopy)
    call wkvect(variNameJv, 'V V K8', nbVari, jvVariName)

! - Coding material parameters
    call pmmaco(materDefi, nbMaterDefi, codi)
    call jeveut(codi//'.CODI', 'L', jvMaterCode)

! - Initializations of material parameters
    call initParaCell(fami, jvMaterCode, materPara)
    call initParaPoin(kpg, ksp, materPara)

! - Initializations
    call pminit(tablName, tablNbPara, tablType, &
                tablParaName, tablParaType, tablVale, &
                pgl, lRota, &
                epsiPrev, sigmPrev, &
                nbVari, zr(jvVim), zr(jvVip), &
                loadEpsiType, loadType, loadFunc, coefImpo, coefAdim, &
                typeMatrPred, lMatrElas, matrElas, lPrintMatr, option, &
                zk8(jvVariName), nbVariTabl, &
                sddisc, ds_conv, ds_algopara, sderro, materPara)

! - Message if PETIT_REAC
    if (loadEpsiType .gt. 0) then
        if (defoComp .eq. 'PETIT_REAC') then
            call utmess('I', 'COMPOR2_93')
        end if
    end if

! - CREATION DE LA SD POUR ARCHIVAGE DES INFORMATIONS DE CONVERGENCE
    call nmcrcv(sdcrit)
    numeInst = 1
!
!==================================
!     BOUCLE SUR lES INSTANTS
!==================================
!
200 continue
!
    liccvg = 0

! - Get times
    timePrev = diinst(sddisc, numeInst-1)
    timeCurr = diinst(sddisc, numeInst)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              comporList, carcriList, &
                              timePrev, timeCurr, &
                              materPara, BEHInteg)

! - Set main parameters for behaviour (on point)
    call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! - Compute external state variables
    call vrcinp(2, timePrev, timeCurr)

! - Prepare stress/strain to impose
    if (loadEpsiType .lt. 2) then
        lLoadGrad = ASTER_FALSE
        do i = 1, 6
            call fointe('F', loadFunc(i), 1, ['INST'], [timeCurr], valeImpo(i), ier)
            if (loadType(i) .eq. 0) then
                valeImpo(i) = valeImpo(i)/coefAdim
            end if
        end do
        ASSERT(defoComp .eq. 'PETIT')
    else if (loadEpsiType .eq. 2) then
        lLoadGrad = ASTER_TRUE
        do i = 1, 9
            call fointe('F', loadFunc(i), 1, ['INST'], [timeCurr], valeImpo(i), ier)
        end do
    end if

! - Compose rotation for load application
    if (lRota) then
        call tnsvec(6, ndim, matrRotaLoad, valeImpo, 1.d0)
        call utbtab('ZERO', 3, 3, matrRotaLoad, pgl, work, vectRotaLoad)
        call tnsvec(3, ndim, vectRotaLoad, valeImpo, 1.d0)
    end if
    if (loadEpsiType .lt. 2) then
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, valeImpo(4), b_incx)
    end if

! - Initialisation of behaviour datastructure - Special for SIMU_POINT_MAT
    call behaviourInitPoint(relaComp, BEHInteg)

! - Newton algorithm ?
    typeMatrTang = nint(carcriList(TYPE_MATR_T))
    lNewtonAlgo = .not. ((loadEpsiType .ge. 1) .and. (typeMatrTang .eq. 0))

! - Given strain: no Newton
    if (.not. lNewtonAlgo) then
        optionNoNewton = 'RAPH_MECA'
        if (lPrintMatr) then
            optionNoNewton = 'FULL_MECA'
        end if
        if (loadEpsiType .eq. 1) then
            ncmp = 6
            do i = 1, ncmp
                epsiIncr(i) = valeImpo(i)-epsiPrev(i)
            end do
        else if (loadEpsiType .eq. 2) then
            ncmp = 9
            call matinv('S', 3, epsiPrev, fem, jm)
            epsiIncr = reshape(matmul(reshape(valeImpo, (/3, 3/)), reshape(fem, (/3, 3/))), (/9/))
            call lcdetf(3, epsiIncr, jd)
            jp = jm*jd
        end if
        b_n = to_blas_int(nbVari)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jvVim), b_incx, zr(jvVimCopy), b_incy)
        sigmCurr = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, optionNoNewton, typmod, &
                    timePrev, timeCurr, &
                    comporList, carcriList, multComp, &
                    ncmp, epsiPrev, epsiIncr, &
                    6, sigmPrev, &
                    zr(jvVimCopy), &
                    sigmCurr, zr(jvVip), &
                    6*ncmp, dsidep, iret)
        if (defoComp .eq. 'SIMO_MIEHE') then
            b_n = to_blas_int(2*ndim)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/jp, sigmCurr, b_incx)
        end if
        call pmimpr(0, &
                    timeCurr, iterNewt, &
                    loadType, valeImpo, &
                    epsiPrev, sigmPrev, nbVari, zr(jvVim), resi)
        if (iret .ne. 0) then
            liccvg(2) = 1
            goto 500
        end if
        goto 550
    end if

! ==================================================================================================
!
! Newton algorithm
!
! ==================================================================================================
    ASSERT(lNewtonAlgo)

! - Prepare unknowns (YM)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigmPrev, b_incx, ym, b_incy)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    call dscal(b_n, 1.d0/coefAdim, ym, b_incx)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsiPrev, b_incx, ym(7), b_incy)

! - Prediction (Euler)
    iterNewt = 0
    if (typeMatrPred .eq. 1) then
        dy = 0.d0
        epsiIncr = 0.d0
        optionPred = 'RIGI_MECA_TANG'
        b_n = to_blas_int(nbVari)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jvVim), b_incx, zr(jvVipCopy), b_incy)
        sigmCurrPred = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, optionPred, typmod, &
                    timePrev, timeCurr, &
                    comporList, carcriList, multComp, &
                    6, epsiPrev, epsiIncr, &
                    6, sigmPrev, &
                    zr(jvVim), &
                    sigmCurrPred, zr(jvVipCopy), &
                    36, dsidep, &
                    iret)
        if (iret .ne. 0) then
            typeMatrPred = 0
        else
            call pmdrdy(dsidep, coefAdim, coefImpo, valeImpo, ym, &
                        sigmPrev, resi, drdy)
        end if

    else if ((typeMatrPred .eq. 0) .or. ((typeMatrPred .eq. -1) .and. (numeInst .eq. 1))) then
        dy = 0.d0
        epsiIncr = 0.d0
        call pmdrdy(matrElas, coefAdim, coefImpo, valeImpo, ym, &
                    sigmPrev, resi, drdy)
    end if

!        SAUVEGARDE DE R(DY0) POUR TEST DE CONVERGENCE
    b_n = to_blas_int(12)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, resi, b_incx, resiInit, b_incy)
    call pmimpr(0, &
                timeCurr, iterNewt, &
                loadType, valeImpo, &
                epsiPrev, sigmPrev, nbVari, zr(jvVim), resi)

! - Correction (Newton)
300 continue
!
    iterNewt = iterNewt+1
!
    if ((iterNewt .eq. 1) .and. (typeMatrPred .eq. -1) .and. (numeInst .gt. 1)) then
!   prediction='extrapole'
        coefextra = (timeCurr-timePrev)/(timePrev-diinst(sddisc, numeInst-2))
!       dy = dy * (ti - ti-1)/(ti-1 - ti-2)
        b_n = to_blas_int(12)
        b_incx = to_blas_int(1)
        call dscal(b_n, coefextra, dy, b_incx)
    else
!
        b_n = to_blas_int(12)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, resi, b_incx, ddy, b_incy)
!
!      RESOLUTION DE DRDY*DDY = - R(Y)  CARGAU = 'NCSP'
        cargau = 'NCWP'
        call mgauss(cargau, drdy, ddy, 12, 12, &
                    1, r8b, iret)
        if (iret .ne. 0) then
            liccvg(5) = 1
            conver = ASTER_FALSE
            goto 500
        end if
!
!      REACTUALISATION DE DY = DY + DDY
        b_n = to_blas_int(12)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, ddy, b_incx, dy, b_incy)
    end if
!
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, dy(7), b_incx, epsiIncr, b_incy)
!
!           POUR LE CALCUL DE LA MATRICE TANGENTE PAR PERTURBATION
400 continue
!
!           CALCUL DU RESIDU
    liccvg(2) = 0
    b_n = to_blas_int(nbVari)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jvVim), b_incx, zr(jvVimCopy), b_incy)
    sigmCurr = 0.d0
    call nmcomp(BEHInteg, &
                ndim, option, typmod, &
                timePrev, timeCurr, &
                comporList, carcriList, multComp, &
                6, epsiPrev, epsiIncr, &
                6, sigmPrev, &
                zr(jvVimCopy), &
                sigmCurr, zr(jvVip), &
                36, dsidep, &
                iret)
!
    call pmimpr(1, &
                timeCurr, iterNewt, &
                loadType, valeImpo, &
                epsiIncr, sigmCurr, nbVari, zr(jvVip), resi)

    if (iret .ne. 0) then
        conver = ASTER_FALSE
        liccvg(2) = 1
        goto 500
    end if
!
!           CALCUL EVENTUEL DE LA MATRICE TGTE PAR PERTURBATION
    call pmvtgt(option, carcriList, epsiIncr, sigmCurr, zr(jvVip), &
                nbVari, epsilo, varia, matper, dsidep, &
                smatr, sdeps, sigmCurrPred, zr(jvVipCopy), itgt)
    if (itgt .ne. 0) then
        goto 400
    end if
!
    b_n = to_blas_int(12)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, ym, b_incx, y, b_incy)
    b_n = to_blas_int(12)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, dy, b_incx, y, b_incy)
    if (lMatrElas) then
        call pmdrdy(matrElas, coefAdim, coefImpo, valeImpo, y, &
                    sigmCurr, resi, drdy)
    else
        call pmdrdy(dsidep, coefAdim, coefImpo, valeImpo, y, &
                    sigmCurr, resi, drdy)
    end if

! - Check convergence
    call pmconv(resi, resiInit, resiEval, &
                ds_conv, &
                timeCurr, iterNewt, &
                coefAdim, sigmCurr, &
                conver, lIterNewtMaxi)
!
!           VERIFICATION DES EVENT-DRIVEN
500 continue
    call pmsta1(sigmPrev, sigmCurr, epsiIncr, &
                nbVari, nbVariTabl, &
                zr(jvVim), zr(jvVip), &
                tablNbPara, tablType, &
                tablParaName, tablParaType, tablVale, &
                lLoadGrad, zk8(jvVariName), sddisc, &
                liccvg, lIterNewtMaxi, conver, newtLoopAction)
!
!           ON CONTINUE NEWTON
    if (newtLoopAction .eq. 2) goto 300
!
! ======================================================================
!     FIN DES ITERATIONS DE NEWTON
! ======================================================================
!
!        GESTION DE LA DECOUPE DU PAS DE TEMPS
!        EN L'ABSENCE DE CONVERGENCE ON CHERCHE A SUBDIVISER LE PAS
!        DE TEMPS SI L'UTILISATEUR A FAIT LA DEMANDE
    call pmactn(sddisc, ds_conv, iterNewt, numeInst, lIterNewtMaxi, &
                sderro, liccvg, newtLoopAction, action)
!
! ---    ACTION
!          0 ARRET DU CALCUL
!          1 NOUVEAU PAS DE TEMPS
!          2 ON FAIT DES ITERATIONS DE NEWTON EN PLUS
!          3 ON FINIT LE PAS DE TEMPS
    if (action .eq. 1) then
        goto 600
    else if (action .eq. 2) then
        goto 300
    else if (action .eq. 3) then
        goto 550
    else if (action .eq. 0) then
        goto 550
    end if
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
550 continue
!
! - Adaptation of next time step
    call nmfinp(sddisc, numeInst, lastTimeStep)
    if (.not. lastTimeStep) then
        call pmadat(sddisc, numeInst, iterNewt)
    end if
    numeInst = numeInst+1

! - Save values in table
    call pmstab(sigmPrev, sigmCurr, epsiPrev, epsiIncr, &
                nbVari, zr(jvVim), zr(jvVip), &
                timePrev, timeCurr, iterNewt, &
                tablName, tablType, tablNbPara, &
                tablParaName, tablVale, &
                lLoadGrad, valeImpo, lPrintMatr, dsidep, zk8(jvVariName), &
                nbVariTabl)

! - Print
    call pmimpr(2, &
                timeCurr, iterNewt, &
                loadType, valeImpo, &
                epsiIncr, sigmPrev, nbVari, zr(jvVip), resi)
!
600 continue
!
! --- DERNIER INSTANT DE CALCUL ? -> ON SORT DE STAT_NON_LINE
!
    if (lastTimeStep .or. (action .eq. 0)) then
        goto 900
    end if
    goto 200
!==================================
!     FIN BOUCLE SUR LES INSTANTS
!==================================
!
900 continue
!
!     GESTION DES VARIABLES DE COMMANDE
    call vrcinp(0, timePrev, timeCurr)

    call detrsd('FONCTION', '&&CPM_F0')

! - De-Allocate big objects
    deallocate (tablParaType)
    deallocate (tablParaName)
    deallocate (tablVale)
!
    call jedema()
end subroutine
