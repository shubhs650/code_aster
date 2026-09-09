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
subroutine vdxnlr(plateCara, plateOrie, &
                  BEHInteg, &
                  option, nomte, nodeCoor, &
                  matrTang, codret)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btkb.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/epseff.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrc2.h"
#include "asterfort/matrkb.h"
#include "asterfort/moytpg.h"
#include "asterfort/nmcomp.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/trndgl.h"
#include "asterfort/trnflg.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "asterfort/vexpan.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=16), intent(in) :: option, nomte
    real(kind=8), intent(in) :: nodeCoor(3, 9)
    real(kind=8), intent(out) :: matrTang(51, 51)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 3
    character(len=16), parameter :: multComp = " "
    integer(kind=8) :: nb1, nb2, nddle, npgsr, npgsn, itab(8)
    integer(kind=8) :: cod, ksp
    real(kind=8) :: vectBase(9, 3, 3)
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42), wmatcb(5, 42), ktildi(42, 42)
    real(kind=8) :: ktild(42, 42)
    real(kind=8) :: ctor, eptot, kappa
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/'E ', 'NU'/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    real(kind=8) :: rotfcm(9), rotfcp(9)
    real(kind=8) :: deplm(42), deplp(42)
    real(kind=8) :: epsi(5), depsi(5), eps2d(4), deps2d(4)
    real(kind=8) :: dtild(5, 5), sgmtd(5), effint(42), vecl(48), vecll(51)
    real(kind=8) :: sign(4), sigma(4), dsidep(6, 6)
    real(kind=8) :: matrElas(5, 5), tempMoye
    integer(kind=8) :: i, ib, jvCarcri, icontm, icontp, iLayer
    integer(kind=8) :: ideplm, ideplp, iinstm, iinstp, inte, intsn
    integer(kind=8) :: intsr, iret, ivarim, ivarip, ivarix, ivectu, j
    integer(kind=8) :: jcrf, k1, k2, kpgs, kwgt, lgpg
    integer(kind=8) :: lzi, lzr, nbLayer, nbvari, nddlet, ndimv
    real(kind=8) :: coef, crf, gxz, gyz, hLayer
    real(kind=8) :: x(1), zic, zmin
    real(kind=8) :: ksi3s2
    aster_logical :: lVect, lMatr, lVari, lSigm
    blas_int :: b_incx, b_incy, b_n
    real(kind=8) :: cisail
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8), parameter :: ndimLdc = 2
    character(len=8), parameter :: typmod(2) = (/"C_PLAN  ", "        "/)
    character(len=16), pointer :: compor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0

! - Get shell parameters
    nbLayer = plateCara%nbLayer
    ASSERT(nbLayer .gt. 0)
    eptot = plateCara%thick
    kappa = plateCara%shearCoef
    ctor = plateCara%coefRigiDRZ
    zmin = -eptot/2.d0
    hLayer = eptot/nbLayer

! - Geometry
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMP', 'L', ivarix)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=itab)
    if (itab(6) .le. 1) then
        lgpg = itab(7)
    else
        lgpg = itab(6)*itab(7)
    end if

! - Set main parameters for behaviour (on cell)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Properties of behaviour
    read (compor(NVAR), '(I16)') nbvari

! - Get elastic properties
    if (BEHInteg%materPara%elasID .ne. ELAS_ISOT .and. &
        BEHInteg%materPara%elasID .ne. ELAS_ORTH) then
        call utmess('F', 'PLATE1_12', sk=BEHInteg%materPara%elasKeyword)
    end if

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    nddle = 5*nb1+2
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get output fields
    if (option .eq. 'RAPH_MECA') then
        call jevech('PCACO3D', 'L', jcrf)
        crf = zr(jcrf)
    else
        call jevech('PCACO3D', 'E', jcrf)
    end if
    ivarip = ivarim
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    ndimv = lgpg*npgsn
    b_n = to_blas_int(ndimv)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)

! - Change coordinates of displacements/rotations
    call trndgl(nb2, plateOrie%vectNorm, plateOrie%vectTang, zr(ideplm), &
                deplm, rotfcm)
    call trndgl(nb2, plateOrie%vectNorm, plateOrie%vectTang, zr(ideplp), &
                deplp, rotfcp)
!
    ktild = 0.d0
    effint = 0.d0
    kwgt = 0
    kpgs = 0
    do iLayer = 1, nbLayer
        do inte = 1, npge
            if (inte .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
                coef = 1.d0/3.d0
            else if (inte .eq. 2) then
                zic = zmin+hLayer/2.d0+(iLayer-1)*hLayer
                coef = 4.d0/3.d0
            else
                zic = zmin+hLayer+(iLayer-1)*hLayer
                coef = 1.d0/3.d0
            end if
            ksi3s2 = zic/hLayer

            do intsr = 1, npgsr
                call mahsms(plateOrie, &
                            0, nb1, &
                            nodeCoor, ksi3s2, intsr, &
                            zr(lzr), hLayer, &
                            vectBaseKpg, vectTangKpg, &
                            hsfm, hss)
                call hsj1ms(hLayer, vectTangKpg, vectBaseKpg, hsfm, hss, &
                            hsj1m, hsj1s)
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            hLayer, plateOrie%vectTang, hsj1m, hsj1s, btdm, &
                            btds)
            end do

            do intsn = 1, npgsn
!
!     CALCUL DE BTDFN : F=FLEXION , N=NORMAL
!     ET DEFINITION DE WGT=PRODUIT DES POIDS ASSOCIES AUX PTS DE GAUSS
!                          (NORMAL) ET DU DETERMINANT DU JACOBIEN
!
                call mahsf(plateOrie, &
                           1, nb1, &
                           nodeCoor, ksi3s2, intsn, &
                           zr(lzr), hLayer, &
                           vectBaseKpg, vectTangKpg, &
                           hsf)
                call hsj1f(intsn, zr(lzr), hLayer, vectTangKpg, vectBaseKpg, &
                           hsf, kwgt, hsj1fx, wgt)
!
!     PRODUIT DU POIDS DES PTS DE GAUSS DANS L'EPAISSEUR ET DE WGT
!
                wgt = coef*wgt
!
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), hLayer, plateOrie%vectTang, hsj1fx, btdf)
!
!     CALCUL DE BTDMN, BTDSN
!     ET
!     FORMATION DE BTILD
!
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
!
!     CALCULS DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES :
!     EPSXX, EPSYY, EPSXY, EPSXZ, EPSYZ (CE SONT LES COMPOSANTES TILDE)
                kpgs = kpgs+1
                call epseff('DEFORM', nb1, deplm, btild, x, &
                            epsi, wgt, x)
                eps2d(1) = epsi(1)
                eps2d(2) = epsi(2)
                eps2d(3) = 0.d0
                eps2d(4) = epsi(3)/rac2
!
                call epseff('DEFORM', nb1, deplp, btild, x, &
                            depsi, wgt, x)
                deps2d(1) = depsi(1)
                deps2d(2) = depsi(2)
                deps2d(3) = 0.d0
                deps2d(4) = depsi(3)/rac2
!
                gxz = epsi(4)+depsi(4)
                gyz = epsi(5)+depsi(5)
!
                k1 = 6*((intsn-1)*npge*nbLayer+(iLayer-1)*npge+inte-1)
                k2 = lgpg*(intsn-1)+(npge*(iLayer-1)+inte-1)*nbvari
                do i = 1, 3
                    sign(i) = zr(icontm-1+k1+i)
                end do
                sign(4) = zr(icontm-1+k1+4)*rac2

                cisail = 0.d0

! ------------- Index of "sub"-point
                ksp = (iLayer-1)*npge+inte

! ------------- Set main parameters for behaviour (on point)
                call behaviourSetParaPoin(intsn, ksp, BEHInteg)

! ------------- Integrator
                if (BEHInteg%materPara%elasID .eq. ELAS_ISOT) then
                    sigma = 0.d0
                    call nmcomp(BEHInteg, &
                                ndimLdc, option, typmod, &
                                zr(iinstm), zr(iinstp), &
                                compor, zr(jvCarcri), multComp, &
                                4, eps2d, deps2d, &
                                4, sign, &
                                zr(ivarim+k2), &
                                sigma, zr(ivarip+k2), &
                                36, dsidep, cod)

                    call rcvalb(BEHInteg%materPara%schemePara%fami, &
                                BEHInteg%materPara%schemePara%kpg, &
                                BEHInteg%materPara%schemePara%ksp, &
                                '+', BEHInteg%materPara%jvMaterCode, &
                                ' ', BEHInteg%materPara%elasKeyword, &
                                0, ' ', [0.d0], &
                                nbProp, propName, propVale, &
                                propCode, 1)
                    cisail = propVale(1)/(1.d0+propVale(2))

!           COD=1 : ECHEC INTEGRATION LOI DE COMPORTEMENT
!           COD=3 : C_PLAN DEBORST SIGZZ NON NUL
                    if (cod .ne. 0) then
                        if (codret .ne. 1) then
                            codret = cod
                        end if
                        if (cod .eq. 1) goto 999
                    end if
                else if (BEHInteg%materPara%elasID .eq. ELAS_ORTH) then
                    call moytpg('RIGI', intsn, 3, '+', tempMoye, iret)
                    call matrc2(plateOrie, vectBaseKpg, tempMoye, kappa, matrElas)
                end if
!
!    CALCULS DE LA MATRICE TANGENTE : BOUCLE SUR L'EPAISSEUR
                if (lMatr) then
                    if (BEHInteg%materPara%elasID .eq. ELAS_ISOT) then
                        dtild(1, 1) = dsidep(1, 1)
                        dtild(1, 2) = dsidep(1, 2)
                        dtild(1, 3) = dsidep(1, 4)/rac2
                        dtild(1, 4) = 0.d0
                        dtild(1, 5) = 0.d0
                        dtild(2, 1) = dsidep(2, 1)
                        dtild(2, 2) = dsidep(2, 2)
                        dtild(2, 3) = dsidep(2, 4)/rac2
                        dtild(2, 4) = 0.d0
                        dtild(2, 5) = 0.d0
                        dtild(3, 1) = dsidep(4, 1)/rac2
                        dtild(3, 2) = dsidep(4, 2)/rac2
                        dtild(3, 3) = dsidep(4, 4)/2.d0
                        dtild(3, 4) = 0.d0
                        dtild(3, 5) = 0.d0
                        dtild(4, 1) = 0.d0
                        dtild(4, 2) = 0.d0
                        dtild(4, 3) = 0.d0
                        dtild(4, 4) = cisail*kappa/2.d0
                        dtild(4, 5) = 0.d0
                        dtild(5, 1) = 0.d0
                        dtild(5, 2) = 0.d0
                        dtild(5, 3) = 0.d0
                        dtild(5, 4) = 0.d0
                        dtild(5, 5) = cisail*kappa/2.d0

                    else if (BEHInteg%materPara%elasID .eq. ELAS_ORTH) then
                        dtild(1, 1) = matrElas(1, 1)
                        dtild(1, 2) = matrElas(1, 2)
                        dtild(1, 3) = matrElas(1, 3)
                        dtild(1, 4) = 0.d0
                        dtild(1, 5) = 0.d0
                        dtild(2, 1) = matrElas(2, 1)
                        dtild(2, 2) = matrElas(2, 2)
                        dtild(2, 3) = matrElas(2, 3)
                        dtild(2, 4) = 0.d0
                        dtild(2, 5) = 0.d0
                        dtild(3, 1) = matrElas(3, 1)
                        dtild(3, 2) = matrElas(3, 2)
                        dtild(3, 3) = matrElas(3, 3)
                        dtild(3, 4) = 0.d0
                        dtild(3, 5) = 0.d0
                        dtild(4, 1) = 0.d0
                        dtild(4, 2) = 0.d0
                        dtild(4, 3) = 0.d0
                        dtild(4, 4) = matrElas(4, 4)
                        dtild(4, 5) = matrElas(4, 5)
                        dtild(5, 1) = 0.d0
                        dtild(5, 2) = 0.d0
                        dtild(5, 3) = 0.d0
                        dtild(5, 4) = matrElas(5, 4)
                        dtild(5, 5) = matrElas(5, 5)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
!
                    b_n = to_blas_int(25)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, wgt, dtild, b_incx)
!
                    call btkb(5, 42, nddle, dtild, btild, &
                              wmatcb, ktildi)
!
                    do i = 1, nddle
                        do j = 1, nddle
                            ktild(i, j) = ktild(i, j)+ktildi(i, j)
                        end do
                    end do
                end if
!
                if (lSigm) then
                    ASSERT(lVect)
                    if (BEHInteg%materPara%elasID .eq. ELAS_ISOT) then
                        do i = 1, 3
                            zr(icontp-1+k1+i) = sigma(i)
                        end do
                        zr(icontp-1+k1+4) = sigma(4)/rac2
                        zr(icontp-1+k1+5) = cisail*kappa*gxz/2.d0
                        zr(icontp-1+k1+6) = cisail*kappa*gyz/2.d0
!
!    CALCULS DES EFFORTS INTERIEURS
                        sgmtd(1) = zr(icontp-1+k1+1)
                        sgmtd(2) = zr(icontp-1+k1+2)
                        sgmtd(3) = zr(icontp-1+k1+4)
                        sgmtd(4) = cisail*kappa*gxz/2.d0
                        sgmtd(5) = cisail*kappa*gyz/2.d0
!
                    else if (BEHInteg%materPara%elasID .eq. ELAS_ORTH) then
                        zr(icontp-1+k1+1) = (epsi(1)+depsi(1))*matrElas(1, 1)+ &
                                            (epsi(2)+depsi(2))*matrElas(1, 2)+ &
                                            (epsi(3)+depsi(3))*matrElas(1, 3)
                        zr(icontp-1+k1+2) = (epsi(1)+depsi(1))*matrElas(2, 1)+ &
                                            (epsi(2)+depsi(2))*matrElas(2, 2)+ &
                                            (epsi(3)+depsi(3))*matrElas(2, 3)
                        zr(icontp-1+k1+3) = 0.d0
                        zr(icontp-1+k1+4) = (epsi(1)+depsi(1))*matrElas(3, 1)+ &
                                            (epsi(2)+depsi(2))*matrElas(3, 2)+ &
                                            (epsi(3)+depsi(3))*matrElas(3, 3)
                        zr(icontp-1+k1+5) = matrElas(4, 4)*gxz+matrElas(4, 5)*gyz
                        zr(icontp-1+k1+6) = matrElas(5, 4)*gxz+matrElas(5, 5)*gyz
!
!    CALCULS DES EFFORTS INTERIEURS
                        sgmtd(1) = zr(icontp-1+k1+1)
                        sgmtd(2) = zr(icontp-1+k1+2)
                        sgmtd(3) = zr(icontp-1+k1+4)
                        sgmtd(4) = dtild(4, 4)*gxz
                        sgmtd(5) = dtild(5, 5)*gyz
                    end if
                    call epseff('EFFORI', nb1, x, btild, sgmtd, x, wgt, effint)
                end if
            end do
        end do
    end do
!
    if (lMatr) then
!
!     EXPANSION DE LA MATRICE : AJOUTER DE LA ROTATION FICTIVE
!
        nddlet = 6*nb1+3
        call matrkb(nb1, 42, 51, nddlet, ktild, &
                    ctor, matrTang, crf)
        zr(jcrf) = crf
!
!     AJOUTER DES 3 TRANSLATIONS FICTIVES ASSOCIEES AU NOEUD INTERNE
!     LES 3 TERMES DE RAIDEUR (FICTIVE) ASSOCIEES ONT POUR VALEUR CELLE
!     DES ROTATION FICTIVE
!
    end if
!
    if (lVect) then
!
        call vexpan(nb1, effint, vecl)
!
        do i = 1, 6*nb1
            vecll(i) = vecl(i)
        end do
        vecll(6*nb1+1) = effint(5*nb1+1)
        vecll(6*nb1+2) = effint(5*nb1+2)
        vecll(6*nb1+3) = 0.d0
!
!     CONTRIBUTION DES DDL DE LA ROTATION FICTIVE DANS EFFINT
!
        do i = 1, nb1
            vecll(6*i) = crf*(rotfcm(i)+rotfcp(i))
        end do
        i = nb2
        vecll(6*nb1+3) = crf*(rotfcm(nb2)+rotfcp(nb2))

! ----- Fuse tangents and normal in same object
        do ib = 1, nb2
            vectBase(ib, 1:2, 1:3) = plateOrie%vectTang(ib, 1:2, 1:3)
            vectBase(ib, 3, 1:3) = plateOrie%vectNorm(ib, 1:3)
        end do

        call trnflg(nb2, vectBase, vecll, zr(ivectu))
    end if
!
999 continue
end subroutine
