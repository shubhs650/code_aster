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
! aslint: disable=W0413
!
subroutine enelpg(materPara, time, &
                  relaName, defoComp, &
                  f, sigmEner, &
                  nbVari, vari, &
                  enerElas)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/d1mamc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nrsmt1.h"
#include "asterfort/nrsmtb.h"
#include "asterfort/nrsmtt.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/zerop3.h"
#include "blas/dcopy.h"
!
    type(Material_Para), intent(inout) :: materPara
    real(kind=8), intent(in) :: time
    character(len=16), intent(in) :: relaName, defoComp
    real(kind=8), intent(in) :: f(3, 3), sigmEner(6)
    integer(kind=8), intent(in) :: nbVari
    real(kind=8), intent(in) :: vari(nbVari)
    real(kind=8), intent(out) :: enerElas
!
! --------------------------------------------------------------------------------------------------
!
! CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!
! DETERMINEE PAR L'EXPRESSION SUIVANTE :
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
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbsig, nsol, i, isig, jsig
    real(kind=8) :: c1, c2, trt
    real(kind=8) :: sol(3), jzero, uzero, mzero, epsi(6)
    real(kind=8) :: mjac, ujac, wbe, be(6), e, nu
    real(kind=8) :: mu, troisk, jac, tau(6), trtau, eqtau, dvtau(6), tlog(6)
    real(kind=8) :: trbe, epsthe, d1(36)
    blas_int :: b_incx, b_incy, b_n
    real(kind=8), parameter :: pdtsca(6) = (/1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 2.d0/)
    real(kind=8), parameter :: kr(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!
    nbsig = nbsigm()
    ASSERT(nbsig .le. 6)
    enerElas = 0.d0

! --- CAS EN GRANDES DEFORMATIONS SIMO_MIEHE
    if ((defoComp .eq. 'SIMO_MIEHE') .and. &
        ((relaName(1:9) .eq. 'VMIS_ISOT') .or. (relaName .eq. 'ELAS'))) then

! ----- Get elastic parameters
        if (materPara%elasID .ne. ELAS_ISOT) then
            call utmess("F", "ENERGY1_2", sk=materPara%elasKeyword)
        end if
        call get_elas_para(materPara%schemePara%fami, &
                           materPara%jvMaterCode, '+', &
                           materPara%schemePara%kpg, &
                           materPara%schemePara%ksp, &
                           materPara%elasID, materPara%elasKeyword, &
                           e_=e, nu_=nu)
        mu = e/(2.d0*(1.d0+nu))
        troisk = e/(1.d0-2.d0*nu)

! ----- Compute jacobian of transformation
        jac = f(1, 1)*(f(2, 2)*f(3, 3)-f(2, 3)*f(3, 2))- &
              f(2, 1)*(f(1, 2)*f(3, 3)-f(1, 3)*f(3, 2))+ &
              f(3, 1)*(f(1, 2)*f(2, 3)-f(1, 3)*f(2, 2))

! ----- CALCUL DE TAU TEL QUE TAU=JAC*SIGMA
        tau = 0.d0
        do i = 1, nbsig
            tau(i) = jac*sigmEner(i)
        end do

! ----- CALCUL DE LA TRACE DE TAU- TAU EQUIVALENT ET TAU DEVIATORIQUE
        trtau = tau(1)+tau(2)+tau(3)
        eqtau = 0.d0
        do i = 1, 6
            dvtau(i) = tau(i)-kr(i)*trtau/3.d0
            eqtau = eqtau+pdtsca(i)*(dvtau(i)**2.d0)
        end do
        eqtau = sqrt(1.5d0*eqtau)

! ----- CALCUL DE LA TRACE DES DEFORMATIONS ELASTIQUES BE
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vari(3), b_incx, be, b_incy)
        trbe = be(1)+be(2)+be(3)
        trbe = jac**(-2.d0/3.d0)*(3.d0-2.d0*trbe)

! ----- DEFORMATION THERMIQUE AU POINT D'INTEGRATION COURANT :
        call verift(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', materPara%jvMaterCode, &
                    epsth_=epsthe)

! ----- ATTENTION, EN PRESENCE DE THERMIQUE, CA MET LE BAZAR...
        sol = 0.d0
        jzero = 0.d0
        uzero = 0.d0
        mzero = 0.d0
        nsol = 0
        mjac = 0.d0
        if (epsthe .ne. 0) then
            call zerop3(-3.d0*epsthe, -1.d0, -3.d0*epsthe, sol, nsol)
            jzero = sol(1)
            call nrsmt1(troisk/3.d0, jzero, uzero)
            call nrsmtt(troisk, jzero, epsthe, mzero)
            call nrsmtt(troisk, jac, epsthe, mjac)
        end if

! ----- CALCUL DES TERMES DE L'ENERGIE
        ujac = 0.d0
        call nrsmt1(troisk/3.d0, jac, ujac)
        wbe = 0.d0
        call nrsmtb(mu, trbe, wbe)
        enerElas = ujac+wbe+mjac-uzero-mzero

    else if ((defoComp(1:8) .eq. 'GDEF_LOG')) then
! ----- Get elastic parameters
        if (materPara%elasID .ne. ELAS_ISOT) then
            call utmess("F", "ENERGY1_2", sk=materPara%elasKeyword)
        end if
        call get_elas_para(materPara%schemePara%fami, &
                           materPara%jvMaterCode, '+', &
                           materPara%schemePara%kpg, &
                           materPara%schemePara%ksp, &
                           materPara%elasID, materPara%elasKeyword, &
                           e_=e, nu_=nu)
        mu = e/(2.d0*(1.d0+nu))
        troisk = e/(1.d0-2.d0*nu)
!
        c1 = (1.d0+nu)/e
        c2 = nu/e
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vari(nbVari-5), b_incx, tlog, b_incy)

        if (lteatt('DIM_TOPO_MAILLE', '3')) then
            trt = tlog(1)+tlog(2)+tlog(3)
            enerElas = 0.5d0*(tlog(1)*(c1*tlog(1)-c2*trt)+ &
                              tlog(2)*(c1*tlog(2)-c2*trt)+ &
                              tlog(3)*(c1*tlog(3)-c2*trt)+ &
                              (tlog(4)*c1*tlog(4)+ &
                               tlog(5)*c1*tlog(5)+ &
                               tlog(6)*c1*tlog(6)))

        else if (lteatt('C_PLAN', 'OUI')) then
            trt = tlog(1)+tlog(2)
            enerElas = 0.5d0*(tlog(1)*(c1*tlog(1)-c2*trt)+ &
                              tlog(2)*(c1*tlog(2)-c2*trt)+ &
                              2.d0*tlog(4)*c1*tlog(4))
        else
            trt = tlog(1)+tlog(2)+tlog(3)
            enerElas = 0.5d0*(tlog(1)*(c1*tlog(1)-c2*trt)+ &
                              tlog(2)*(c1*tlog(2)-c2*trt)+ &
                              tlog(3)*(c1*tlog(3)-c2*trt)+ &
                              2.d0*tlog(4)*c1*tlog(4))
        end if

    else if (defoComp(1:5) .eq. 'PETIT') then

! ----- CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE
        call d1mamc(materPara, '+', time, nbSig, d1)

! ----- DENSITE D'ENERGIE POTENTIELLE ELASTIQUE AU POINT D'INTEGRATION COURANT
        epsi = 0.d0
        do isig = 1, nbsig
            do jsig = 1, nbsig
                epsi(isig) = epsi(isig)+d1(nbsig*(isig-1)+jsig)*sigmEner(jsig)
            end do
            enerElas = enerElas+0.5d0*sigmEner(isig)*epsi(isig)
        end do

    else
        call utmess('F', "ENERGY1_1", sk=defoComp)
    end if
!
end subroutine
