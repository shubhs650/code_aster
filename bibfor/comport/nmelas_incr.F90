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
! aslint: disable=W1306
!
subroutine nmelas_incr(BEHInteg, &
                       fami, kpg, ksp, typmod, &
                       jvMaterCode, deps, sigm, option, sigp, &
                       vip, dsidep)
!
    use Behaviour_type
    use tenseur_dime_module, only: sph_norm, deviator, kron, voigt, proten, identity
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/verifepsa.h"
#include "asterfort/verifh.h"
#include "asterfort/verifs.h"
#include "asterfort/verift.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: jvMaterCode
    real(kind=8), intent(in) :: deps(:), sigm(:)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(:), vip(1), dsidep(:, :)
! --------------------------------------------------------------------------------------------------
!     REALISE LA LOI DE VON MISES ISOTROPE ET ELASTIQUE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
! --------------------------------------------------------------------------------------------------
! IN  KPG,KSP : NUMERO DU (SOUS)POINT DE GAUSS
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  DEPS    : INCREMENT DE DEFORMATION
!               SI C_PLAN DEPS(3) EST EN FAIT INCONNU (ICI:0)
!                 =>  ATTENTION LA PLACE DE DEPS(3) EST ALORS UTILISEE.
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: elasID = ELAS_ISOT
    character(len=16), parameter :: elasKeyword = 'ELAS'
! --------------------------------------------------------------------------------------------------
    aster_logical :: cplan, resi, rigi
    integer(kind=8)       :: ndimsi, k, l
    real(kind=8)  :: kr(size(deps))
    real(kind=8)  :: em, num, lambdam, deuxmum, troiskm
    real(kind=8)  :: ep, nup, lambdap, deuxmup, troiskp
    real(kind=8)  :: e, nu, lambda, deuxmu, troisk
    real(kind=8)  :: dEpsiTher, dEpsiHydr, dEpsiSech, dEpsiAnel(6), dEpsiVarc(size(deps))
    real(kind=8)  :: dEpsiMeca(size(deps)), sigmp(size(deps))
! --------------------------------------------------------------------------------------------------

    ! Initialisation
    ndimsi = size(deps)
    kr = kron(ndimsi)
    cplan = typmod(1) .eq. 'C_PLAN'
    rigi = option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA'
    resi = option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA' .or. &
           option .eq. 'RIGI_MECA_IMPLEX'

    ! Caracteristiques elastiques t- et t+

    call get_elas_para(fami, jvMaterCode, '-', kpg, ksp, elasID, elasKeyword, e_=em, nu_=num, &
                       BEHInteg=BEHInteg)
    lambdam = em*num/((1-2*num)*(1+num))
    deuxmum = em/(1+num)
    troiskm = em/(1-2*num)

    call get_elas_para(fami, jvMaterCode, '+', kpg, ksp, elasID, elasKeyword, e_=ep, nu_=nup, &
                       BEHInteg=BEHInteg)
    lambdap = ep*nup/((1-2*nup)*(1+nup))
    deuxmup = ep/(1+nup)
    troiskp = ep/(1-2*nup)

    e = merge(ep, em, resi)
    nu = merge(nup, num, resi)
    lambda = merge(lambdap, lambdam, resi)
    deuxmu = merge(deuxmup, deuxmum, resi)
    troisk = merge(troiskp, troiskm, resi)

    ! Contrainte initiale corrigee
    sigmp = deuxmup/deuxmum*deviator(sigm)+troiskp/troiskm*dot_product(kr, sigm)/3.d0*kr

! - Get increment of external state variables
    call verift(fami, kpg, ksp, 'T', jvMaterCode, epsth_=dEpsiTher)
    call verifh(fami, kpg, ksp, 'T', jvMaterCode, dEpsiHydr)
    call verifs(fami, kpg, ksp, 'T', jvMaterCode, dEpsiSech)
    call verifepsa(fami, kpg, ksp, 'T', dEpsiAnel)
    dEpsiVarc = (dEpsiTher+dEpsiHydr+dEpsiSech)*kr+dEpsiAnel(1:ndimsi)*voigt(ndimsi)

    ! Calcul de la contrainte

    if (resi) then
        dEpsiMeca = deps-dEpsiVarc
        if (cplan) dEpsiMeca(3) = -(lambda*(dEpsiMeca(1)+dEpsiMeca(2))+sigmp(3))/(lambda+deuxmu)

        sigp = sigmp+lambda*dot_product(kr, dEpsiMeca)*kr+deuxmu*dEpsiMeca
        vip(1) = 0.d0

    else
        ! Pour la prediction (loi old)
        sigp = sigmp

    end if

    ! Matrice tangente (ancienne formulation de la prediction)

    if (rigi) then
        dsidep = lambda*proten(kr, kr)+deuxmu*identity(ndimsi)

        if (cplan) then
            do k = 1, ndimsi
                if (k .ne. 3) then
                    do l = 1, ndimsi
                        if (l .ne. 3) then
                            dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep(3, l)
                        end if
                    end do
                end if
            end do
        end if
    end if

end subroutine
