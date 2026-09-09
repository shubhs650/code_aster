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
subroutine lcmate(materPara, &
                  carcri, relaComp, typmod1, &
                  nmat, tempd, tempf, tref, rungeKutta, &
                  typma, hsr, materd, materf, matcst, &
                  nbcomm, cpmono, pgl, itmax, &
                  toler, ndt, ndi, nr, &
                  nvi, vind, nfs, nsg, toutms, &
                  nhsr, numhsr, sigd, multComp_)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/cvmmat.h"
#include "asterfort/haymat.h"
#include "asterfort/hbrmat.h"
#include "asterfort/irrmat.h"
#include "asterfort/lcmatt.h"
#include "asterfort/lcmmap.h"
#include "asterfort/lcmmat.h"
#include "asterfort/lglmat.h"
#include "asterfort/lkimat.h"
#include "asterfort/matect.h"
#include "asterfort/rslmat.h"
#include "asterfort/rsvmat.h"
#include "asterfort/srimat.h"
#include "asterfort/vecmat.h"
!
    type(Material_Para), intent(in) :: materPara
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: relaComp
    character(len=8), intent(in) :: typmod1
    integer(kind=8), intent(in) :: nvi
    character(len=16), optional, intent(in) :: multComp_
!
! --------------------------------------------------------------------------------------------------
!
!       RECUPERATION DU MATERIAU A TEMPF ET TEMPD
!
! --------------------------------------------------------------------------------------------------
!
!           NMAT   :  DIMENSION 1 DE MATER
!           TEMPD  :  TEMPERATURE A T
!           TEMPF  :  TEMPERATURE A T + DT
!           IMPEXP : 0 IMPLICITE, 1 EXPLICITE
!           SIGD   :  ETAT DE CONTRAINTES A T
!       OUT MATERD :  COEFFICIENTS MATERIAU A T    (TEMPD )
!           MATERF :  COEFFICIENTS MATERIAU A T+DT (TEMPF )
!                     MATER(*,I) = CARACTERISTIQUES MATERIAU
!                                    I = 1  CARACTERISTIQUES ELASTIQUES
!                                    I = 2  CARACTERISTIQUES PLASTIQUES
!           MATCST :  'OUI' SI  MATERIAU A T = MATERIAU A T+DT
!                     'NON' SINON OU 'NAP' SI NAPPE DANS 'VECMAT.F'
!           NBCOMM : POSITION DES COEF POUR CHAQUE LOI DE CHAQUE SYSTEME
!           CPMONO : NOMS DES LOIS POUR CHAQUE FAMILLE DE SYSTEME
!           PGL    : MATRICE DE PASSAGE
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!           NR     :  NB DE COMPOSANTES SYSTEME NL
!           NVI    :  NB DE VARIABLES INTERNES
!           TOUTMS :  TOUS LES TENSEURS MS
!           HSR    : MATRICE D'INTERACTION POUR L'ECROUISSAGE ISOTROPE
!                    UTILISEE SEULEMENT POUR LE MONOCRISTAL IMPLICITE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ::  nmat, ndt, ndi, nr, i, itmax, rungeKutta
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), tempd, tempf, tref
    real(kind=8) :: vind(*), pgl(3, 3), toler, sigd(6)
    character(len=8) :: typma
    character(len=3) :: matcst
    integer(kind=8) :: numhsr(*), nbcomm(*), nfs, nsg, nhsr
    real(kind=8) :: hsr(*), toutms(*)
    character(len=24) :: cpmono(*)
    character(len=8) :: fami
    integer(kind=8) :: jvMaterCode, kpg, ksp
    character(len=16) :: multComp
!
! --------------------------------------------------------------------------------------------------
!
    do i = 1, nmat
        materd(i, 1) = 0.d0
        materd(i, 2) = 0.d0
        materf(i, 1) = 0.d0
        materf(i, 2) = 0.d0
    end do
! - For number of phases when is not a crystal behaviour (issue30310)
    nbcomm(1) = 1

! - Access to material parameters
    jvMaterCode = materPara%jvMaterCode
    fami = materPara%schemePara%fami
    kpg = materPara%schemePara%kpg
    ksp = materPara%schemePara%ksp

!
    multComp = ' '
    if (present(multComp_)) then
        multComp = multComp_
    end if
    if (relaComp .eq. 'ROUSS_PR') then
        call rslmat(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, materd, materf, matcst, ndt, &
                    ndi, nr, nvi, vind)
!
    else if (relaComp .eq. 'ROUSS_VISC') then
        call rsvmat(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, materd, materf, matcst, ndt, &
                    ndi, nr, nvi, vind)
!
    else if (relaComp .eq. 'VISCOCHAB') then
        call cvmmat(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, materd, materf, matcst, typma, &
                    ndt, ndi, nr, carcri, vind, &
                    nvi, sigd)
!
    else if (relaComp .eq. 'VENDOCHAB' .or. relaComp .eq. 'VISC_ENDO_LEMA') then
        call vecmat(fami, kpg, ksp, typmod1, relaComp, &
                    jvMaterCode, nmat, materd, materf, matcst, &
                    typma, ndt, ndi, nr, nvi)
!
    else if (relaComp(1:6) .eq. 'LAIGLE') then
        call lglmat(typmod1, jvMaterCode, nmat, tempd, materd, &
                    materf, matcst, ndt, ndi, nr, &
                    nvi)
!
    elseif ((relaComp .eq. 'HOEK_BROWN') .or. (relaComp .eq. 'HOEK_BROWN_EFF')) then
        call hbrmat(typmod1, jvMaterCode, nmat, tempd, materd, &
                    materf, matcst, ndt, ndi, nr, &
                    nvi)
!
    else if (relaComp .eq. 'MONOCRISTAL') then
        ASSERT(multComp .ne. ' ')
        call lcmmat(materPara, &
                    multComp, typmod1, &
                    nmat, pgl, materd, &
                    materf, matcst, nbcomm, cpmono, ndt, &
                    ndi, nr, nvi, hsr, nfs, &
                    nsg, toutms, vind, rungeKutta)
        typma = 'COHERENT'
        if (typmod1 .ne. '3D') then
            sigd(5) = 0.d0
            sigd(6) = 0.d0
        end if
!
    else if (relaComp .eq. 'POLYCRISTAL') then
        ASSERT(multComp .ne. ' ')
        call lcmmap(materPara, &
                    multComp, typmod1, &
                    nmat, pgl, materd, &
                    materf, matcst, nbcomm, cpmono, ndt, &
                    ndi, nr, nvi, nfs, nsg, &
                    nhsr, numhsr, hsr)
        typma = 'COHERENT'

    else if (relaComp .eq. 'IRRAD3M') then
        call irrmat(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, itmax, toler, materd, materf, &
                    matcst, ndt, ndi, nr, nvi)

    else if (relaComp .eq. 'LETK') then
        call lkimat(typmod1, jvMaterCode, nmat, materd, materf, &
                    matcst, ndt, ndi, nvi, nr)
        typma = 'COHERENT'

    else if (relaComp .eq. 'LKR') then
        call srimat(typmod1, jvMaterCode, nmat, tempd, tempf, tref, materd, materf, &
                    matcst, ndt, ndi, nvi, nr)
        typma = 'COHERENT'

    else if (relaComp .eq. 'HAYHURST') then
        call haymat(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, '-', materd(1, 1), materd(1, 2), nvi, &
                    nr)
        call haymat(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, '+', materf(1, 1), materf(1, 2), nvi, &
                    nr)
        call matect(materd, materf, nmat, matcst)
        typma = 'COHERENT'

    else
        call lcmatt(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, '-', relaComp, materd(1, 1), materd(1, 2), &
                    typma, ndt, ndi, nr, nvi)
        call lcmatt(fami, kpg, ksp, typmod1, jvMaterCode, &
                    nmat, '+', relaComp, materf(1, 1), materf(1, 2), &
                    typma, ndt, ndi, nr, nvi)
        call matect(materd, materf, nmat, matcst)

    end if
!
!     - DANS LCPLNL ON DIMENSIONNE DES TABLES AVEC (NDT+NVI) QUI SONT
!       ENSUITE UTILISEES PAR NEWTON
!     - LA DIMENSION DU SYSTEME DIFFERENTIEL EST NR
!     ==> IL FAUT DONC NDT+NVI >= NR
    ASSERT((ndt+nvi) .ge. nr)
!
end subroutine
