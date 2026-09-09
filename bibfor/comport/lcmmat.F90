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
subroutine lcmmat(materPara, &
                  multComp, typmod1, &
                  nmat, pgl, materd, &
                  materf, matcst, nbcomm, cpmono, ndt, &
                  ndi, nr, nvi, hsr, nfs, &
                  nsg, toutms, vind, impexp)
!
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/calcmm.h"
#include "asterfort/d1ma3d.h"
#include "asterfort/dmat3d.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lcmaec.h"
#include "asterfort/lcmaei.h"
#include "asterfort/lcmafl.h"
#include "asterfort/lcmmjv.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/matrot.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: multComp
    character(len=8), intent(in)  :: typmod1
!
! --------------------------------------------------------------------------------------------------
!
!       MONOCRISTAL : RECUPERATION DU MATERIAU A T(TEMPD) ET T+DT(TEMPF)
!                    NB DE CMP DIRECTES/CISAILLEMENT , NB VAR. INTERNES
!                    MATER(*,1) = E , NU , ALPHA
!                    MATER(*,2) = COEFFICIENT DE CHAQUE COMPORTEMENT
!                    VARIABLES INTERNES :
!                     EPSVP(6)+ALPHA,GAMMA,P PAR SYSTEME DE GLISSEMENT
!
! --------------------------------------------------------------------------------------------------
!
!       IN  FAMI   : FAMILLE DE POINT DE GAUSS
!            KPG   : NUMERO DU POINT DE GAUSS
!            KSP   : NUMERO DU SOUS-POINT DE GAUSS
! In  rela_comp: RELATION for comportment
!           IMAT   :  ADRESSE DU MATERIAU CODE
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION  MAXIMUM DE MATER
!          ANGMAS  :  TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
!       OUT MATERD :  COEFFICIENTS MATERIAU A T
!           PGL    : MATRICE DE PASSAGE GLOBAL LOCAL
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!                     MATER(*,1) = CARACTERISTIQUES   ELASTIQUES
!                     MATER(*,2) = CARACTERISTIQUES   PLASTIQUES
!           MATCST :  'OUI' SI  MATERIAU A T = MATERIAU A T+DT
!                     'NON' SINON
!           NBCOMM : POSITION DES COEF POUR CHAQUE LOI DE CHAQUE SYSTEME
!           CPMONO : NOMS DES LOIS POUR CHAQUE FAMILLE DE SYSTEME
!
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!           NR     :  NB DE COMPOSANTES SYSTEME NL
!           NVI    :  NB DE VARIABLES INTERNES
!           HSR    : MATRICE D'INTERACTION POUR L'ECROUISSAGE ISOTROPE
!                    UTILISEE SEULEMENT POUR LE MONOCRISTAL IMPLICITE
!           TOUTMS : TOUS LES TENSEURS D'ORIENTATION POUR TOUS LES
!                    SYSTEMES DE GLISSEMENT
!           IMPEXP : 0 IMPLICITE, 1 EXPLICITE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nfs, nsg, nmat
    integer(kind=8) :: ndt, ndi, nr, nvi, nbcomm(nmat, 3)
    integer(kind=8) :: irota, impexp, nbval, nvini
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: hooke(6, 6), hookef(6, 6)
    real(kind=8) :: ekooh(6, 6), fekooh(6, 6)
    real(kind=8) :: epsi, pgl(3, 3), q(3, 3)
    real(kind=8) :: valres(nmat), ms(6), ng(3), lg(3), vind(*)
    real(kind=8) :: hsr(nsg, nsg), toutms(nfs, nsg, 6)
    character(len=8) :: propName(3)
    integer(kind=8) :: cerr(3)
    character(len=3) :: matcst
    character(len=16) :: nmater, necoul, necris, necrci
    character(len=16) :: nomfam
    character(len=24) :: cpmono(5*nmat+1)
    integer(kind=8) :: i, nbfsys, ifa, j, dimtms, itbint
    integer(kind=8) :: nbsyst, nbsys
    character(len=8) :: fami
    integer(kind=8) :: jvMaterCode, kpg, ksp
    character(len=16) :: elasKeyword
    integer(kind=8) :: elasID
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to material parameters
    jvMaterCode = materPara%jvMaterCode
    elasKeyword = materPara%elasKeyword
    elasID = materPara%elasID
    fami = materPara%schemePara%fami
    kpg = materPara%schemePara%kpg
    ksp = materPara%schemePara%ksp

! - NB DE COMPOSANTES / VARIABLES INTERNES
    if (typmod1(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
    else if (typmod1(1:6) .eq. 'D_PLAN' .or. typmod1(1:4) .eq. 'AXIS') then
        ndt = 6
        ndi = 3
    else if (typmod1(1:6) .eq. 'C_PLAN') then
        ndt = 6
        ndi = 3
    end if
    dimtms = nfs*nsg*6
    call r8inir(dimtms, 0.d0, toutms, 1)
    materd = 0.d0
    materf = 0.d0
!
    call lcmmjv(multComp, nmat, cpmono, nbfsys, irota, &
                itbint, nsg, hsr)
    if (impexp .eq. 1) then
        if (irota .ne. 0) then
            call utmess('F', 'COMPOR2_11')
        end if
    end if

! - Compute matrix for rotation
    call matrot(materPara%lcsPara%lcsAngle, pgl)
!
    do i = 1, nmat
        do j = 1, 3
            nbcomm(i, j) = 0
        end do
    end do
    nbcomm(1, 1) = 1
!
    do ifa = 1, nbfsys
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        call lcmmsg(nomfam, nbsys, 0, pgl, ms, &
                    ng, lg, 0, q)
!
        nmater = cpmono(5*(ifa-1)+2) (1:16)
        necoul = cpmono(5*(ifa-1)+3) (1:16)
        necris = cpmono(5*(ifa-1)+4) (1:16)
        necrci = cpmono(5*(ifa-1)+5) (1:16)
!
!        COEFFICIENTS MATERIAUX LIES A L'ECOULEMENT
        call lcmafl(fami, kpg, ksp, '-', nmater, &
                    jvMaterCode, necoul, nbval, valres, nmat, &
                    itbint, nfs, nsg, hsr, nbsys)
        nvini = nbcomm(ifa, 1)
        if (necoul .eq. 'MONO_DD_KR') then
            nbval = nbval+1
!           une seule matrice d'interaction pour le monocristal
            valres(nbval) = 1
        end if
        do i = 1, nbval
            materd(nvini-1+i, 2) = valres(i)
        end do
        nbcomm(ifa, 2) = nvini+nbval
!
!        COEFFICIENTS MATERIAUX LIES A L'ECROUISSAGE CINEMATIQUE
        call lcmaec(fami, kpg, ksp, '-', nmater, &
                    jvMaterCode, necrci, nbval, valres, nmat)
        nvini = nbcomm(ifa, 2)
        do i = 1, nbval
            materd(nvini-1+i, 2) = valres(i)
        end do
        nbcomm(ifa, 3) = nvini+nbval
!
!        COEFFICIENTS MATERIAUX LIES A L'ECROUISSAGE ISOTROPE
        call lcmaei(fami, kpg, ksp, '-', nmater, &
                    jvMaterCode, necris, necoul, nbval, valres, &
                    nmat, itbint, nfs, nsg, hsr, &
                    ifa, nomfam, nbsys)
        nbval = nbval+1
!        une seule matrice d'interaction pour le monocristal
        valres(nbval) = 1
        nvini = nbcomm(ifa, 3)
        do i = 1, nbval
            materd(nvini-1+i, 2) = valres(i)
        end do
        nbcomm(ifa+1, 1) = nvini+nbval
!
!
    end do
!     ON STOCKE A LA FIN LE NOMBRE TOTAL DE COEF MATERIAU
    nbcomm(nmat, 2) = nbfsys
    nbcomm(nmat, 3) = nbcomm(nbfsys+1, 1)+1
    nbcomm(1, 1) = 1
!
    nbsyst = 0
!
    do ifa = 1, nbfsys
!
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        call lcmmsg(nomfam, nbsys, 0, pgl, ms, &
                    ng, lg, 0, q)
        nmater = cpmono(5*(ifa-1)+2) (1:16)
        necoul = cpmono(5*(ifa-1)+3) (1:16)
        necris = cpmono(5*(ifa-1)+4) (1:16)
        necrci = cpmono(5*(ifa-1)+5) (1:16)
!
        nbsyst = nbsyst+nbsys
!
        call lcmafl(fami, kpg, ksp, '+', nmater, &
                    jvMaterCode, necoul, nbval, valres, nmat, &
                    itbint, nfs, nsg, hsr, nbsys)
        nvini = nbcomm(ifa, 1)
        if (necoul .eq. 'MONO_DD_KR') then
            nbval = nbval+1
!           une seule matrice d'interaction pour le monocristal
            valres(nbval) = 1
        end if
        do i = 1, nbval
            materf(nvini-1+i, 2) = valres(i)
        end do
        nbcomm(ifa, 2) = nvini+nbval
!
        call lcmaec(fami, kpg, ksp, '+', nmater, &
                    jvMaterCode, necrci, nbval, valres, nmat)
        nvini = nbcomm(ifa, 2)
        do i = 1, nbval
            materf(nvini-1+i, 2) = valres(i)
        end do
        nbcomm(ifa, 3) = nvini+nbval
!
        call lcmaei(fami, kpg, ksp, '+', nmater, &
                    jvMaterCode, necris, necoul, nbval, valres, &
                    nmat, itbint, nfs, nsg, hsr, &
                    ifa, nomfam, nbsys)
        nvini = nbcomm(ifa, 3)
        nbval = nbval+1
!        une seule matrice d'interaction pour le monocristal
        valres(nbval) = 1
        do i = 1, nbval
            materf(nvini-1+i, 2) = valres(i)
        end do
        nbcomm(ifa+1, 1) = nvini+nbval
    end do

    if (elasID .eq. ELAS_ISOT) then
        propName(1) = 'E'
        propName(2) = 'NU'
        propName(3) = 'ALPHA'

! -     RECUPERATION MATERIAU A TEMPD (T)
        call rcvalb(fami, kpg, ksp, '-', jvMaterCode, &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    2, propName(1), materd(1, 1), cerr(1), 1)
        call rcvalb(fami, kpg, ksp, '-', jvMaterCode, &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    1, propName(3), materd(3, 1), cerr(3), 0)
        if (cerr(3) .ne. 0) materd(3, 1) = 0.d0
        materd(nmat, 1) = 0

! -     RECUPERATION MATERIAU A TEMPF (T+DT)
        call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    2, propName(1), materf(1, 1), cerr(1), 1)
        call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    1, propName(3), materf(3, 1), cerr(3), 0)
        if (cerr(3) .ne. 0) materf(3, 1) = 0.d0
        materf(nmat, 1) = 0

    else if (elasID .eq. ELAS_ORTH) then
        call dmat3d(materPara, '-', r8vide(), hooke)
        call d1ma3d(materPara, '-', r8vide(), ekooh)
        do j = 4, 6
            do i = 1, 6
                hooke(i, j) = hooke(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                hooke(i, j) = hooke(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 4, 6
            do i = 1, 6
                ekooh(i, j) = ekooh(i, j)/sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                ekooh(i, j) = ekooh(i, j)/sqrt(2.d0)
            end do
        end do
        do i = 1, 6
            do j = 1, 6
                materd(6*(j-1)+i, 1) = hooke(i, j)
                materd(36+6*(j-1)+i, 1) = ekooh(i, j)
            end do
        end do
        materd(nmat, 1) = 1
        propName(1) = 'ALPHA_L'
        propName(2) = 'ALPHA_T'
        propName(3) = 'ALPHA_N'
        call rcvalb(fami, kpg, ksp, '-', jvMaterCode, &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    3, propName, materd(73, 1), cerr, 0)
        if (cerr(1) .ne. 0) materd(73, 1) = 0.d0
        if (cerr(2) .ne. 0) materd(74, 1) = 0.d0
        if (cerr(3) .ne. 0) materd(75, 1) = 0.d0

        call dmat3d(materPara, '+', r8vide(), hookef)
        call d1ma3d(materPara, '+', r8vide(), fekooh)
        do j = 4, 6
            do i = 1, 6
                hookef(i, j) = hookef(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                hookef(i, j) = hookef(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 4, 6
            do i = 1, 6
                fekooh(i, j) = fekooh(i, j)/sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                fekooh(i, j) = fekooh(i, j)/sqrt(2.d0)
            end do
        end do
        do i = 1, 6
            do j = 1, 6
                materf(6*(j-1)+i, 1) = hookef(i, j)
                materf(36+6*(j-1)+i, 1) = fekooh(i, j)
            end do
        end do
        materf(nmat, 1) = 1
        call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    3, propName, materf(73, 1), cerr, 0)
        if (cerr(1) .ne. 0) materf(73, 1) = 0.d0
        if (cerr(2) .ne. 0) materf(74, 1) = 0.d0
        if (cerr(3) .ne. 0) materf(75, 1) = 0.d0
    else
        ASSERT(ASTER_FALSE)
    end if
!
    nr = ndt+nbsyst
    call calcmm(nbcomm, cpmono, nmat, pgl, nfs, &
                nsg, toutms, nvi, vind, &
                irota)

! - MATERIAU CONSTANT ?
    matcst = 'OUI'
    epsi = 1.d-3
    do i = 1, nmat
        if (abs(materd(i, 1)-materf(i, 1)) .gt. epsi*materd(i, 1)) then
            matcst = 'NON'
            goto 999
        end if
    end do
    do i = 1, nmat
        if (abs(materd(i, 2)-materf(i, 2)) .gt. epsi*materd(i, 2)) then
            matcst = 'NON'
            call utmess('F', 'COMPOR1_27')
            goto 999
        end if
    end do
!
999 continue
!
!     ON STOCKE A LA FIN LE NOMBRE TOTAL DE COEF MATERIAU
!      MATERD(NMAT,2)=NBCOMM(NMAT,3)
!      MATERF(NMAT,2)=NBCOMM(NMAT,3)
!
    call jedema()
!
end subroutine
