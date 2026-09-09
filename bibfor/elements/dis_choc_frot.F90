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
subroutine dis_choc_frot(DD, iret)
!
! --------------------------------------------------------------------------------------------------
!
! IN    DD      : voir l'appel
! OUT   iret    : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/dis_choc_frot_syme.h"
#include "asterfort/dis_choc_frot_nosyme.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpnlg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: DD
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jdc, irep, imat, ivarim, ii, ivitp, idepen, iviten, neq, igeom, ivarip, fdnc
    integer(kind=8) :: iretlc, ifono, icontm, icontp, iiter, iterat
!
!   Les variables internes
    integer(kind=8), parameter :: nbvari = 10
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
!
    real(kind=8) :: dvl(12), dpe(12), dve(12), ulp(12), fl(12), force(3)
    real(kind=8) :: klc(144), klv(78), kgv(78)
!
    real(kind=8) :: r8bid
    character(len=8) :: k8bid
!
    aster_logical :: IsSymetrique, IsDynamique, IsStatique, Prediction
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTMR', 'L', icontm)
!   on recupere le no de l'iteration de newton
    call jevech('PITERAT', 'L', iiter)
    iterat = zi(iiter)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   absolu vers local ? ---
!   irep = 1 = matrice en repère global
!       klv : matrice dans le repère local
!       kgv : matrice dans le repère global
    if (irep .eq. 1) then
! Matrice dans le repère global
        b_n = to_blas_int(DD%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, kgv, b_incy)
! Matrice dans le repère local
        if (DD%ndim .eq. 3) then
            call utpsgl(DD%nno, DD%nc, DD%pgl, kgv, klv)
        else
            call ut2mgl(DD%nno, DD%nc, DD%pgl, kgv, klv)
        end if
    else
! Matrice dans le repère local
        b_n = to_blas_int(DD%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
! Matrice dans le repère global
        if (DD%ndim .eq. 3) then
            call utpslg(DD%nno, DD%nc, DD%pgl, klv, kgv)
        else
            call ut2mlg(DD%nno, DD%nc, DD%pgl, klv, kgv)
        end if
    end if
    call jevech('PMATERC', 'L', imat)
    call jevech('PVARIMR', 'L', ivarim)
!
    do ii = 1, nbvari
        varmo(ii) = zr(ivarim+ii-1)
        varpl(ii) = varmo(ii)
    end do
!
    call tecach('ONO', 'PVITPLU', 'L', iretlc, iad=ivitp)
    if (iretlc .eq. 0) then
        IsDynamique = ASTER_TRUE
        if (DD%ndim .eq. 3) then
            call utpvgl(DD%nno, DD%nc, DD%pgl, zr(ivitp), dvl)
        else
            call ut2vgl(DD%nno, DD%nc, DD%pgl, zr(ivitp), dvl)
        end if
    else
        IsDynamique = ASTER_FALSE
        dvl(:) = 0
    end if
    IsStatique = .not. IsDynamique
!
    call tecach('ONO', 'PDEPENT', 'L', iretlc, iad=idepen)
    if (iretlc .eq. 0) then
        if (DD%ndim .eq. 3) then
            call utpvgl(DD%nno, DD%nc, DD%pgl, zr(idepen), dpe)
        else
            call ut2vgl(DD%nno, DD%nc, DD%pgl, zr(idepen), dpe)
        end if
    else
        dpe(:) = 0.0d0
    end if
!
    call tecach('ONO', 'PVITENT', 'L', iretlc, iad=iviten)
    if (iretlc .eq. 0) then
        if (DD%ndim .eq. 3) then
            call utpvgl(DD%nno, DD%nc, DD%pgl, zr(iviten), dve)
        else
            call ut2vgl(DD%nno, DD%nc, DD%pgl, zr(iviten), dve)
        end if
    else
        dve(:) = 0.d0
    end if
    !
    ulp(:) = DD%ulm(:)+DD%dul(:)
    !
!   Relation de comportement de choc
    call hasSymmetricTangentMatrix(DD, IsSymetrique)
    if (IsSymetrique) then
        ! Prédiction en dynamique, on retourne les efforts précédents
        Prediction = IsDynamique .and. (iterat .eq. 1) .and. (DD%option .eq. 'RAPH_MECA')
        ! Formulation symétrique
        call dis_choc_frot_syme(DD, zi(imat), ulp, zr(igeom), klv, &
                                kgv, dpe, Prediction, &
                                force, varmo, varpl)
    else
        ! Formulation non symétrique
        call dis_choc_frot_nosyme(DD, zi(imat), ulp, zr(igeom), klv, &
                                  dpe, varmo, force, varpl)
    end if
!   Actualisation de la matrice tangente
    if (DD%lMatr) then
        if (IsSymetrique) then
            call jevech('PMATUUR', 'E', imat)
            if (DD%ndim .eq. 3) then
                call utpslg(DD%nno, DD%nc, DD%pgl, klv, zr(imat))
            else
                call ut2mlg(DD%nno, DD%nc, DD%pgl, klv, zr(imat))
            end if
        else
            call jevech('PMATUNS', 'E', imat)
            call utpnlg(DD%nno, DD%nc, DD%pgl, klv, zr(imat))
        end if
    end if
!
    neq = DD%nno*DD%nc
    fdnc = DD%nc
!   calcul des efforts généralisés, des forces nodales et des variables internes
    if (DD%lVect .or. DD%lSigm) then
        if (IsSymetrique) then
!           Demi-matrice klv (triangulaire supérieure) transformée en matrice pleine klc
            call vecma(klv, DD%nbt, klc, neq)
!           Calcul de fl = klc.dul (DD%incrément d'effort)
            call pmavec('ZERO', neq, klc, DD%dul, fl)
        else
!           klv est déjà la matrice pleine
            call pmavec('ZERO', neq, klv, DD%dul, fl)
        end if
    end if
!   calcul des efforts généralisés
    if (DD%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (DD%nno .eq. 1) then
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (DD%nno .eq. 2) then
            do ii = 1, fdnc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+fdnc) = fl(ii+fdnc)+zr(icontm-1+ii+fdnc)
            end do
        end if
        if (DD%nno .eq. 1) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+2) = force(2)
            if (DD%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
            end if
        else if (DD%nno .eq. 2) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+1+fdnc) = force(1)
            zr(icontp-1+2) = force(2)
            zr(icontp-1+2+fdnc) = force(2)
            if (DD%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
                zr(icontp-1+3+fdnc) = force(3)
            end if
        end if
!       Si force de contact petite ==> Tout est nul
        if (abs(force(1)) .lt. r8prem()) then
            do ii = 1, neq
                zr(icontp-1+ii) = 0.0d0
            end do
        end if
    end if
! calcul des efforts généralisés et des forces nodales
    if (DD%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (DD%nno .eq. 1) then
            do ii = 1, neq
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (DD%nno .eq. 2) then
            do ii = 1, fdnc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+fdnc) = fl(ii+fdnc)+zr(icontm-1+ii+fdnc)
            end do
        end if
        if (DD%nno .eq. 1) then
            fl(1) = force(1)
            fl(2) = force(2)
            if (DD%ndim .eq. 3) then
                fl(3) = force(3)
            end if
        else if (DD%nno .eq. 2) then
            fl(1) = -force(1)
            fl(1+fdnc) = force(1)
            fl(2) = -force(2)
            fl(2+fdnc) = force(2)
            if (DD%ndim .eq. 3) then
                fl(3) = -force(3)
                fl(3+fdnc) = force(3)
            end if
        end if
!       Si force de contact petite ==> Tout est nul
        if (abs(force(1)) .lt. r8prem()) then
            fl(1:neq) = 0.0d0
        end if
!       forces nodales aux noeuds 1 et 2 (repère global)
        if (DD%ndim .eq. 3) then
            call utpvlg(DD%nno, DD%nc, DD%pgl, fl, zr(ifono))
        else
            call ut2vlg(DD%nno, DD%nc, DD%pgl, fl, zr(ifono))
        end if
    end if
!
!   Mise à jour des variables internes
    if (DD%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        do ii = 1, nbvari
            zr(ivarip+ii-1) = varpl(ii)
            if (DD%nno .eq. 2) zr(ivarip+ii-1+nbvari) = varpl(ii)
        end do
    end if
end subroutine
