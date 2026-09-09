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
subroutine te0039(option, nomte)
!
    use te0047_type
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dikpkt.h"
#include "asterfort/dis_choc_frot_syme.h"
#include "asterfort/dis_choc_frot_nosyme.h"
#include "asterfort/dis_elas_para_klfl.h"
#include "asterfort/discret_sief.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/lteatt.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "blas/dcopy.h"
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! IN OPTION    : K16 :  OPTION DE CALCUL
!     'SIEQ_ELNO'  'SIEQ_ELGA'
!     'FORC_NODA'  'REFE_FORC_NODA'
!     'FONL_NOEU'
! IN NOMTE     : K16 : NOM DU TYPE ELEMENT
!     DISCRETS :
!        'MECA_DIS_T_N'      'MECA_DIS_T_L'     'MECA_DIS_TR_N'
!        'MECA_DIS_TR_L'     'MECA_2D_DIS_T_N'  'MECA_2D_DIS_T_L'
!        'MECA_2D_DIS_TR_N'  'MECA_2D_DIS_TR_L'
! --------------------------------------------------------------------------------------------------
!
    type(te0047_dscr) :: for_discret
!
!   Les variables internes
    integer(kind=8), parameter :: nbvari = 10
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
!
    integer(kind=8) :: lorien, lmater, ii, jj, i, j, k
    integer(kind=8) :: ivectu, jvSief
    integer(kind=8) :: iplouf, infodi, ibid
    integer(kind=8) :: igeom, jdc, irep, ifono, ilogic, jvDisp
!
    real(kind=8) :: pgl(3, 3), force(3)
    real(kind=8) :: fs(12), ugp(12), ulp(12), dpe(12)
    real(kind=8) :: sim(12), sip(12), fono(12)
    real(kind=8) :: forref, momref
    real(kind=8) :: r8bid
    real(kind=8) :: utotxyz(3), flp(3)
    real(kind=8) :: kp, kt1, kt2
    real(kind=8), pointer :: klvp(:) => null()
    real(kind=8), pointer :: klv(:) => null()
    real(kind=8), pointer :: kgv(:) => null()
!
    character(len=8) :: k8bid
    character(len=16) :: kmess(5)
    character(len=16), pointer :: compor(:) => null()
!
    aster_logical, parameter :: Predic = ASTER_FALSE
    aster_logical :: lMatrTangSyme, IsSymetrique
    blas_int :: b_incx, b_incy, b_n
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
    infodi = 1
!   On verifie que les caracteristiques ont ete affectees
!   Le code du discret
    call infdis('CODE', ibid, r8bid, nomte)
!   Le code stoké dans la carte
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    end if
!   Discret de type raideur
    call infdis('DISK', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_27', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    end if

!   Récupere les informations sur les elements
    for_discret%option = option
    for_discret%nomte = nomte
!   Matrice de raideur symetrique ou pas, pour les discrets
    call infdis('SYMK', for_discret%syme, r8bid, k8bid)
    call getDiscretInformations(for_discret)
!
    if (option(1:14) .eq. 'REFE_FORC_NODA') then
        call refe%Init(nomte)
        call jevech('PVECTUR', 'E', ivectu)
        if (lteatt('MODELI', 'DTR')) then
            forref = refe%GetRef('EFFORT')
            momref = refe%GetRef('MOMENT')
            do ii = 1, for_discret%nno
                do jj = 1, 3
                    zr(ivectu+(ii-1)*for_discret%nc+jj-1) = forref
                end do
                do jj = 4, for_discret%nc
                    zr(ivectu+(ii-1)*for_discret%nc+jj-1) = momref
                end do
            end do
        else if (lteatt('MODELI', '2DT')) then
            forref = refe%GetRef('EFFORT')
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc) = forref
                zr(ivectu+(ii-1)*for_discret%nc+1) = forref
            end do
        else if (lteatt('MODELI', '2TR')) then
            forref = refe%GetRef('EFFORT')
            momref = refe%GetRef('MOMENT')
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc) = forref
                zr(ivectu+(ii-1)*for_discret%nc+1) = forref
                zr(ivectu+(ii-1)*for_discret%nc+2) = momref
            end do
        else if (lteatt('MODELI', 'DIT')) then
            forref = refe%GetRef('EFFORT')
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc) = forref
                zr(ivectu+(ii-1)*for_discret%nc+1) = forref
                zr(ivectu+(ii-1)*for_discret%nc+2) = forref
            end do
        else
            kmess(1) = option
            kmess(2) = nomte
            kmess(3) = 'TE0039'
            call utmess('F', 'DISCRETS_15', nk=2, valk=kmess)
        end if
        call refe%Check()
    else if (option .eq. 'FONL_NOEU') then
        AS_ALLOCATE(vr=klv, size=for_discret%nbt)
        klv = 0.d0
        AS_ALLOCATE(vr=kgv, size=for_discret%nbt)
        kgv = 0.d0

        call jevech('PGEOMER', 'L', igeom)
        call jevech('PDEPLAR', 'L', jvDisp)
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PMATERC', 'L', lmater)
        if (lteatt('MODELI', 'DTR') .or. lteatt('MODELI', 'DIT')) then
!   PARAMETRES EN ENTREE
            call jevech('PCAORIE', 'L', lorien)
            call matrot(zr(lorien), for_discret%pgl)
! DÉPLACEMENTS DANS LE REPÈRE GLOBAL
!   UGM = DEPLACEMENT PRECEDENT
!   UGP = DEPLACEMENT COURANT
            do ii = 1, for_discret%neq
                ugp(ii) = zr(jvDisp+ii-1)
            end do
! DÉPLACEMENTS DANS LE REPÈRE LOCAL
!   ULM = DEPLACEMENT PRECEDENT    = PLG * UGM
!   DUL = INCREMENT DE DEPLACEMENT = PLG * DUG
!   ULP = DEPLACEMENT COURANT      = PLG * UGP
            if (for_discret%ndim .eq. 3) then
                call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugp, ulp)
            else if (for_discret%ndim .eq. 2) then
                call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugp, ulp)
            end if
! Seul le cas symetrique est traite
            call infdis('SYMK', iplouf, r8bid, k8bid)
            if (iplouf .ne. 1) then
                kmess(1) = option
                kmess(2) = nomte
                kmess(3) = 'TE0039'
                kmess(4) = ' '
                call utmess('F', 'DISCRETS_12', nk=4, valk=kmess)
            end if
            !
            call jevech('PCADISK', 'L', jdc)
            call infdis('REPK', irep, r8bid, k8bid)
! irep = 1 = matrice en repère global
            if (irep .eq. 1) then
! Matrice dans le repère global
                b_n = to_blas_int(for_discret%nbt)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(jdc), b_incx, kgv, b_incy)
! Matrice dans le repère local
                call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, kgv, klv)
            else
! Matrice dans le repère local
                b_n = to_blas_int(for_discret%nbt)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
! Matrice dans le repère global
                call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, kgv)
            end if
            !
            if (compor(RELA_NAME) .eq. 'DIS_CHOC') then
                varmo(:) = 0.0; dpe(:) = 0.0
! Relation de comportement de choc : forces nodales
                call jevech('PVECTUR', 'E', ifono)
                do ii = 1, for_discret%neq
                    zr(ifono+ii-1) = 0.0
                    sim(ii) = 0.0
                end do
                !
                ilogic = 0; force(1:3) = 0.0
                call discret_sief(for_discret, klv, ulp, sim, ilogic, &
                                  sip, fono, force)

                ! Verification du caractère symétrique de la matrice tangente
                call hasSymmetricTangentMatrix(for_discret, lMatrTangSyme)

                !   Calcul des efforts
                if (lMatrTangSyme) then
                    call dis_choc_frot_syme(for_discret, zi(lmater), ulp, zr(igeom), klv, &
                                            kgv, dpe, Predic, &
                                            force, varmo, varpl)
                else
                    call dis_choc_frot_nosyme(for_discret, zi(lmater), ulp, zr(igeom), klv, &
                                              dpe, varmo, force, varpl)
                end if

                ! Ajout de la contribution d'un discret élastique en parallèle
                ! --- Déplacement total
                utotxyz(1:3) = 0.d0
                if (for_discret%nno .eq. 2) then
                    utotxyz(1) = ulp(1+for_discret%nc)-ulp(1)+dpe(1+for_discret%nc)-dpe(1)
                    utotxyz(2) = ulp(2+for_discret%nc)-ulp(2)+dpe(2+for_discret%nc)-dpe(2)
                    if (for_discret%ndim .eq. 3) then
                        utotxyz(3) = ulp(3+for_discret%nc)-ulp(3)+dpe(3+for_discret%nc)-dpe(3)
                    end if
                else
                    utotxyz(1) = ulp(1)+dpe(1)
                    utotxyz(2) = ulp(2)+dpe(2)
                    if (for_discret%ndim .eq. 3) then
                        utotxyz(3) = ulp(3)+dpe(3)
                    end if
                end if

                AS_ALLOCATE(vr=klvp, size=for_discret%nsym)
                klvp = 0.d0

                ! --- Raideurs du discret élastique
                call dikpkt(zi(lmater), 'DIS_CONTACT', kp, kt1, kt2)
                ! --- Matrice tangente et vecteur force du discret élastique (repère local)
                call dis_elas_para_klfl(for_discret, kp, kt1, kt2, utotxyz, klvp, flp)
                ! --- Ajout de la contribution du discret en parallèle
                call hasSymmetricTangentMatrix(for_discret, IsSymetrique)
                if (IsSymetrique) then
                    ASSERT(for_discret%nbt .eq. for_discret%nsym)
                    klv = klv+klvp
                else
                    ASSERT(for_discret%nbt .eq. for_discret%nfull)
                    k = 1
                    do i = 1, for_discret%neq
                        do j = 1, i
                            ii = for_discret%neq*(i-1)
                            jj = for_discret%neq*(j-1)
                            !
                            klv(ii+j) = klv(ii+j)+klvp(k)
                            if (i /= j) then
                                klv(jj+i) = klv(jj+i)+klvp(k)
                            end if
                            k = k+1
                        end do
                    end do
                end if
                !
                force = force+flp
                !
                ilogic = 2
                call discret_sief(for_discret, klv, ulp, sim, ilogic, &
                                  sip, zr(ifono), force)
                do ii = 1, for_discret%neq
                    zr(ifono+ii-1) = zr(ifono+ii-1)-fono(ii)
                end do
                if (for_discret%nno .eq. 2) then
                    do ii = 1, for_discret%nc
                        zr(ifono+ii-1) = 0.0
                    end do
                end if
                AS_DEALLOCATE(vr=klvp)
            end if
        end if
        AS_DEALLOCATE(vr=klv)
        AS_DEALLOCATE(vr=kgv)
    else if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PVECTUR', 'E', ivectu)
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%neq
                fs(ii) = zr(jvSief+ii-1)
            end do
        else
            do ii = 1, for_discret%nc
                fs(ii) = -zr(jvSief+ii-1)
                fs(ii+for_discret%nc) = zr(jvSief+ii+for_discret%nc-1)
            end do
        end if
        call jevech('PCAORIE', 'L', lorien)
        !
        call matrot(zr(lorien), pgl)
        if (for_discret%ndim .eq. 3) then
            call utpvlg(for_discret%nno, for_discret%nc, pgl, fs, zr(ivectu))
        else if (for_discret%ndim .eq. 2) then
            call ut2vlg(for_discret%nno, for_discret%nc, pgl, fs, zr(ivectu))
        end if
    else
        ASSERT(.false.)
    end if

!
end subroutine
