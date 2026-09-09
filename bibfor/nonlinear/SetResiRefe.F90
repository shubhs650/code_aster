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

subroutine SetResiRefe(mesh, ds_conv)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/r8vide.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/getelem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"

    character(len=8), intent(in):: mesh
    type(NL_DS_Conv), intent(inout) :: ds_conv
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Set values for reference residual (by name)
!
! --------------------------------------------------------------------------------------------------
! IN  mesh             : mesh name
! IO  ds_conv          : datastructure for convergence management
! --------------------------------------------------------------------------------------------------
    character(len=16), parameter :: motclf = 'CONVERGENCE_REFE'
    character(len=24), parameter :: lst_mail = '&&SRR.MAIL'
    character(len=19), parameter :: cresicmp = '&&SRR.CRESICMP'
    character(len=19), parameter :: cresiref = '&&SRR.CRESIREF'
! --------------------------------------------------------------------------------------------------
    character(len=8):: k8b
    integer(kind=8):: gd, nocc, iocc, mc, nb_refe, nb_neut, nbma, nb
    real(kind=8) :: ref
    character(len=8), pointer:: cmp_names(:) => null()
    character(len=8), pointer:: neut_names(:) => null()
    character(len=8), pointer:: cresicmp_names(:) => null()
    character(len=8), pointer:: cresicmp_values(:) => null()
    character(len=8), pointer:: cresiref_names(:) => null()
    real(kind=8), pointer:: cresiref_values(:) => null()
    integer(kind=8), pointer:: numa(:) => null()
! --------------------------------------------------------------------------------------------------

    call jemarq()

    ASSERT(mesh .ne. ' ')

    ! Noms des cartes
    ds_conv%cresicmp = cresicmp
    ds_conv%cresiref = cresiref

    ! Noms des grandeurs de référence
    call jenonu(jexnom('&CATA.GD.NOMGD', 'RESIREF'), gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nb_refe)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', vk8=cmp_names)

    ! Carte des noms des grandeurs de référence
    call jenonu(jexnom('&CATA.GD.NOMGD', 'RESICMP'), gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nb_neut)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', vk8=neut_names)
    ASSERT(nb_neut .eq. nb_refe)

    call detrsd('CARTE', cresicmp)
    call alcart('V', cresicmp, mesh, 'RESICMP')
    call jeveuo(cresicmp//'.NCMP', 'E', vk8=cresicmp_names)
    call jeveuo(cresicmp//'.VALV', 'E', vk8=cresicmp_values)

    cresicmp_names = neut_names
    cresicmp_values = cmp_names
    call nocart(cresicmp, 1, nb_refe)
    call jedetr(cresicmp//'.NCMP')
    call jedetr(cresicmp//'.VALV')

    ! Carte des valeurs de référence par défaut
    call detrsd('CARTE', cresiref)
    call alcart('V', cresiref, mesh, 'RESIREF')
    call jeveuo(cresiref//'.NCMP', 'E', vk8=cresiref_names)
    call jeveuo(cresiref//'.VALV', 'E', vr=cresiref_values)

    cresiref_names = cmp_names
    cresiref_values = r8vide()
    call nocart(cresiref, 1, nb_refe)

    ! Carte des valeurs de référence fournies par l'utilisateur
    call getfac(motclf, nocc)
    if (nocc .le. 0) goto 800
    do iocc = 1, nocc

        ! Lecture des valeurs de référence présentes
        do mc = 1, nb_refe
            call getvr8(motclf, cmp_names(mc), iocc=iocc, scal=ref, nbret=nb)
            cresiref_values(mc) = merge(r8vide(), ref, nb .eq. 0)
        end do

        ! Lieu d'affectation
        call getvtx(motclf, 'TOUT', iocc=iocc, scal=k8b, nbret=nb)
        if (nb .ne. 0) then
            call nocart(cresiref, 1, nb_refe)
!
        else
            call getelem(mesh, motclf, iocc, ' ', lst_mail, nbma)
            if (nbma .ne. 0) then
                call jeveuo(lst_mail, 'L', vi=numa)
                call nocart(cresiref, 3, nb_refe, mode='NUM', nma=nbma, &
                            limanu=numa)
                call jedetr(lst_mail)
            end if
        end if
    end do
800 continue

    call jedetr(cresiref//'.NCMP')
    call jedetr(cresiref//'.VALV')
!
    call jedema()
end subroutine
