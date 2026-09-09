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

subroutine te0004(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ppgan2.h"
#include "asterfort/tecach.h"
    character(len=16) :: option, nomte
!
!     BUT:
!       !! ROUTINE GENERIQUE !!
!       ASSURER LE PASSAGE ELGA -> ELNO
!
!.......................................................................
!
!
!
    integer(kind=8) :: ndim, nno, nnos, npg, nbsp
    integer(kind=8) :: ipoids, ivf, idfde, jgano
    integer(kind=8) :: iret, nbcmp, itabin(7), itabou(7)
    integer(kind=8) :: iinpg, ioutno
!
    character(len=8) :: fami
!
! ----------------------------------------------------------------------
!
    if (option .eq. 'DERA_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDERAPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDERANO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'DIFF_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDIFFPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDIFFNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'DISS_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDISSPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDISSNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EFGE_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PEFGAR', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PEFFORR', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'ENDO_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PTRIAPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PTRIANO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'ENEL_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PENERPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PENERNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPFD_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPFP_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPME_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPMG_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPSG_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
    else if (option .eq. 'EPSL_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPSI_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPSP_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'EPVC_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEFOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEFONO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'ETOT_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PENERPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PENERNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'FLUX_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PFLUXPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PFLUXNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'GRAT_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PGRATPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PGRATNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'HYDR_ELNO') then
        fami = 'MASS'
        call tecach('OOO', 'PHYDRPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PHYDRNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'HYGR_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PHYGRPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PHYGRNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'SAUT_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PDEPOPG', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PDEPSNO', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'SIGM_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PSIEFNOR', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'SIEF_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PSIEFNOR', 'E', iret, nval=7, &
                    itab=itabou)
!
    else if (option .eq. 'VARC_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PVARCGR', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PVARCNR', 'E', iret, nval=7, &
                    itab=itabou)
    else if (option .eq. 'VARI_ELNO') then
        fami = 'RIGI'
        call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, &
                    itab=itabin)
        call tecach('OOO', 'PVARINR', 'E', iret, nval=7, &
                    itab=itabou)
!
    else
        ASSERT(.false.)
    end if
!
!    -------------------------------------------------------------------
!    -- COMMUN A TOUTES LES OPTIONS
!    -------------------------------------------------------------------
!
! --- IL FAUDRAIT PENSER A RECUPERER FAMI DE MANIERE AUTOMATIQUE ET SURE
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    iinpg = itabin(1)
    ioutno = itabou(1)
!
! --- NOMBRE DE COMPOSANTE
    nbcmp = itabin(2)/itabin(3)
    ASSERT(nbcmp .gt. 0)
    ASSERT(nbcmp .eq. itabou(2)/itabou(3))
!
! --- ON ESSAIE UNE AUTRE FAMILLE DE POINT DE GAUSS
    if (itabin(3) .ne. npg) then
        fami = 'MASS'
        call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    end if
    ASSERT(itabin(3) .eq. npg)
!
! --- NOMBRE DE SOUS-POINT
    nbsp = itabin(7)
    ASSERT(nbsp .eq. itabou(7))
!
! --- NOMBRE DE COMPOSANTE DYNAMIQUE
    if (option .eq. 'VARI_ELNO') then
        nbcmp = itabin(6)
        ASSERT(nbcmp .eq. itabou(6))
    end if
!
    call ppgan2(jgano, nbsp, nbcmp, zr(iinpg), zr(ioutno))
!
end subroutine
