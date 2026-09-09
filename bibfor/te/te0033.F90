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
subroutine te0033(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    use plateMaterial_module, only: chckMultiLayer
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/cosiro.h"
#include "asterfort/dkqedg.h"
#include "asterfort/dkqsie.h"
#include "asterfort/dktedg.h"
#include "asterfort/dktsie.h"
#include "asterfort/dsqedg.h"
#include "asterfort/dsqsie.h"
#include "asterfort/dstedg.h"
#include "asterfort/dstsie.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxsiro.h"
#include "asterfort/dxsit2.h"
#include "asterfort/dxsit3.h"
#include "asterfort/dxsith.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/q4gedg.h"
#include "asterfort/q4gsie.h"
#include "asterfort/t3gedg.h"
#include "asterfort/t3gsie.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT/DKTG/DST/Q4G/Q4GG
!
! Options: SIEF_ELGA/EPSI_ELGA/DEGE_ELNO/DEGE/ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, npg, nbLayer
    integer(kind=8) :: jvDisp, jeffg, jvGeom, jvMaterc, jsigm
    integer(kind=8) :: multic
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: depl(24)
    real(kind=8) :: effgt(32), effpg(32)
    character(len=8) :: fami
    aster_logical :: lComposite
    type(Material_Para) :: materPara
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    if (option(6:9) .eq. 'ELNO') then
        fami = 'NOEU'
    else
        fami = 'RIGI'
    end if
    call elrefe_info(fami=fami, &
                     ndim=ndim, nno=nno, npg=npg)
!
    if (option .ne. 'SIEF_ELGA' .and. option .ne. 'EPSI_ELGA' .and. &
        option .ne. 'DEGE_ELNO' .and. option .ne. 'DEGE_ELGA') then
        ASSERT(ASTER_FALSE)
    end if
!
    effgt = 0.d0

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Material parameters
    lComposite = ASTER_FALSE
    if (option .eq. 'SIEF_ELGA' .or. option .eq. 'EPSI_ELGA') then
! ----- Get material parameters
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), materPara)
        lComposite = materPara%elasID .eq. ELAS_COMPOSITE

! ----- Set local coordinate system from user
        call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

    end if

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - For multi-layers
    if (option .eq. 'SIEF_ELGA' .or. option .eq. 'EPSI_ELGA') then
        if (lComposite) then
            call chckMultiLayer(materPara, plateCara)
        end if
    end if
    nbLayer = plateCara%nbLayer
    ASSERT(nbLayer .ge. 1)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry and displacements
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)
    call utpvgl(nno, 6, pgl, zr(jvDisp), depl)
!
    if (option(1:9) .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jsigm)
        if (nomte .eq. 'MEDKTR3') then
            call dktsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEDSTR3') then
            call dstsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEDKQU4') then
            call dkqsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEDSQU4') then
            call dsqsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MET3TR3') then
            call t3gsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else
            ASSERT(ASTER_FALSE)
        end if

        if (materPara%elasID .eq. ELAS_ISOT .or. &
            materPara%elasID .eq. ELAS_ORTH .or. &
            materPara%elasID .eq. ELAS_ISTR) then
            call dxsith(plateCara, &
                        materPara, zr(jsigm))

        else if (materPara%elasID .eq. ELAS_COMPOSITE) then
            call dxsit2(plateCara, plateOrie, &
                        zr(jsigm))

        elseif (materPara%elasID .eq. ELAS_SHELL) then
            call dxsit3(plateCara, plateOrie, &
                        zi(jvMaterc), zr(jsigm))

        else
            call utmess('F', 'PLATE1_1', nk=2, valk=[option, materPara%elasKeyword])
        end if

    else if (option(1:9) .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jsigm)
        if (nomte .eq. 'MEDKTR3') then
            call dktsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEDSTR3') then
            call dstsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEDKQU4') then
            call dkqsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEDSQU4') then
            call dsqsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        else if (nomte .eq. 'MET3TR3') then
            call t3gsie(plateCara, plateOrie, &
                        option, fami, xyzl, depl, &
                        zr(jsigm))
        end if
        call dxsiro(npg*plateCara%nbLayer*3, plateOrie%t2iu, zr(jsigm), zr(jsigm))

    else if (option(1:9) .eq. 'DEGE_ELNO') then
        call jevech('PDEFOGR', 'E', jeffg)
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effgt, multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effgt)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effgt)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effgt)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effgt)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effgt)
        end if
        call dxefro(nno, plateOrie%t2iu, effgt, zr(jeffg))

    else if (option(1:9) .eq. 'DEGE_ELGA') then
        call jevech('PDEFOPG', 'E', jeffg)
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effpg, multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effpg)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effpg)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effpg)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effpg)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(plateCara, plateOrie, &
                        xyzl, option, depl, &
                        effpg)
        end if
        call dxefro(npg, plateOrie%t2iu, effpg, zr(jeffg))

    end if
!
    if (option .eq. 'SIEF_ELGA') then
        call cosiro(plateCara, plateOrie, &
                    'PCONTRR', 'E', 'IU', 'G', &
                    jsigm)
    end if
!
end subroutine
