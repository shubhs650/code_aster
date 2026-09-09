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
subroutine te0123(option, nomte)
!
    use calcul_module, only: ca_jelvoi_, ca_jptvoi_, ca_jrepe_
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/massup.h"
#include "asterfort/nmplgs.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/voiuti.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN_GRAD_SIGM
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          MASS_MECA_*
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: famiRigi = "RIGI"
    integer(kind=8), parameter :: nvoima = 12, nscoma = 4
    integer(kind=8) :: dlns
    integer(kind=8) :: nno, nnob, nnos, npg, imatuu, lgpg, lgpg2
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom, jvMaterc
    integer(kind=8) :: ivfb, idfdeb
    integer(kind=8) :: icontm, ivarim
    integer(kind=8) :: iinstm, iinstp, idplgm, iddplg, jvCarcri
    integer(kind=8) :: ivectu, icontp, ivarip
    integer(kind=8) :: ivarix
    integer(kind=8) :: jtab(7), iadzi, iazk24, icoret, codret
    integer(kind=8) :: ndim, iret, ntrou
    character(len=16) :: codvoi
    integer(kind=8) :: nbvois
    integer(kind=8) :: livois(1:nvoima), tyvois(1:nvoima), nbnovo(1:nvoima)
    integer(kind=8) :: nbsoco(1:nvoima), lisoco(1:nvoima, 1:nscoma, 1:2)
    integer(kind=8) :: numa
    integer(kind=8) :: icodr1(1)
    character(len=8) :: typmod(2), lielrf(10)
    character(len=16) :: relaComp, defoComp
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: lVect, lMatr, lVari, lSigm, lMass
    blas_int :: b_incx, b_incy, b_n
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    ivarip = 1
    imatuu = 1
    ivectu = 1
    ivarix = 1
    icoret = 1
    codret = 0
!
    lMass = option(1:9) .eq. 'MASS_MECA'

! - Get element parameters
    call elref2(nomte, 10, lielrf, ntrou)
    ASSERT(ntrou .ge. 2)
    if (lMass) then
        call elrefe_info(elrefe=lielrf(1), fami='MASS', ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    else
        call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)
        call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nnob, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivfb, jdfde=idfdeb)
    end if

! - Type of finite element
    if (ndim .eq. 2 .and. lteatt('C_PLAN', 'OUI')) then
        typmod(1) = 'C_PLAN  '
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN  '
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        ASSERT(ndim .eq. 3)
    end if
    typmod(2) = 'GRADSIGM'

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMaterc)
!
    if (lMass) then
        call jevech('PMATUUR', 'E', imatuu)
        if (ndim .eq. 2) then
! - 2 DEPLACEMENTS + 4 DEF
            dlns = 6
        else if (ndim .eq. 3) then
! - 3 DEPLACEMENTS + 6 DEF
            dlns = 9
        else
            ASSERT(ndim .eq. 3)
        end if
        call massup(zi(jvMaterc), &
                    option, ndim, dlns, nno, nnos, &
                    npg, ipoids, idfde, &
                    zr(jvGeom), zr(ivf), imatuu, icodr1, jvGeom, &
                    ivf)
    else
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', idplgm)
        call jevech('PDEPLPR', 'L', iddplg)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
        ASSERT(jtab(1) .eq. ivarim)
        lgpg = max(jtab(6), 1)*jtab(7)

! ----- Get fields for behaviour
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', jvCarcri)
        relaComp = compor(RELA_NAME)
        defoComp = compor(DEFO)
        if (relaComp .ne. 'ENDO_HETEROGENE') then
            call utmess('F', 'COMPOR2_13')
        end if
        if (defoComp .ne. 'PETIT') then
            call utmess('F', 'ELEMENTS3_16', sk=defoComp)
        end if

! ----- Initialisation of behaviour datastructure
        call behaviourInit(BEHInteg)

! ----- Initializations of material parameters on current cell
        call initParaCell(famiRigi, zi(jvMaterc), materPara)

! ----- Set local coordinate system from user
        call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  zr(iinstm), zr(iinstm), &
                                  materPara, BEHInteg)

! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

! ----- Get output fields
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
        if (lVari) then
            call tecach('OOO', 'PVARIPR', 'E', iret, nval=7, itab=jtab)
            lgpg2 = max(jtab(6), 1)*jtab(7)
            call jevech('PVARIPR', 'E', ivarip)
            call jevech('PVARIMP', 'L', ivarix)
            b_n = to_blas_int(npg*lgpg2)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
        end if
        if (lSigm) then
            call jevech('PCONTPR', 'E', icontp)
        end if
        if (lVari) then
            ASSERT(lgpg .eq. lgpg2)
        end if

! ----- HYPO-ELASTICITE
        call tecael(iadzi, iazk24)
        numa = zi(iadzi-1+1)
        codvoi = 'A2'
!
        call voiuti(numa, codvoi, nvoima, nscoma, ca_jrepe_, &
                    ca_jptvoi_, ca_jelvoi_, nbvois, livois, tyvois, &
                    nbnovo, nbsoco, lisoco)
! ----- Compute
        call nmplgs(BEHInteg, &
                    ndim, nno, nnob, npg, &
                    zr(ivf), idfde, &
                    zr(ivfb), idfdeb, &
                    ipoids, zr(jvGeom), &
                    typmod, option, compor, zr(jvCarcri), &
                    zr(iinstm), zr(iinstp), &
                    zr(idplgm), zr(iddplg), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    zr(imatuu), zr(ivectu), &
                    codret, &
                    livois, nbvois, &
                    numa, lisoco, nbsoco, &
                    lVari, lSigm, lMatr, lVect)

! ----- Save return code
        if (lSigm) then
            call jevech('PCODRET', 'E', icoret)
            zi(icoret) = codret
        end if
    end if
end subroutine
