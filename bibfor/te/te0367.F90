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
subroutine te0367(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gedisc.h"
#include "asterfort/jevech.h"
#include "asterfort/nmtstm.h"
#include "asterfort/nmspfm.h"
#include "asterfort/matrot.h"
#include "asterfort/spmats.h"
#include "asterfort/rccome.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Element: 3D_INTERF_POU
!
! Options: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ndim = 3
    character(len=4), parameter :: fami = "RIGI"
    integer(kind=8) :: nno, nno_s, nno_p, npg, nddl, nddl_s, nddl_p, nddlsym
    integer(kind=8) :: ipoids, ivf, ivfs, ivfp, icoopg
    integer(kind=8) :: iorie, igeom, imater, icarcr, idepm, iddep, icoret
    integer(kind=8) :: icontm, icontp, ivect, imatr
    integer(kind=8) :: ivarim, ivarip, iinstm, iinstp
    integer(kind=8) :: lgpg, codret, icodret
    real(kind=8) :: coopg(4, 4)
    real(kind=8) :: pgl(3, 3)
    character(len=8), parameter :: typmod(2) = (/'3D      ', 'INTSOLPI'/)
    character(len=16), dimension(COMPOR_SIZE) :: compor
    aster_logical :: matsym
    aster_logical :: lVect, lMatr, lVari, lSigm, lElas
    type(Behaviour_Integ) :: BEHinteg
    type(Material_Para) :: materPara
    character(len=8) :: matint, matpou
!
! --------------------------------------------------------------------------------------------------
!
    ivarip = 1
    icoret = 1
    icontp = 1
    ivect = 1
    icoret = 1
    imatr = 1
    icarcr = 1
    icontm = 1
    idepm = 1
    iddep = 1
    ivarim = 1
    iinstp = 1
    lgpg = 1

    ! - Get general pointer info on element
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf)

    nno_s = 8
    nno_p = 2
    nddl_s = 3*nno_s
    nddl_p = 6*nno_p
    nddl = nddl_s+nddl_p
    ivfs = ivf
    ivfp = icoopg

    ASSERT(nno .eq. 27)
    ASSERT(npg .eq. 2 .or. npg .eq. 3 .or. npg .eq. 4)

! - Get input fields
    call jevech('PCAORIE', 'L', iorie)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imater)
    call jevech('PINSTR', 'L', iinstm)

! - Get coordinates of gauss points
    call gedisc(3, nno, npg, zr(ivf), zr(igeom), &
                coopg)

! - Get multiple materials
    call spmats(imater, matint, matpou)
    call rccome(matint, 'INTERF_POU_ELAS', icodret)
    ASSERT(icodret .eq. 0)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(imater), materPara)
    materPara%matname = matint

! - Set local coordinate system from user
    ! call getUserLCS(ndim, nno, igeom, materPara%lcsPara)
    call initLCSNone(materPara)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Force behaviour
    compor(1:COMPOR_SIZE) = 'VIDE'
    compor(DEFO_LDC) = 'TOTALE'
    compor(RELA_NAME) = 'INTERF_POU_ELAS'
    compor(DEFO) = 'PETIT'
    compor(NUME) = '1'
    compor(NVAR) = '1'

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(icarcr), &
                              zr(iinstm), zr(iinstp), &
                              materPara, BEHinteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)
    lElas = ASTER_TRUE

! - Get output fields
    ASSERT(lMatr)
    call nmtstm(zr(icarcr), imatr, matsym)
    if (matsym) then
        nddlsym = nddl*(nddl+1)/2
    else
        nddlsym = nddl*nddl
    end if
    zr(imatr:imatr-1+nddlsym) = 0.d0

! - Get orientation
    call matrot(zr(iorie), pgl)

! - Main calculation
    call nmspfm(BEHinteg, typmod, ndim, nno, nddl, nddlsym, &
                nno_p, nno_s, nddl_p, nddl_s, npg, lgpg, &
                zr(ipoids), zr(ivfs), zr(ivfp), &
                pgl, zr(igeom), 3, zi(imater), matint, matpou, option, &
                zr(idepm), zr(iddep), zr(icontm), zr(icontp), zr(ivect), &
                zr(imatr), zr(ivarim), zr(ivarip), zr(icarcr), compor, &
                zr(iinstm), zr(iinstp), coopg, matsym, lMatr, lVect, lSigm, lElas, &
                codret)

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', icoret)
        zi(icoret) = codret
    end if
!
end subroutine
