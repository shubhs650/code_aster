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
subroutine te0139(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmdlog.h"
#include "asterfort/nmgpfi.h"
#include "asterfort/nmgrla.h"
#include "asterfort/nmplxd.h"
#include "asterfort/nmtstm.h"
#include "asterfort/tecach.h"
#include "asterfort/tgveri_use.h"
#include "asterfort/tgveri.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "FE_module.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
!           3D_SI (HEXA20)
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuad
    type(FE_basis) :: FEBasis
!
    character(len=8) :: typmod(2)
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8) :: sz_tens, ndim
    integer(kind=8) :: nno, npg, imatuu, lgpg, iret
    integer(kind=8) :: jvGeom, jvMaterc, iuse
    integer(kind=8) :: icontm, ivarim
    integer(kind=8) :: jvInstmr, jvInstpr, ideplm, ideplp, jvCarcri
    integer(kind=8) :: ivectu, icontp, ivarip
    integer(kind=8) :: ivarix
    integer(kind=8) :: jtab(7)
    aster_logical :: matsym
    character(len=16), pointer :: jvCompor(:) => null(), mulcom(:) => null()
    character(len=16) :: multComp, defoComp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: codret
    integer(kind=8) :: jv_codret
!     POUR TGVERI
    real(kind=8) :: sdepl(3*MT_NNOMAX3D), svect(3*MT_NNOMAX3D), scont(6*MT_NNOMAX3D)
    real(kind=8) :: epsilo, disp_curr(MAX_BV_CG)
    real(kind=8), pointer :: varia(:) => null(), smatr(:) => null()
    blas_int :: b_incx, b_incy, b_n
    type(Behaviour_Integ) :: BEHInteg
    type(Material_Para) :: materPara
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    ivarip = 1
    imatuu = 1
    ivectu = 1
    ivarix = 1
    jv_codret = 1
    codret = 0

! - Set objects for finite element
    call FECell%init()
    nno = FECell%nbnodes
    ASSERT(nno .le. 27)
    ndim = FECell%ndim
    sz_tens = 2*ndim

! - Type of finite element
    if (ndim == 3) then
        typmod(1) = '3D'
    elseif (ndim == 2) then
        if (lteatt('AXIS', 'OUI')) then
            typmod(1) = 'AXIS'
        else if (lteatt('C_PLAN', 'OUI')) then
            typmod(1) = 'C_PLAN'
        else if (lteatt('D_PLAN', 'OUI')) then
            typmod(1) = 'D_PLAN'
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
    end if
    typmod(2) = ' '

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Properties of behaviour
    call jevech('PCOMPOR', 'L', vk16=jvCompor)
    call jevech('PCARCRI', 'L', jvCarcri)
    call jevech('PMULCOM', 'L', vk16=mulcom)
    multComp = mulcom(1)
    defoComp = jvCompor(DEFO)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              jvCompor, zr(jvCarcri), &
                              zr(jvInstmr), zr(jvInstpr), &
                              materPara, BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, jvCompor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Update displacements
    if (defoComp == "PETIT_REAC") then
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ideplm), b_incx, disp_curr, b_incy)
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, zr(ideplp), b_incx, disp_curr, b_incy)
        call FECell%updateCoordinates(disp_curr)
    end if

! - Init quadrature
    call FEQuad%initCell(FECell, fami)
    npg = FEQuad%nbQuadPoints

! - Init cell
    call FEBasis%initCell(FECell)

! - Get output fields
    if (lMatr) then
        call nmtstm(zr(jvCarcri), imatuu, matsym)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if
    if (option .eq. 'RIGI_MECA_IMPLEX') then
        call jevech('PCONTXR', 'E', icontp)
        b_n = to_blas_int(npg*sz_tens)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(icontm), b_incx, zr(icontp), b_incy)
    end if

! - Calcul de la matrice TGTE par PERTURBATION
    call tgveri_use(option, zr(jvCarcri), jvCompor, iuse)
    if (iuse == 1) then
        allocate (varia(2*3*MT_NNOMAX3D*3*MT_NNOMAX3D))
        allocate (smatr(3*MT_NNOMAX3D*3*MT_NNOMAX3D))
    end if

! - Update displacements
    if (defoComp == "PETIT_REAC") then
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ideplm), b_incx, disp_curr, b_incy)
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, zr(ideplp), b_incx, disp_curr, b_incy)
        call FECell%updateCoordinates(disp_curr)
    end if
!
500 continue
!
    if (defoComp .eq. 'PETIT') then
        call nmplxd(FECell, FEBasis, FEQuad, &
                    nno, npg, ndim, &
                    typmod, option, &
                    jvCompor, zr(jvCarcri), multComp, &
                    BEHInteg, &
                    zr(jvInstmr), zr(jvInstpr), &
                    zr(ideplm), zr(ideplp), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    matsym, zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999

    else if (defoComp .eq. 'PETIT_REAC') then
        call nmplxd(FECell, FEBasis, FEQuad, &
                    nno, npg, ndim, &
                    typmod, option, &
                    jvCompor, zr(jvCarcri), multComp, &
                    BEHInteg, &
                    zr(jvInstmr), zr(jvInstpr), &
                    zr(ideplm), zr(ideplp), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    matsym, zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999

    else if (defoComp .eq. 'SIMO_MIEHE') then
        call nmgpfi(BEHInteg, &
                    typmod, option, &
                    nno, npg, ndim, zr(jvGeom), &
                    jvCompor, zr(jvCarcri), multComp, &
                    zr(jvInstmr), zr(jvInstpr), &
                    zr(ideplm), zr(ideplp), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    zr(ivectu), zr(imatuu), codret)
        if (codret .ne. 0) goto 999

    else if (defoComp .eq. 'GREEN_LAGRANGE') then
        call nmgrla(FECell, FEBasis, FEQuad, &
                    nno, npg, ndim, &
                    typmod, option, &
                    jvCompor, zr(jvCarcri), multComp, &
                    BEHInteg, &
                    zr(jvInstmr), zr(jvInstpr), &
                    zr(ideplm), zr(ideplp), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    matsym, zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999

    else if (defoComp .eq. 'GDEF_LOG') then
        call nmdlog(FECell, FEBasis, FEQuad, &
                    nno, npg, ndim, &
                    typmod, option, &
                    jvCompor, zr(jvCarcri), multComp, &
                    BEHInteg, &
                    zr(jvInstmr), zr(jvInstpr), &
                    zr(ideplm), zr(ideplp), &
                    lgpg, zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    matsym, zr(ivectu), zr(imatuu), &
                    codret)
        if (codret .ne. 0) goto 999

    else
        ASSERT(ASTER_FALSE)
    end if

! - Calcul eventuel de la matrice TGTE par PERTURBATION
    call tgveri(option, zr(jvCarcri), jvCompor, nno, zr(jvGeom), &
                ndim, ndim*nno, zr(ideplp), sdepl, zr(ivectu), &
                svect, sz_tens*npg, zr(icontp), scont, npg*lgpg, &
                zr(ivarip), zr(ivarix), zr(imatuu), smatr, matsym, &
                epsilo, varia, iret)
    if (iret .ne. 0) then
        goto 500
    end if
!
999 continue

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', jv_codret)
        zi(jv_codret) = codret
    end if

! - Free large arrays
    if (iuse == 1) then
        deallocate (smatr)
        deallocate (varia)
    end if
!
end subroutine
