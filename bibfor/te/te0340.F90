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
subroutine te0340(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cgfint.h"
#include "asterfort/cginit.h"
#include "asterfort/cgtang.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: CABLE_GAINE
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
    character(len=8), parameter :: fami = "RIGI"
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: comporKit(COMPOR_SIZE)
    character(len=8):: typmod(2)
    character(len=8) :: lielrf(10)
    integer(kind=8) :: nno1, nno2, npg, imatuu, lgpg, lgpg1, lgpg2
    integer(kind=8) :: iw, ivf1, idf1, jvGeom, jvMaterc
    integer(kind=8) :: npgn, idf1n
    integer(kind=8) :: ivf2, ino, i, nddl1
    integer(kind=8) :: ivarim, ivarip, jvInstmr, jvInstpr
    integer(kind=8) :: iddlm, iddld, jvCarcri
    integer(kind=8) :: ivectu, icontp
    integer(kind=8) :: ivarix
    integer(kind=8) :: jtab(7), jcret, codret
    integer(kind=8) :: ndim, iret, ntrou
    integer(kind=8) :: iu(3, 3), iuc(3), im(3), isect, icontm
    real(kind=8) :: tang(3, 3), a, geom(3, 3)
    character(len=16) :: defoComp
    aster_logical :: lVect, lMatr, lVari, lSigm
    blas_int :: b_incx, b_incy, b_n
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    imatuu = 1
    ivectu = 1
    icontp = 1
    ivarip = 1

    ! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))

    ! - FONCTIONS DE FORME
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', jvf=ivf1, jdfde=idf1)
    call elrefe_info(elrefe=lielrf(1), fami='NOEU', nno=nno1, npg=npgn, jdfde=idf1n)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', nno=nno2, npg=npg, jpoids=iw, &
                     jvf=ivf2)
    ndim = 3
    nddl1 = 5

! - DECALAGE D'INDICE POUR LES ELEMENTS D'INTERFACE
    call cginit(nomte, iu, iuc, im)

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', iddlm)
    call jevech('PDEPLPR', 'L', iddld)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PCAGNBA', 'L', isect)
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Get fields for non-linear behaviour
    call jevech('PCARCRI', 'L', jvCarcri)
    call jevech('PCOMPOR', 'L', vk16=compor)

! - Properties of behaviour
    defoComp = compor(DEFO)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Create compor kit
    comporKit(1:COMPOR_SIZE) = compor(1:COMPOR_SIZE)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              comporKit, zr(jvCarcri), &
                              zr(jvInstmr), zr(jvInstpr), &
                              materPara, BEHInteg)

!   MISE A JOUR EVENTUELLE DE LA GEOMETRIE
    if (defoComp .eq. 'PETIT') then
        do ino = 1, nno1
            do i = 1, ndim
                geom(i, ino) = zr(jvGeom-1+(ino-1)*ndim+i)
            end do
        end do

    else if (defoComp .eq. 'PETIT_REAC') then
        do ino = 1, nno1
            do i = 1, ndim
                geom(i, ino) = zr(jvGeom-1+(ino-1)*ndim+i)+ &
                               zr(iddlm-1+(ino-1)*nddl1+i)+ &
                               zr(iddld-1+(ino-1)*nddl1+i)
            end do
        end do

    else
        call utmess('F', 'CABLE0_6', sk=defoComp)

    end if

!   DEFINITION DES TANGENTES
    call cgtang(3, nno1, npgn, geom, zr(idf1n), tang)

!   SECTION DE LA BARRE
    a = zr(isect)

! - ON VERIFIE QUE PVARIMR ET PVARIPR ONT LE MEME NOMBRE DE V.I. :
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg1 = max(jtab(6), 1)*jtab(7)
!
    if (lVari) then
        call tecach('OOO', 'PVARIPR', 'E', iret, nval=7, itab=jtab)
        lgpg2 = max(jtab(6), 1)*jtab(7)
        ASSERT(lgpg1 .eq. lgpg2)
    end if
    lgpg = lgpg1

! - Get output fields
    if (lMatr) then
        call jevech('PMATUNS', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if

! - FORCES INTERIEURES ET MATRICE TANGENTE
    call cgfint(BEHInteg, &
                typmod, option, &
                comporKit, zr(jvCarcri), &
                ndim, nno1, nno2, npg, &
                zr(iw), zr(ivf1), zr(ivf2), zr(idf1), &
                geom, tang, &
                zr(jvInstmr), zr(jvInstpr), zr(iddlm), zr(iddld), &
                iu, iuc, im, a, &
                lgpg, zr(icontm), zr(ivarim), &
                zr(icontp), zr(ivarip), &
                zr(imatuu), zr(ivectu), &
                codret)
!
    if (lSigm) then
        zi(jcret) = codret
    end if
!
end subroutine
