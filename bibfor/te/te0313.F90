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
subroutine te0313(option, nomte)
!
    use THM_type
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/ismaem.h"
#include "asterf_types.h"
#include "asterfort/aseihm.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/caeihm.h"
#include "asterfort/fneihm.h"
#include "asterfort/jevech.h"
#include "asterfort/poeihm.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/thmGetElemModel.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_JOINT_HYME
!           PLAN_JOINT_HYME
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          FORC_NODA, VARI_ELNO, SIEF_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: jgano, imatuu, ndim, jvMaterc, jvInstmr, jcret, nb_strain_meca
    integer(kind=8) :: iret, ichg, ichn, itabin(7), itabou(7), jvCarcri
    integer(kind=8) :: ivf2
    integer(kind=8) :: idf2, npi, npg
    integer(kind=8) :: codret
    integer(kind=8) :: ipoids, ivf1, idf1, jvGeom
    integer(kind=8) :: jvInstpr, ideplm, ideplp
    integer(kind=8) :: icontm, ivarip, ivarim, ivectu, icontp, jvSief
    integer(kind=8) :: mecani(8), press1(9), press2(9), tempe(5), dimuel
    integer(kind=8) :: dimdef, dimcon, nbvari, nb_vari_meca
    integer(kind=8) :: nno1, nno2
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: iu(3, 18), ip(2, 9), ipf(2, 2, 9), iq(2, 2, 9)
    real(kind=8) :: r(22)
    character(len=3) :: intgType
    character(len=8) :: typmod(2)
    type(THM_DS) :: ds_thm
    integer(kind=8) :: li
    aster_logical :: axi
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: comporCopy(COMPOR_SIZE)
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!

! - Modelisation
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))

! - Get model of finite element
    call thmGetElemModel(ds_thm)

! - Preparation
    call caeihm(ds_thm, nomte, axi, mecani, press1, &
                press2, tempe, dimdef, dimcon, ndim, &
                nno1, nno2, npi, npg, dimuel, &
                ipoids, ivf1, idf1, ivf2, idf2, &
                jgano, iu, ip, ipf, iq, &
                intgType)

! - Get material parameters
    if ((option .eq. 'FORC_NODA') .or. (option(1:9) .eq. 'RIGI_MECA') .or. &
        (option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA')) then
        call jevech('PMATERC', 'L', jvMaterc)
        call initParaCell(fami, zi(jvMaterc), materPara)
        call initLCSPg(ndim, nno2, materPara)
    end if

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

    if ((option(1:9) .eq. 'RIGI_MECA') .or. (option(1:9) .eq. 'RAPH_MECA') .or. &
        (option(1:9) .eq. 'FULL_MECA')) then

! ----- Input fields
        call jevech('PGEOMER', 'L', jvGeom)
        call jevech('PINSTMR', 'L', jvInstmr)
        call jevech('PINSTPR', 'L', jvInstpr)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PCONTMR', 'L', icontm)

! ----- Get fields for non-linear behaviour
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', jvCarcri)

! ----- Force DEFO_LDC="MECANIQUE" for THM
        comporCopy(1:COMPOR_SIZE) = compor(1:COMPOR_SIZE)
        comporCopy(DEFO_LDC) = "MECANIQUE"

! ----- Properties of behaviour
        read (comporCopy(NVAR), '(I16)') nbvari

! ----- Select objects to construct from option name
        call behaviourOption(option, comporCopy, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  zr(jvInstmr), zr(jvInstpr), &
                                  materPara, BEHInteg)
        ds_thm%ds_behaviour%BEHInteg = BEHInteg

! ----- Output fields
        imatuu = ismaem()
        ivectu = ismaem()
        icontp = ismaem()
        ivarip = ismaem()
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
        end if
        if (lVari) then
            call jevech('PVARIPR', 'E', ivarip)
        end if
        if (lSigm) then
            call jevech('PCONTPR', 'E', icontp)
            call jevech('PCODRET', 'E', jcret)
            zi(jcret) = 0
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
! ----- Integration
        codret = 0
        if (option(1:9) .eq. 'RIGI_MECA') then
            call aseihm(ds_thm, option, &
                        lSigm, lVari, lMatr, lVect, &
                        axi, ndim, nno1, nno2, &
                        npi, npg, dimuel, dimdef, dimcon, &
                        nbvari, iu, ip, ipf, &
                        iq, mecani, press1, press2, tempe, &
                        zr(ivf1), zr(ivf2), zr(idf2), zr(jvInstmr), zr(jvInstpr), &
                        zr(ideplm), zr(ideplm), zr(icontm), zr(icontp), zr(ivarim), &
                        zr(ivarim), zr(ipoids), zr(jvGeom), &
                        comporCopy, zr(ivectu), zr(imatuu), &
                        codret)
        else
            do li = 1, dimuel
                zr(ideplp+li-1) = zr(ideplm+li-1)+zr(ideplp+li-1)
            end do
            call aseihm(ds_thm, option, &
                        lSigm, lVari, lMatr, lVect, &
                        axi, ndim, nno1, nno2, &
                        npi, npg, dimuel, dimdef, dimcon, &
                        nbvari, iu, ip, ipf, &
                        iq, mecani, press1, press2, tempe, &
                        zr(ivf1), zr(ivf2), zr(idf2), zr(jvInstmr), zr(jvInstpr), &
                        zr(ideplm), zr(ideplp), zr(icontm), zr(icontp), zr(ivarim), &
                        zr(ivarip), zr(ipoids), zr(jvGeom), &
                        comporCopy, zr(ivectu), zr(imatuu), &
                        codret)
            if (lSigm) then
                zi(jcret) = codret
            end if
        end if
    end if

    if (option .eq. 'FORC_NODA') then
        call jevech('PGEOMER', 'L', jvGeom)
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PVECTUR', 'E', ivectu)
        BEHInteg%materPara = materPara
        ds_thm%ds_behaviour%BEHInteg = BEHInteg
        call fneihm(ds_thm, &
                    nno1, nno2, &
                    npi, npg, zr(ipoids), iu, ip, &
                    ipf, iq, zr(ivf1), zr(ivf2), zr(idf2), &
                    zr(jvGeom), zr(jvSief), r, zr(ivectu), &
                    mecani, press1, press2, dimdef, &
                    dimcon, dimuel, ndim, axi)
!
    end if

    if (option .eq. 'SIEF_ELNO') then
        call jevech('PCONTRR', 'L', ichg)
        call jevech('PSIEFNOR', 'E', ichn)
        nb_strain_meca = mecani(6)
        call poeihm(nomte, option, intgType, jgano, nno1, &
                    nno2, dimcon, nb_strain_meca, zr(ichg), zr(ichn))
    end if

    if (option .eq. 'VARI_ELNO') then
        call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, itab=itabin)
        call tecach('OOO', 'PVARINR', 'E', iret, nval=7, itab=itabou)
        ichg = itabin(1)
        ichn = itabou(1)
!
        call jevech('PCOMPOR', 'L', vk16=compor)
        read (compor(NVAR), '(I16)') nbvari
        read (compor(MECA_NVAR), '(I16)') nb_vari_meca
        call poeihm(nomte, option, intgType, jgano, nno1, &
                    nno2, nbvari, nb_vari_meca, zr(ichg), zr(ichn))
    end if
!
end subroutine
