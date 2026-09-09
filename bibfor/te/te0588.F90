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
subroutine te0588(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterc/ismaem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elref1.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/xasshm.h"
#include "asterfort/xcaehm.h"
#include "asterfort/xfnohm.h"
#include "asterfort/xhmddl.h"
#include "asterfort/xhmini.h"
#include "asterfort/xpeshm.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_HM_*, D_PLAN_HM_* for XFEM
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          FORC_NODA, CHAR_MECA_PESA_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno, imatuu, ndim, jvMaterc, jvInstmr, jcret
    integer(kind=8) :: dimmat, npi, npg, li, ibid, yaenrm
    integer(kind=8) :: codret, icodre(1)
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom
    integer(kind=8) :: jvInstpr, ideplm, ideplp, jvCarcri, ipesa
    integer(kind=8) :: icontm, ivarip, ivarim, ivectu, icontp
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), dimuel
    integer(kind=8) :: dimdef, dimcon, nbvari, nddls, nddlm
    integer(kind=8) :: nmec, np1, np2, nnos
    integer(kind=8) :: nnom
    real(kind=8) :: defgep(13), defgem(13)
    real(kind=8) :: dfdi(20, 3), dfdi2(20, 3)
    real(kind=8) :: drds(25, 11+5), drdsr(25, 11+5), dsde(11+5, 25)
    real(kind=8) :: r(25), sigbar(25), c(25), ck(25), cs(25)
    real(kind=8) :: angnau(3)
    real(kind=8), dimension(:, :), pointer :: work1 => null(), work2 => null(), b => null()
    character(len=3) :: modint
    character(len=8) :: typmod(2)
    character(len=16) :: thmDiffusionKeyword, elref
    real(kind=8) :: rho(1), rbid(1)
    aster_logical :: axi
    type(THM_DS) :: ds_thm
! =====================================================================
!  CETTE ROUTINE FAIT UN CALCUL EN HM AVEC XFEM
!  25 = (9 DEF MECA) + (9 DEF HEAV MECA) + 4 POUR P1 + 3 pour P1 HEAV
!  16 = 12 MECA + 4 POUR P1
! =====================================================================
!  POUR LES TABLEAUX DEFGEP ET DEFGEM ON A DANS L'ORDRE :
!                                      (PARTIE CLASSIQUE)
!                                      DX DY DZ
!                                      EPXX EPYY EPZZ EPXY EPXZ EPYZ
!                                      PRE1 P1DX P1DY P1DZ
!                                      (PARTIE ENRICHIE)
!                                      H1X  H1Y H1Z  H2X  H2Y  H2Z
!                                      H3X  H3Y H3Z
!                                      H1PRE1  H2PRE1  H3PRE1
!            EPSXY = RAC2/2*(DU/DY+DV/DX)
! =====================================================================
!    POUR LES CHAMPS DE CONTRAINTE
!                                      SIXX SIYY SIZZ SIXY SIXZ SIYZ
!                                      SIPXX SIPYY SIPZZ SIPXY SIPXZ SIPYZ
!                                      M11 FH11X FH11Y FH11Z
!
!        SIXY EST LE VRAI DE LA MECANIQUE DES MILIEUX CONTINUS
!        DANS EQUTHM ON LE MULITPLIERA PAR RAC2
! =====================================================================
!   POUR L'OPTION FORCNODA
!  SI LES TEMPS PLUS ET MOINS SONT PRESENTS
!  C'EST QUE L'ON APPELLE DEPUIS STAT NON LINE  : FNOEVO = VRAI
!  ET ALORS LES TERMES DEPENDANT DE DT SONT EVALUES
!  SI LES TEMPS PLUS ET MOINS NE SONT PAS PRESENTS
!  C'EST QUE L'ON APPELLE DEPUIS CALCNO  : FNOEVO = FAUX
!  ET ALORS LES TERMES DEPENDANT DE DT NE SONT PAS EVALUES
! =====================================================================
! AXI       AXISYMETRIQUE?
! TYPMOD    MODELISATION (D_PLAN, AXI, 3D ?)
! MODINT    METHODE D'INTEGRATION (CLASSIQUE,LUMPEE(D),REDUITE(R) ?)
! NNO       NB DE NOEUDS DE L'ELEMENT
! NNOS      NB DE NOEUDS SOMMETS DE L'ELEMENT
! NNOM      NB DE NOEUDS MILIEUX DE L'ELEMENT
! NDDLS     NB DE DDL SUR LES SOMMETS
! NDDLM     NB DE DDL SUR LES MILIEUX
! NPI       NB DE POINTS D'INTEGRATION DE L'ELEMENT
! NPG       NB DE POINTS DE GAUSS     POUR CLASSIQUE(=NPI)
!                 SOMMETS             POUR LUMPEE   (=NPI=NNOS)
!                 POINTS DE GAUSS     POUR REDUITE  (<NPI)
! NDIM      DIMENSION DE L'ESPACE
! DIMUEL    NB DE DDL TOTAL DE L'ELEMENT
! DIMCON    DIMENSION DES CONTRAINTES GENERALISEES ELEMENTAIRES
! DIMDEF    DIMENSION DES DEFORMATIONS GENERALISEES ELEMENTAIRES
! IVF       FONCTIONS DE FORMES QUADRATIQUES
! =====================================================================
    character(len=8), parameter :: famiMater = 'RIGI'
    integer(kind=8) :: nfh, nfiss, jfisno, ddlc, contac
    integer(kind=8) :: ddld, ddlm, ddlp, nnop, nnops, nnopm
    integer(kind=8) :: enrmec(3), nenr, dimenr, enrhyd(3)
    integer(kind=8) :: jpintt, jcnset, jheavt, jpmilt, jheavn
    integer(kind=8) :: jlonch, jlst, jstno
    character(len=8) :: enr
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=16) :: comporCopy(COMPOR_SIZE)
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: itabin(2), iSigm, iret
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    angnau = 0.d0
    imatuu = ismaem()
    ivectu = ismaem()
    icontp = ismaem()
    ivarip = ismaem()
    allocate (work1(11+5, 52*20))
    allocate (work2(25, 52*20))
    allocate (b(25, 52*20))

! - Get model of finite element
    call thmGetElemModel(ds_thm)
!
    call xhmini(nomte, nfh, ddld, ddlm, ddlp, nfiss, ddlc, contac)
    call xcaehm(ds_thm, nomte, axi, typmod, modint, &
                mecani, press1, press2, tempe, dimdef, &
                dimcon, nmec, np1, np2, ndim, &
                nno, nnos, nnom, npi, npg, &
                nddls, nddlm, dimuel, ipoids, ivf, &
                idfde, ddld, ddlm, ddlp, enrmec, nenr, &
                dimenr, nnop, nnops, nnopm, enrhyd, ddlc, nfh)

! - Get finite element
    call elref1(elref)
    call teattr('S', 'XFEM', enr, ibid)
    ASSERT(enr(1:2) .eq. 'XH')

! - Input fields for XFEM
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PLST', 'L', jlst)
    call jevech('PSTANO', 'L', jstno)
    call jevech('PHEA_NO', 'L', jheavn)
    if ((ibid .eq. 0) .and. (enr(1:2) .eq. 'XH') .and. .not. iselli(elref)) then
        call jevech('PPMILTO', 'L', jpmilt)
    end if
    if (nfiss .gt. 1) then
        call jevech('PFISNO', 'L', jfisno)
    end if

! - Generic input fields
    call jevech('PGEOMER', 'L', jvGeom)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(famiMater, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Transfer type of elasticity
    ds_thm%ds_material%elas%id = materPara%elasID
    ds_thm%ds_material%elas%keyword = materPara%elasKeyword

    if ((option(1:9) .eq. 'RIGI_MECA') .or. (option(1:9) .eq. 'RAPH_MECA') .or. &
        (option(1:9) .eq. 'FULL_MECA')) then
        call jevech('PINSTMR', 'L', jvInstmr)
        call jevech('PINSTPR', 'L', jvInstpr)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PCONTMR', 'L', icontm)

! ----- Get fields for non-linear behaviour
        call jevech('PCARCRI', 'L', jvCarcri)
        call jevech('PCOMPOR', 'L', vk16=compor)

! ----- Force DEFO_LDC="MECANIQUE" for THM
        comporCopy = compor
        if (option(1:9) .eq. 'RIGI_MECA') then
            comporCopy(DEFO_LDC) = "MECANIQUE"
        end if
        read (comporCopy(NVAR), '(I16)') nbvari

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  zr(jvInstmr), zr(jvInstpr), &
                                  materPara, ds_thm%ds_behaviour%BEHInteg)

! ----- Select objects to construct from option name
        call behaviourOption(option, comporCopy, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

! ----- Output fields
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
        if (lVari) then
            call jevech('PVARIPR', 'E', ivarip)
        end if
        if (lSigm) then
            call jevech('PCONTPR', 'E', icontp)
            call tecach('OOO', 'PCONTMR', 'L', iret=iret, nval=2, itab=itabin)
            do iSigm = 1, itabin(2)
                zr(icontp-1+iSigm) = zr(icontm-1+iSigm)
            end do
            call jevech('PCODRET', 'E', jcret)
        end if
! ----- Compute
        codret = 0
        dimmat = nddls*nnop
        if (option(1:9) .eq. 'RIGI_MECA') then
            call xasshm(ds_thm, &
                        nno, npg, npi, ipoids, ivf, &
                        idfde, jvGeom, zr(jvGeom), zr(jvCarcri), zr(ideplm), &
                        zr(ideplm), zr(icontm), zr(icontp), zr(ivarim), zr(ivarim), &
                        defgem, defgep, drds, drdsr, dsde, &
                        b, dfdi, dfdi2, r, sigbar, &
                        c, ck, cs, zr(imatuu), zr(ivectu), &
                        zr(jvInstmr), zr(jvInstpr), option, mecani, &
                        press1, press2, tempe, dimdef, dimcon, &
                        dimuel, nbvari, nddls, nddlm, nmec, &
                        np1, ndim, comporCopy, axi, modint, &
                        codret, nnop, nnops, nnopm, enrmec, &
                        dimenr, zi(jheavt), zi(jlonch), zi(jcnset), jpintt, &
                        jpmilt, jheavn, dimmat, enrhyd, nfiss, nfh, jfisno, &
                        work1, work2, lVect, lMatr, lVari, lSigm)
        else
            do li = 1, dimuel
                zr(ideplp+li-1) = zr(ideplm+li-1)+zr(ideplp+li-1)
            end do
            call xasshm(ds_thm, &
                        nno, npg, npi, ipoids, ivf, &
                        idfde, jvGeom, zr(jvGeom), zr(jvCarcri), zr(ideplm), &
                        zr(ideplp), zr(icontm), zr(icontp), zr(ivarim), zr(ivarip), &
                        defgem, defgep, drds, drdsr, dsde, &
                        b, dfdi, dfdi2, r, sigbar, &
                        c, ck, cs, zr(imatuu), zr(ivectu), &
                        zr(jvInstmr), zr(jvInstpr), option, mecani, &
                        press1, press2, tempe, dimdef, dimcon, &
                        dimuel, nbvari, nddls, nddlm, nmec, &
                        np1, ndim, comporCopy, axi, modint, &
                        codret, nnop, nnops, nnopm, enrmec, &
                        dimenr, zi(jheavt), zi(jlonch), zi(jcnset), jpintt, &
                        jpmilt, jheavn, dimmat, enrhyd, nfiss, nfh, jfisno, &
                        work1, work2, lVect, lMatr, lVari, lSigm)
        end if
        if (lSigm) then
            zi(jcret) = codret
        end if

        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, zr(imatuu), &
                    zr(ivectu), nddlm, nfiss, jfisno, .false._1, contac)
    end if

    if (option .eq. 'CHAR_MECA_PESA_R') then
        call jevech('PPESANR', 'L', ipesa)
        call jevech('PVECTUR', 'E', ivectu)

        call rccoma(zi(jvMaterc), 'THM_DIFFU', 1, thmDiffusionKeyword, icodre(1))
        call rcvalb('FPG1', 1, 1, '+', zi(jvMaterc), &
                    ' ', thmDiffusionKeyword, 0, ' ', [0.d0], &
                    1, 'RHO', rho(1), icodre, 1)

!       INDICATEUR POUR SAVOIR SI ON A DE L'ENRICHISSEMENT
        yaenrm = enrmec(1)
!
        call xpeshm(nno, nnop, nnops, ndim, nddls, &
                    nddlm, npg, jvGeom, jpintt, jpmilt, jheavn, &
                    ivf, ipoids, idfde, ivectu, ipesa, &
                    zi(jheavt), zi(jlonch), zi(jcnset), rho(1), axi, &
                    yaenrm, nfiss, nfh, jfisno)
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, rbid, &
                    zr(ivectu), nddlm, nfiss, jfisno, .false._1, contac)

    end if

    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', icontm)
        call jevech('PVECTUR', 'E', ivectu)
        ds_thm%ds_behaviour%BEHInteg%materPara = materPara
        call xfnohm(ds_thm, &
                    nno, npg, ipoids, &
                    ivf, idfde, zr(jvGeom), zr(icontm), b, &
                    dfdi, dfdi2, r, zr(ivectu), &
                    mecani, press1, dimcon, nddls, nddlm, &
                    dimuel, nmec, np1, ndim, axi, &
                    dimenr, nnop, nnops, nnopm, jvGeom, &
                    jpintt, jpmilt, jheavn, zi(jlonch), zi(jcnset), zi(jheavt), &
                    enrmec, enrhyd, nfiss, nfh, jfisno)
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, rbid, &
                    zr(ivectu), nddlm, nfiss, jfisno, .false._1, contac)
    end if
!
    deallocate (work1)
    deallocate (work2)
    deallocate (b)
end subroutine
