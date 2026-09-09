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
! aslint: disable=W1504
!
subroutine get_elas_para(fami, jvMaterCode, poum, ipg, ispg, &
                         elas_id, elas_keyword, &
                         time, temp, &
                         e_, nu_, g_, &
                         e1_, e2_, e3_, &
                         nu12_, nu13_, nu23_, &
                         g1_, g2_, g3_, &
                         BEHInteg, &
                         ei_, nui_, gi_, &
                         e1i_, e2i_, e3i_, &
                         nu12i_, nu13i_, nu23i_, &
                         g1i_, g2i_, g3i_)
!
    use Behaviour_type
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/hypmat.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvalc.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: ipg, ispg
    integer(kind=8), intent(in) :: elas_id
    character(len=16), intent(in) :: elas_keyword
    real(kind=8), optional, intent(in) :: time
    real(kind=8), optional, intent(in) :: temp
    real(kind=8), optional, intent(out) :: e_, nu_, g_
    real(kind=8), optional, intent(out) :: ei_, nui_, gi_
    real(kind=8), optional, intent(out) :: e1_, e2_, e3_
    real(kind=8), optional, intent(out) :: e1i_, e2i_, e3i_
    real(kind=8), optional, intent(out) :: nu12_, nu13_, nu23_
    real(kind=8), optional, intent(out) :: nu12i_, nu13i_, nu23i_
    real(kind=8), optional, intent(out) :: g1_, g2_, g3_
    real(kind=8), optional, intent(out) :: g1i_, g2i_, g3i_
    type(Behaviour_Integ), optional, intent(in) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility
!
! Get elastic parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  jvMaterCode      : coded material address
! In  time             : current time
! In  time             : current temperature
! In  poum             : '-' or '+' for parameters evaluation (previous or current temperature)
! In  ipg              : current point gauss
! In  ispg             : current "sous-point" gauss
! In  elas_id          : Type of elasticity
!                 1 - Isotropic
!                 2 - Orthotropic
!                 3 - Transverse isotropic
!                    or viscoelasticity
!                 4 - Isotropic
!                 5 - Orthotropic
!                 6 - Transverse isotropic
! In  elas_keyword     : keyword factor linked to type of elasticity parameters
! Out e                : real Young modulus (isotropic)
! Out ei               : imaginary Young modulus (isotropic)
! Out nu               : real Poisson ratio (isotropic)
! Out nui              : imaginary Poisson ratio (isotropic)
! Out e1               : real Young modulus - Direction 1 (Orthotropic/Transverse isotropic)
! Out e2               : real Young modulus - Direction 2 (Orthotropic)
! Out e3               : real Young modulus - Direction 3 (Orthotropic/Transverse isotropic)
! Out e1i              : imaginary Young modulus - Direction 1 (Orthotropic/Transverse isotropic)
! Out e2i              : imaginary Young modulus - Direction 2 (Orthotropic)
! Out e3i              : imaginary Young modulus - Direction 3 (Orthotropic/Transverse isotropic)
! Out nu12             : real Poisson ratio - Coupling 1/2 (Orthotropic/Transverse isotropic)
! Out nu13             : real Poisson ratio - Coupling 1/3 (Orthotropic/Transverse isotropic)
! Out nu23             : real Poisson ratio - Coupling 2/3 (Orthotropic)
! Out nu12i            : imaginary Poisson ratio - Coupling 1/2 (Orthotropic/Transverse isotropic)
! Out nu13i            : imaginary Poisson ratio - Coupling 1/3 (Orthotropic/Transverse isotropic)
! Out nu23i            : imaginary Poisson ratio - Coupling 2/3 (Orthotropic)
! Out g1               : real shear ratio (Orthotropic)
! Out g2               : real shear ratio (Orthotropic)
! Out g3               : real shear ratio (Orthotropic)
! Out g1i              : imaginary shear ratio (Orthotropic)
! Out g2i              : imaginary shear ratio (Orthotropic)
! Out g3i              : imaginary shear ratio (Orthotropic)
! Out g                : real shear ratio (isotropic/Transverse isotropic)
! Out gi               : imaginary shear ratio (isotropic)
! In  BEHInteg         : parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbresm = 9
    integer(kind=8) :: icodre(nbresm)
    character(len=16) :: nomres(nbresm)
    real(kind=8) :: valres(nbresm)
    complex(kind=8) :: valresc(nbresm)
!
    character(len=8) :: para_name(5)
    real(kind=8) :: para_vale(5)
    integer(kind=8) :: nbres, nb_para
    real(kind=8), parameter :: un = 1.d0
    real(kind=8) :: c10, c01, c20, k
    real(kind=8) :: er, nur, gr
    real(kind=8) :: ei, nui, gi
    real(kind=8) :: e1r, e2r, e3r
    real(kind=8) :: e1i, e2i, e3i
    real(kind=8) :: nu12r, nu13r, nu23r
    real(kind=8) :: nu12i, nu13i, nu23i
    real(kind=8) :: g1r, g2r, g3r
    real(kind=8) :: g1i, g2i, g3i
    complex(kind=8) :: Gc, nuc
!
! --------------------------------------------------------------------------------------------------
!
    er = 0.d0
    nur = 0.d0
    gr = 0.d0
    ei = 0.d0
    nui = 0.d0
    gi = 0.d0
    e1r = 0.d0
    e2r = 0.d0
    e3r = 0.d0
    e1i = 0.d0
    e2i = 0.d0
    e3i = 0.d0
    g1r = 0.d0
    g2r = 0.d0
    g3r = 0.d0
    g1i = 0.d0
    g2i = 0.d0
    g3i = 0.d0
!
    nb_para = 0
    para_name = ' '
    para_vale = 0.d0
    if (present(time)) then
        if (time .ne. r8vide()) then
            nb_para = nb_para+1
            para_name(nb_para) = 'INST'
            para_vale(nb_para) = time
        end if
    end if
    if (present(temp)) then
        nb_para = nb_para+1
        para_name(nb_para) = 'TEMP'
        para_vale(nb_para) = temp
    end if
    if (present(BEHInteg)) then
        if (.not. BEHInteg%behavESVA%lGeomInESVA .and. (fami .ne. "XFEM")) then
            ASSERT(ipg <= VARC_GEOM_NBMAXI)
            nb_para = nb_para+1
            para_name(nb_para) = 'X'
            para_vale(nb_para) = BEHInteg%behavESVA%behavESVAGeom%coorElga(ipg, 1)
            nb_para = nb_para+1
            para_name(nb_para) = 'Y'
            para_vale(nb_para) = BEHInteg%behavESVA%behavESVAGeom%coorElga(ipg, 2)
            nb_para = nb_para+1
            para_name(nb_para) = 'Z'
            para_vale(nb_para) = BEHInteg%behavESVA%behavESVAGeom%coorElga(ipg, 3)
        end if
    end if
!
! - Get elastic parameters
!
    if (elas_id .eq. 1) then
        if (elas_keyword .eq. 'ELAS_HYPER' .or. elas_keyword .eq. 'ELAS_HYPER_VISC') then
            call hypmat(fami, ipg, ispg, poum, jvMaterCode, &
                        c10, c01, c20, k)
            nur = (3.d0*k-4.0d0*(c10+c01))/(6.d0*k+4.0d0*(c10+c01))
            er = 4.d0*(c10+c01)*(un+nur)
            gr = er/(2.d0*(1.d0+nur))
        else
            nomres(1) = 'E'
            nomres(2) = 'NU'
            nbres = 2
            call rcvalb(fami, ipg, ispg, poum, jvMaterCode, &
                        ' ', elas_keyword, nb_para, para_name, [para_vale], &
                        nbres, nomres, valres, icodre, 1)
            er = valres(1)
            nur = valres(2)
            gr = er/(2.d0*(1.d0+nur))
        end if

    elseif (elas_id .eq. 2) then
        nomres(1) = 'E_L'
        nomres(2) = 'E_T'
        nomres(3) = 'E_N'
        nomres(4) = 'NU_LT'
        nomres(5) = 'NU_LN'
        nomres(6) = 'NU_TN'
        nomres(7) = 'G_LT'
        nomres(8) = 'G_LN'
        nomres(9) = 'G_TN'
        nbres = 9
        call rcvalb(fami, ipg, ispg, poum, jvMaterCode, &
                    ' ', elas_keyword, nb_para, para_name, [para_vale], &
                    nbres, nomres, valres, icodre, 1)
        e1r = valres(1)
        e2r = valres(2)
        e3r = valres(3)
        nu12r = valres(4)
        nu13r = valres(5)
        nu23r = valres(6)
        g1r = valres(7)
        g2r = valres(8)
        g3r = valres(9)
    elseif (elas_id .eq. 3) then
        nomres(1) = 'E_L'
        nomres(2) = 'E_N'
        nomres(3) = 'NU_LT'
        nomres(4) = 'NU_LN'
        nomres(5) = 'G_LN'
        nbres = 5
        call rcvalb(fami, ipg, ispg, poum, jvMaterCode, &
                    ' ', elas_keyword, nb_para, para_name, [para_vale], &
                    nbres, nomres, valres, icodre, 1)
        e1r = valres(1)
        e3r = valres(2)
        nu12r = valres(3)
        nu13r = valres(4)
        gr = valres(5)
    elseif (elas_id .eq. 4) then
        nomres(1) = 'G'
        nomres(2) = 'NU'
        nbres = 2
        call rcvalc(jvMaterCode, elas_keyword, nbres, nomres, valresc, icodre, 1)
        Gc = valresc(1)
        nuc = valresc(2)
        nur = real(nuc)
        nui = aimag(nuc)
        er = 2.d0*real(Gc*(1.d0+nuc))
        ei = 2.d0*aimag(Gc*(1.d0+nuc))
        gr = real(Gc)
        gi = aimag(Gc)
    elseif (elas_id .eq. 5) then
        nomres(1) = 'E_L'
        nomres(2) = 'E_T'
        nomres(3) = 'E_N'
        nomres(4) = 'NU_LT'
        nomres(5) = 'NU_LN'
        nomres(6) = 'NU_TN'
        nomres(7) = 'G_LT'
        nomres(8) = 'G_LN'
        nomres(9) = 'G_TN'
        nbres = 9
        call rcvalc(jvMaterCode, elas_keyword, nbres, nomres, valresc, icodre, 1)
        e1r = real(valresc(1))
        e2r = real(valresc(2))
        e3r = real(valresc(3))
        nu12r = real(valresc(4))
        nu13r = real(valresc(5))
        nu23r = real(valresc(6))
        g1r = real(valresc(7))
        g2r = real(valresc(8))
        g3r = real(valresc(9))
        e1i = aimag(valresc(1))
        e2i = aimag(valresc(2))
        e3i = aimag(valresc(3))
        nu12i = aimag(valresc(4))
        nu13i = aimag(valresc(5))
        nu23i = aimag(valresc(6))
        g1i = aimag(valresc(7))
        g2i = aimag(valresc(8))
        g3i = aimag(valresc(9))
    elseif (elas_id .eq. 6) then
        nomres(1) = 'E_L'
        nomres(2) = 'E_N'
        nomres(3) = 'NU_LT'
        nomres(4) = 'NU_LN'
        nomres(5) = 'G_LN'
        nbres = 5
        call rcvalc(jvMaterCode, elas_keyword, nbres, nomres, valresc, icodre, 1)
        e1r = real(valresc(1))
        e3r = real(valresc(2))
        nu12r = real(valresc(3))
        nu13r = real(valresc(4))
        gr = real(valresc(5))
        e1i = aimag(valresc(1))
        e3i = aimag(valresc(2))
        nu12i = aimag(valresc(3))
        nu13i = aimag(valresc(4))
        gi = aimag(valresc(5))
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Output
!
    if (present(e_)) then
        e_ = er
    end if
    if (present(nu_)) then
        nu_ = nur
    end if
    if (present(g_)) then
        g_ = gr
    end if
    if (present(ei_)) then
        ei_ = ei
    end if
    if (present(nui_)) then
        nui_ = nui
    end if
    if (present(gi_)) then
        gi_ = gi
    end if
    if (present(e1_)) then
        e1_ = e1r
    end if
    if (present(e2_)) then
        e2_ = e2r
    end if
    if (present(e3_)) then
        e3_ = e3r
    end if
    if (present(g1_)) then
        g1_ = g1r
    end if
    if (present(g2_)) then
        g2_ = g2r
    end if
    if (present(g3_)) then
        g3_ = g3r
    end if
    if (present(nu12_)) then
        nu12_ = nu12r
    end if
    if (present(nu13_)) then
        nu13_ = nu13r
    end if
    if (present(nu23_)) then
        nu23_ = nu23r
    end if
    if (present(e1i_)) then
        e1i_ = e1i
    end if
    if (present(e2i_)) then
        e2i_ = e2i
    end if
    if (present(e3i_)) then
        e3i_ = e3i
    end if
    if (present(g1i_)) then
        g1i_ = g1i
    end if
    if (present(g2i_)) then
        g2i_ = g2i
    end if
    if (present(g3i_)) then
        g3i_ = g3i
    end if
    if (present(nu12i_)) then
        nu12i_ = nu12i
    end if
    if (present(nu13i_)) then
        nu13i_ = nu13i
    end if
    if (present(nu23i_)) then
        nu23i_ = nu23i
    end if
!
end subroutine
