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

subroutine dtmprep_noli_deci(sd_dtm_, sd_nl_, icomp)
    implicit none
! dtmprep_noli_deci : prepare the calculations for a localized nonlinearity
!                     of type : discrete model with isotropic behavior. This
!                     routine adds a single occurence to sd_nl and increments
!                     NB_NOLI in sd_dtm
!
!             icomp : an integer giving the index of occurence of the
!                     nonlinearity to be treated under the factor kw
!                     COMPORTEMENT of the command DYNA_VIBRA.
!
!
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterc/r8miem.h"
#include "asterfort/angvx.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmsav.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mdchdl.h"
#include "asterfort/nlget.h"
#include "asterfort/nlinivec.h"
#include "asterfort/nlsav.h"
#include "asterfort/nltype.h"
#include "asterfort/posddl.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), intent(in) :: icomp
!
!   -0.2- Local variables
    aster_logical     :: lnoeu2, l_rota
    integer(kind=8) :: i, n1, nbdeci, nbnoli, ii
    integer(kind=8) :: nbmcl, ier, nbma, ino1, jmama
    integer(kind=8) :: ino2, ind1, ind2, nbmode, elem_nume
    integer(kind=8) :: j, neq, mxlevel, nexcit
    integer(kind=8) :: nunoe, nuddl, nbnodes, nbddl, iddl, n2, n3
!
    real(kind=8) :: r8bid, alpha, beta
    real(kind=8) :: axe(3)
    real(kind=8) :: kelas(6), limy(6), kcine(6), puis(6), limu(6)
!
    character(len=8)  :: sd_dtm, sd_nl, mesh, mesh1, mesh2
    character(len=8)  :: nume, nume1, nume2, no1_name, no2_name
    character(len=8)  :: monmot, intk
    character(len=16) :: typnum, typem, limocl(2)
    character(len=16) :: obst_typ, motfac
    character(len=19) :: nomres
    character(len=24) :: nl_title
!
    integer(kind=8), pointer       :: ddlnodes(:) => null()
    real(kind=8), pointer       :: coor_no1(:) => null()
    real(kind=8), pointer       :: coor_no2(:) => null()
    real(kind=8), pointer       :: vale(:) => null()
    real(kind=8), pointer       :: defmod1(:) => null()
    real(kind=8), pointer       :: defmod2(:) => null()
    real(kind=8), pointer       :: ps2del1(:) => null()
    real(kind=8), pointer       :: ps2del2(:) => null()
    real(kind=8), pointer       :: bmodal_v(:) => null()
    real(kind=8), pointer       :: ps1del_v(:) => null()
!
    integer(kind=8), pointer :: elems(:) => null()
!
#define ps1del(m,n) ps1del_v((n-1)*neq+m)
#define bmodal(m,n) bmodal_v((n-1)*neq+m)
!
!   --- 0. Diverse initializations
    call jemarq()
!
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
!
    lnoeu2 = ASTER_FALSE
    !
    motfac = 'COMPORTEMENT'
    call nlget(sd_nl, _MAX_LEVEL, iscal=mxlevel)
    nbnoli = mxlevel+1
    i = mxlevel+1
!
!   --- 1 - Basic information about the mesh and numbering
!
    call dtmget(sd_dtm, _NUM_DDL, kscal=nume)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
    call gettco(nume, typnum)

!
!   --- 1.1 - Case with a simple modal projection (direct calculation)
    if (typnum(1:16) .eq. 'NUME_DDL_SDASTER') then
        call dismoi('NOM_MAILLA', nume, 'NUME_DDL', repk=mesh)
        mesh1 = mesh
        nume1 = nume
        mesh2 = mesh
        nume2 = nume
        call nlsav(sd_nl, _NUMDDL_1, 1, iocc=i, kscal=nume1(1:8))
        call nlsav(sd_nl, _MESH_1, 1, iocc=i, kscal=mesh1)

!   --- 1.2 - Case with double (or triple) projections (sub-structuring case)
!             Not supported for buckling non linearities
    else if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
        call utmess('F', 'ALGORITH5_36')
    else
        ASSERT(.false.)
    end if
!
!
!   --- 2 - Localisation (support nodes) of the buckling non linearity
!
!   --- 2.1 - Definition using elements or element groups (MAILLE/GROUP_MA)

    typem = 'NU_MAILLE'
    nbmcl = 2
    limocl(1) = 'GROUP_MA'
    limocl(2) = 'MAILLE'
    call reliem(' ', mesh1, typem, motfac, icomp, &
                nbmcl, limocl, limocl, sd_nl//'.INDI_SUPP.TEMP', nbma)
!
    if (nbma .ne. 1) call utmess('F', 'ALGORITH13_38')

    call jeveuo(sd_nl//'.INDI_SUPP.TEMP', 'L', vi=elems)
    elem_nume = elems(1)
    call jeveuo(jexnum(mesh1//'.CONNEX', elem_nume), 'L', jmama)
    call jelira(jexnum(mesh1//'.CONNEX', elem_nume), 'LONMAX', nbnodes)
    ASSERT(nbnodes .ge. 1 .and. nbnodes .le. 2)

    no1_name = int_to_char8(zi(jmama))
    call nlsav(sd_nl, _NO1_NAME, 1, iocc=i, kscal=no1_name)
    l_rota = ASTER_TRUE
    call posddl('NUME_DDL', nume1, no1_name, 'DRX', nunoe, &
                nuddl)
    if (nuddl .eq. 0) l_rota = ASTER_FALSE

    if (nbnodes .eq. 2) then
        lnoeu2 = ASTER_TRUE
        no2_name = int_to_char8(zi(jmama+1))
        call nlsav(sd_nl, _NO2_NAME, 1, iocc=i, kscal=no2_name)
        call nlsav(sd_nl, _NUMDDL_2, 1, iocc=i, kscal=nume2(1:8))
        call nlsav(sd_nl, _MESH_2, 1, iocc=i, kscal=mesh2)
        call posddl('NUME_DDL', nume2, no2_name, 'DRX', nunoe, &
                    nuddl)
        if (nuddl .eq. 0) l_rota = ASTER_FALSE
    end if
    call jedetr(sd_nl//'.INDI_SUPP.TEMP')
!
!   --- 3 - Filling up the sd_nl with further information regarding the
!           nonlinearity(ies)
!
    nbddl = 3
    if (l_rota) nbddl = 6
    AS_ALLOCATE(vi=ddlnodes, size=2*nbddl)
!
!   --- Loop over the detected nonlineary in the current COMPORTEMENT
!       occurence
!
    call nlsav(sd_nl, _NL_TYPE, 1, iocc=i, iscal=NL_DIS_ECRO_CINE)
!
!   --- 3.1 - DOF numbering localisation index for the concerned nodes
    call mdchdl(lnoeu2, i, ddlnodes, ier, l_rota)
!
!   --- 3.2 - Coordinates of the nodes
    call jeveuo(mesh1//'.COORDO    .VALE', 'L', vr=vale)
    ino1 = char8_to_int(no1_name)
    ind1 = 1+3*(ino1-1)
    ind2 = ind1+3
    call nlsav(sd_nl, _COOR_NO1, 3, iocc=i, rvect=vale(ind1:ind2))
    if (lnoeu2) then
        if (mesh2 .ne. mesh1) then
            call jeveuo(mesh2//'.COORDO    .VALE', 'L', vr=vale)
        end if
        ino2 = char8_to_int(no2_name)
        ind1 = 1+3*(ino2-1)
        ind2 = ind1+3
        call nlsav(sd_nl, _COOR_NO2, 3, iocc=i, rvect=vale(ind1:ind2))
    end if
!
!   --- 3.3 - Other information are read from the user input
    call codent(i, 'D0', intk)
    nl_title = nltype(NL_DIS_ECRO_CINE)//intk
    call nlsav(sd_nl, _NL_TITLE, 1, iocc=i, kscal=nl_title)

    kelas(:) = 0.d0
    limy(:) = 0.d0
    kcine(:) = 0.D0
    puis(:) = 0.d0
    limu(:) = 0.d0
    ii = 1
    call getvr8(motfac, 'KELA_DX', iocc=icomp, scal=kelas(ii), nbret=n1)
    if (n1 .gt. 0) then
        call getvr8(motfac, 'LIMY_DX', iocc=icomp, scal=limy(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
        call getvr8(motfac, 'KCIN_DX', iocc=icomp, scal=kcine(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
!
        call getvr8(motfac, 'PUIS_DX', iocc=icomp, scal=puis(ii), nbret=n2)
        if (n2 .gt. 0) then
            call getvr8(motfac, 'LIMU_DX', iocc=icomp, scal=limu(ii), nbret=n3)
            ASSERT(n3 .gt. 0)
        else
            puis(ii) = 0.d0
        end if
    end if

    ii = 2
    call getvr8(motfac, 'KELA_DY', iocc=icomp, scal=kelas(ii), nbret=n1)
    if (n1 .gt. 0) then
        call getvr8(motfac, 'LIMY_DY', iocc=icomp, scal=limy(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
        call getvr8(motfac, 'KCIN_DY', iocc=icomp, scal=kcine(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
!
        call getvr8(motfac, 'PUIS_DY', iocc=icomp, scal=puis(ii), nbret=n2)
        if (n2 .gt. 0) then
            call getvr8(motfac, 'LIMU_DY', iocc=icomp, scal=limu(ii), nbret=n3)
            ASSERT(n3 .gt. 0)
        else
            puis(ii) = 0.d0
        end if
    end if

    ii = 3
    call getvr8(motfac, 'KELA_DZ', iocc=icomp, scal=kelas(ii), nbret=n1)
    if (n1 .gt. 0) then
        call getvr8(motfac, 'LIMY_DZ', iocc=icomp, scal=limy(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
        call getvr8(motfac, 'KCIN_DZ', iocc=icomp, scal=kcine(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
!
        call getvr8(motfac, 'PUIS_DZ', iocc=icomp, scal=puis(ii), nbret=n2)
        if (n2 .gt. 0) then
            call getvr8(motfac, 'LIMU_DZ', iocc=icomp, scal=limu(ii), nbret=n3)
            ASSERT(n3 .gt. 0)
        else
            puis(ii) = 0.d0
        end if
    end if

    ii = 4
    call getvr8(motfac, 'KELA_RX', iocc=icomp, scal=kelas(ii), nbret=n1)
    if (n1 .gt. 0) then
        if (.not. l_rota) call utmess('F', 'DISCRETS_46', sk='KELA_RX')
        call getvr8(motfac, 'LIMY_RX', iocc=icomp, scal=limy(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
        call getvr8(motfac, 'KCIN_RX', iocc=icomp, scal=kcine(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
!
        call getvr8(motfac, 'PUIS_RX', iocc=icomp, scal=puis(ii), nbret=n2)
        if (n2 .gt. 0) then
            call getvr8(motfac, 'LIMU_RX', iocc=icomp, scal=limu(ii), nbret=n3)
            ASSERT(n3 .gt. 0)
        else
            puis(ii) = 0.d0
        end if
    end if

    ii = 5
    call getvr8(motfac, 'KELA_RY', iocc=icomp, scal=kelas(ii), nbret=n1)
    if (n1 .gt. 0) then
        if (.not. l_rota) call utmess('F', 'DISCRETS_46', sk='KELA_RY')
        call getvr8(motfac, 'LIMY_RY', iocc=icomp, scal=limy(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
        call getvr8(motfac, 'KCIN_RY', iocc=icomp, scal=kcine(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
!
        call getvr8(motfac, 'PUIS_RY', iocc=icomp, scal=puis(ii), nbret=n2)
        if (n2 .gt. 0) then
            call getvr8(motfac, 'LIMU_RY', iocc=icomp, scal=limu(ii), nbret=n3)
            ASSERT(n3 .gt. 0)
        else
            puis(ii) = 0.d0
        end if
    end if

    ii = 6
    call getvr8(motfac, 'KELA_RZ', iocc=icomp, scal=kelas(ii), nbret=n1)
    if (n1 .gt. 0) then
        if (.not. l_rota) call utmess('F', 'DISCRETS_46', sk='KELA_RZ')
        call getvr8(motfac, 'LIMY_RZ', iocc=icomp, scal=limy(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
        call getvr8(motfac, 'KCIN_RZ', iocc=icomp, scal=kcine(ii), nbret=n2)
        ASSERT(n2 .gt. 0)
!
        call getvr8(motfac, 'PUIS_RZ', iocc=icomp, scal=puis(ii), nbret=n2)
        if (n2 .gt. 0) then
            call getvr8(motfac, 'LIMU_RZ', iocc=icomp, scal=limu(ii), nbret=n3)
            ASSERT(n3 .gt. 0)
        else
            puis(ii) = 0.d0
        end if
    end if

    call nlsav(sd_nl, _ECRCIN_KELA, 6, iocc=i, rvect=kelas)
    call nlsav(sd_nl, _ECRCIN_LIMY, 6, iocc=i, rvect=limy)
    call nlsav(sd_nl, _ECRCIN_KCIN, 6, iocc=i, rvect=kcine)
    call nlsav(sd_nl, _ECRCIN_PUIS, 6, iocc=i, rvect=puis)
    call nlsav(sd_nl, _ECRCIN_LIMU, 6, iocc=i, rvect=limu)
!
!   --- 3.4 - Calculation of geometrical properties :
!             play, orientation, local coordinates, distances
!             Vector x1x2 must not be zero
    if (lnoeu2) then
        call nlget(sd_nl, _COOR_NO1, iocc=i, vr=coor_no1)
        call nlget(sd_nl, _COOR_NO2, iocc=i, vr=coor_no2)
        axe(1) = (coor_no2(1)-coor_no1(1))
        axe(2) = (coor_no2(2)-coor_no1(2))
        axe(3) = (coor_no2(3)-coor_no1(3))

        r8bid = axe(1)**2+axe(2)**2+axe(3)**2
        if (r8bid .le. r8miem()) then
            call utmess('F', 'DISCRETS_43')
        end if

        call angvx(axe, alpha, beta)
        call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=i, rvect=[sin(alpha), cos(alpha)])
        call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=i, rvect=[sin(beta), cos(beta)])
        call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=i, rvect=[0.d0, 1.d0])

        obst_typ = 'BI_PLANY'
    else
        obst_typ = 'MONO'
    end if
    call nlsav(sd_nl, _OBST_TYP, 1, iocc=i, kscal=obst_typ)
!
!   -- 3.6 - Modal displacements of the node(s)
!            Note : if a single node is used, we fill with zeros the
!                   deformations for node_2, this simplifies the
!                   case treatments for calculating the forces
    call dtmget(sd_dtm, _BASE_VEC, vr=bmodal_v)
    call dtmget(sd_dtm, _NB_PHYEQ, iscal=neq)
    call nlinivec(sd_nl, _MODAL_DEPL_NO1, nbddl*nbmode, iocc=i, vr=defmod1)
    if (lnoeu2) call nlinivec(sd_nl, _MODAL_DEPL_NO2, nbddl*nbmode, iocc=i, vr=defmod2)

    do j = 1, nbmode
        do iddl = 1, nbddl
            defmod1(nbddl*(j-1)+iddl) = bmodal(ddlnodes(iddl), j)
        end do

        if (obst_typ(1:2) .eq. 'BI') then
            do iddl = 1, nbddl
                defmod2(nbddl*(j-1)+iddl) = bmodal(ddlnodes(nbddl+iddl), j)
            end do
        end if
    end do
!
!   --- 3.7 - Multi supported systems, filling up of psixdelta for the
!             concerned nodes
    call dtmget(sd_dtm, _MULTI_AP, kscal=monmot)
    if (monmot(1:3) .eq. 'OUI') then
        call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
        call dtmget(sd_dtm, _NB_EXC_T, iscal=nexcit)
        call jeveuo(nomres//'.IPSD', 'E', vr=ps1del_v)

        call nlinivec(sd_nl, _PSI_DELT_NO1, nbddl*nexcit, iocc=i, vr=ps2del1)
        do j = 1, nexcit
            do iddl = 1, nbddl
                ps2del1(nbddl*(j-1)+iddl) = ps1del(ddlnodes(iddl), j)
            end do
        end do
        if (obst_typ(1:2) .eq. 'BI') then
            call nlinivec(sd_nl, _PSI_DELT_NO2, nbddl*nexcit, iocc=i, vr=ps2del2)
            do j = 1, nexcit
                do iddl = 1, nbddl
                    ps2del2(nbddl*(j-1)+iddl) = ps1del(ddlnodes(nbddl+iddl), j)
                end do
            end do
        end if
    end if
!
!   --- 4 - Updating indices for sd_nl and sd_dtm
    call nlsav(sd_nl, _MAX_LEVEL, 1, iscal=nbnoli)
    call dtmsav(sd_dtm, _NB_NONLI, 1, iscal=nbnoli)
!
    call nlget(sd_nl, _NB_DIS_ECRO_CINE, iscal=nbdeci)
    nbdeci = nbdeci+1
    call nlsav(sd_nl, _NB_DIS_ECRO_CINE, 1, iscal=nbdeci)
!
    AS_DEALLOCATE(vi=ddlnodes)
!
    call jedema()
end subroutine
