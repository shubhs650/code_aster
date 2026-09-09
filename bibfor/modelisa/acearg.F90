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

subroutine acearg(nbocc, infdonn, infcarte, zjdlm)
!
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!
!     AFFECTATION DES CARACTERISTIQUES POUR LES ELEMENTS DISCRET PAR RAIDEUR REPARTIE
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
    implicit none
    type(cara_elem_info) :: infdonn
    integer(kind=8) :: nbocc, zjdlm(*)
    type(cara_elem_carte) :: infcarte(*)
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/affdis.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/in_liste_entier.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!&<
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbval
    parameter(nbval=5)
    integer(kind=8) :: jdc(3), jdv(3), ifm
    integer(kind=8) :: jdcinf, jdvinf
    integer(kind=8) :: iv, ioc, ll
    integer(kind=8) :: nbno, ncmp
    integer(kind=8) :: ndim, dimcar
! --------------------------------------------------------------------------------------------------

    character(len=1)  :: kma(3)
    data kma/'K', 'M', 'A'/
    character(len=8)  :: nomu, noma
    character(len=19) :: cart(3), cartdi
!
    aster_logical :: is_uniform, is_quadratic, is_line_or_biquad
    integer(kind=8) :: jdme, npg, posi2, posit(9), ndim2, ipg
    real(kind=8) :: xg(2, 9), wg(9), dff(3, 9), ff(9), ksi(2)
    real(kind=8) :: tan_1(3), tan_2(3), jac
    character(len=8) :: typelem, schema
    character(len=16) :: ktyelm
!
    integer(kind=8) :: n_groups, nbparno, ntopo
    integer(kind=8) :: j_typ, j_grma, j_no, j_grno
    integer(kind=8) :: nb_cells, nb_nodes, i_cell, num_cell, nb_no, i_no, i_noe, i_val
    real(kind=8) :: vale(nbval)
    real(kind=8) :: x(9), y(9), z(9), a(3), b(3), c(3), surf, surtot, surf_cell
    real(kind=8) :: i_x, i_y, x_g, y_g, x_p, y_p, rigi(6), mass(4)
    character(len=8) :: typm, method
    character(len=24) :: grma, gr_seg(nbval), gr_centre(nbval)
    character(len=24) :: magrno, magrma, manoma, matyma
    real(kind=8), pointer :: coeno(:) => null()
    integer(kind=8), pointer :: parno(:) => null()
    real(kind=8), pointer :: coord(:) => null()
    integer(kind=8), pointer :: parcell(:) => null()
    blas_int :: b_1, b_3
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    nomu = infdonn%nomu
    noma = infdonn%maillage
    ndim = infdonn%dimmod
    jdme = infdonn%jmodmail
!   Si c'est un maillage partionné ==> PLOUF
    if (infdonn%IsParaMesh) then
        call utmess('F', 'AFFECARAELEM_99', sk='RIGI_GRILLE')
    endif
    ASSERT(ndim .eq. 3)
!
!   Les cartes sont déjà construites : ace_crea_carte
    cartdi = infcarte(ACE_CAR_DINFO)%nom_carte
    jdcinf = infcarte(ACE_CAR_DINFO)%adr_cmp
    jdvinf = infcarte(ACE_CAR_DINFO)%adr_val
    dimcar = infcarte(ACE_CAR_DINFO)%nbr_cmp
!
    cart(1) = infcarte(ACE_CAR_DISCK)%nom_carte
    jdc(1)  = infcarte(ACE_CAR_DISCK)%adr_cmp
    jdv(1)  = infcarte(ACE_CAR_DISCK)%adr_val
!
    cart(2) = infcarte(ACE_CAR_DISCM)%nom_carte
    jdc(2)  = infcarte(ACE_CAR_DISCM)%adr_cmp
    jdv(2)  = infcarte(ACE_CAR_DISCM)%adr_val
!
    cart(3) = infcarte(ACE_CAR_DISCA)%nom_carte
    jdc(3)  = infcarte(ACE_CAR_DISCA)%adr_cmp
    jdv(3)  = infcarte(ACE_CAR_DISCA)%adr_val
!
    ifm = infdonn%ivr(4)

    magrno = noma//'.GROUPENO'
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    matyma = noma//'.TYPMAIL'
    call jeveuo(matyma, 'L', j_typ)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=coord)

!   Boucle sur les occurrences de rigi_grille
    do ioc = 1, nbocc
!
!       Methode UNIFORME : répartition surfacique uniforme
!               DIAG     : répartition par intégration via les fonctions de forme
        call getvtx('RIGI_GRILLE', 'METHODE', iocc=ioc, nbval=1, vect=method)
        if (method(1:8) .eq. 'UNIFORME') then
            is_uniform = ASTER_TRUE
        else if (method(1:4) .eq. 'DIAG') then
            is_uniform = ASTER_FALSE
        else
            ASSERT(ASTER_FALSE)
        end if
        is_quadratic = ASTER_FALSE
        is_line_or_biquad = ASTER_FALSE

        call getvem(noma, 'GROUP_MA', 'RIGI_GRILLE', 'GROUP_MA', ioc, 1, grma, n_groups)
        call getvem(noma, 'GROUP_MA', 'RIGI_GRILLE', 'GROUP_MA_SEG2', ioc, nbval, gr_seg, n_groups)
        call getvem(noma, 'GROUP_NO', 'RIGI_GRILLE', 'GROUP_NO_CENTRE', ioc, nbval, gr_centre, &
                    n_groups)
        call getvr8('RIGI_GRILLE', 'VALE', iocc=ioc, nbval=nbval, vect=vale)

        call jelira(jexnom(magrma, grma), 'LONUTI', nb_cells)
        ASSERT(nb_cells .ne. 0)
        call jeveuo(jexnom(magrma, grma), 'L', j_grma)

!       calcul des surfaces au support des noeuds

        nb_nodes = 0
        do i_cell = 1, nb_cells
            num_cell = zi(j_grma-1+i_cell)
            call jelira(jexnum(manoma, num_cell), 'LONMAX', nb_no)
            nb_nodes = nb_nodes+nb_no
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(j_typ-1+num_cell)), typm)
            call dismoi('DIM_TOPO', typm, 'TYPE_MAILLE', repi=ntopo)
            if (ntopo .ne. 2) then
                call utmess('F', 'AFFECARAELEM_26', sk=grma, ni=2, vali=[ioc, num_cell])
            end if
            call jenuno(jexnum('&CATA.TE.NOMTE', zi(jdme-1+num_cell)), ktyelm)
            if ((ktyelm .eq. "MEC3TR7H") .or. (ktyelm .eq. "MEC3QU9H")) then
                call utmess('F', 'MODELISA6_39', sk='RIGI_GRILLE')
            end if
            if (nb_no .eq. 6 .or. nb_no .eq. 8) then
                is_quadratic = .true.
            else
                is_line_or_biquad = .true.
            end if
        end do
        if (is_quadratic .and. is_line_or_biquad) then
            call utmess('F', 'MODELISA6_41', sk='RIGI_GRILLE')
        end if
!
        b_1 = to_blas_int(1)
        b_3 = to_blas_int(3)

!       Coefficients des noeuds de l interface
        AS_ALLOCATE(vr=coeno, size=nb_nodes)
!       Participation des noeuds de l interface
        AS_ALLOCATE(vi=parno, size=nb_nodes)

        nbparno = 0
        surtot = 0.d0
        do i_cell = 1, nb_cells
            num_cell = zi(j_grma-1+i_cell)
            call jelira(jexnum(manoma, num_cell), 'LONMAX', nb_no)
            call jeveuo(jexnum(manoma, num_cell), 'L', j_no)
            do i_no = 1, nb_no
                i_noe = zi(j_no-1+i_no)
!               On enregistre le numéro du noeud dans parno, s'il n'y est pas déjà
                if (in_liste_entier(i_noe, parno(1:nbparno), posi2)) then
                    posit(i_no) = posi2
                else
                    nbparno = nbparno+1
                    ASSERT(nbparno .le. nb_nodes)
                    parno(nbparno) = i_noe
                    posit(i_no) = nbparno
                end if
                x(i_no) = coord(3*(i_noe-1)+1)
                y(i_no) = coord(3*(i_noe-1)+2)
                z(i_no) = coord(3*(i_noe-1)+3)
            end do
!
            if (nb_no .eq. 3) then
                typelem = 'TR3'
            else if (nb_no .eq. 4) then
                typelem = 'QU4'
            else if (nb_no .eq. 6) then
                typelem = 'TR6'
            else if (nb_no .eq. 8) then
                typelem = 'QU8'
            else if (nb_no .eq. 9) then
                call jenuno(jexnum('&CATA.TE.NOMTE', zi(jdme-1+num_cell)), ktyelm)
                if (ktyelm .ne. "MECA_FACE9") then
                    call utmess('F', 'MODELISA6_40', sk='RIGI_GRILLE')
                end if
                typelem = 'QU9'
            else
                ASSERT(.false.)
            end if
            if (is_uniform) then
                a(1) = x(3)-x(1)
                a(2) = y(3)-y(1)
                a(3) = z(3)-z(1)
                if (nb_no .eq. 3 .or. nb_no .eq. 6) then
                    b(1) = x(2)-x(1)
                    b(2) = y(2)-y(1)
                    b(3) = z(2)-z(1)
                else if (nb_no .eq. 4 .or. nb_no .eq. 8 .or. nb_no .eq. 9) then
                    b(1) = x(4)-x(2)
                    b(2) = y(4)-y(2)
                    b(3) = z(4)-z(2)
                else
                    ASSERT(.false.)
                end if
                call provec(a, b, c)
                surf_cell = ddot(b_3, c, b_1, c, b_1)
                schema = 'NOEU'
            else
                if (nb_no .eq. 3) then
                    schema = 'FPG3'
                else if (nb_no .eq. 4) then
                    schema = 'FPG4'
                else if (nb_no .eq. 6) then
                    schema = 'FPG6'
                else if (nb_no .eq. 8 .or. nb_no .eq. 9) then
                    schema = 'FPG9'
                else
                    ASSERT(.false.)
                end if
            end if

            call elraga(typelem, schema, ndim2, npg, xg, wg)

            do ipg = 1, npg

                ksi = xg(:, ipg)
                if (is_uniform) then
                    surf = surf_cell / nb_no
                else
                    call elrfdf(typelem, ksi, dff)
                    ! calcul des tangentes au point de Gauss
                    tan_1 = 0.d0
                    tan_2 = 0.d0
                    do i_no = 1, nb_no
                        tan_1(1) = tan_1(1)+dff(1, i_no)*x(i_no)
                        tan_1(2) = tan_1(2)+dff(1, i_no)*y(i_no)
                        tan_1(3) = tan_1(3)+dff(1, i_no)*z(i_no)
                        tan_2(1) = tan_2(1)+dff(2, i_no)*x(i_no)
                        tan_2(2) = tan_2(2)+dff(2, i_no)*y(i_no)
                        tan_2(3) = tan_2(3)+dff(2, i_no)*z(i_no)
                    end do

                    ! calcul du jacobien
                    jac = sqrt(abs(ddot(b_3, tan_1, b_1, tan_1, b_1)* &
                                ddot(b_3, tan_2, b_1, tan_2, b_1) &
                                -ddot(b_3, tan_1, b_1, tan_2, b_1)**2))

                    surf = wg(ipg)*jac
                end if

                call elrfvf(typelem, ksi, ff)
                ! ajout des contributions aux noeuds
                do i_no = 1, nb_no
                    if ((nb_no .eq. 3) .or. (nb_no .eq. 4)  .or. (nb_no .eq. 9)) then
                        coeno(posit(i_no)) = coeno(posit(i_no))+ff(i_no)*surf
                        surtot = surtot+ff(i_no)*surf
                    else
                        ! pour TR6 et QU8
                        ! on prend la fonction de forme au carré pour éviter des contributions
                        ! négative ou nuls
                        coeno(posit(i_no)) = coeno(posit(i_no))+ff(i_no)*ff(i_no)*surf
                        surtot = surtot+ff(i_no)*ff(i_no)*surf
                    end if
                end do
            end do
        end do
        nbno = nbparno
        coeno = coeno/surtot
!
!       Numeros de mailles seg2 associé à parno
        AS_ALLOCATE(vi=parcell, size=nbno)

        do i_val = 1, nbval
            parcell = 0
!           récupération des seg2 associés à parno
            call jelira(jexnom(magrma, gr_seg(i_val)), 'LONUTI', nb_cells)
            if (nb_cells .ne. nbno) then
                call utmess('F', 'AFFECARAELEM_27', sk=gr_seg(i_val), &
                            ni=3, vali=[ioc, nb_cells, nbno])
            end if
            call jeveuo(jexnom(magrma, gr_seg(i_val)), 'L', j_grma)
            do i_cell = 1, nb_cells
                num_cell = zi(j_grma-1+i_cell)
                call jelira(jexnum(manoma, num_cell), 'LONMAX', nb_no)
                if (nb_no .ne. 2) then
                    call utmess('F', 'AFFECARAELEM_28', sk=gr_seg(i_val), &
                                ni=2, vali=[ioc, num_cell])
                end if
                call jeveuo(jexnum(manoma, num_cell), 'L', j_no)
!               Vérification que un des noeuds fait partie de la surface
                do i_no = 1, nbno
                    if (zi(j_no) .eq. parno(i_no) .or. zi(j_no+1) .eq. parno(i_no)) then
                        if (parcell(i_no) .ne. 0) then
                            call utmess('F', 'AFFECARAELEM_29', sk=gr_seg(i_val), &
                                        ni=4, vali=[ioc, parcell(i_no), num_cell, parno(i_no)])
                        end if
                        parcell(i_no) = num_cell
                        goto 22
                    end if
                end do
                call utmess('F', 'AFFECARAELEM_30', sk=gr_seg(i_val), ni=2, vali=[ioc, num_cell])
22              continue
            end do

            if (i_val .gt. 3) then
                ! recup des coordonnées du centre
                call jeveuo(jexnom(magrno, gr_centre(i_val)), 'L', j_grno)
                i_noe = zi(j_grno)
                if (i_val .eq. 4) then
                    y_g = coord(3*(i_noe-1)+2)
                    i_x = 0.d0
                else
                    x_g = coord(3*(i_noe-1)+1)
                    i_y = 0.d0
                end if
                ! calcul du moment
                do i_no = 1, nbno
                    i_noe = parno(i_no)
                    if (i_val .eq. 4) then
                        y_p = coord(3*(i_noe-1)+2)
                        i_x = i_x + coeno(i_no)*(y_p - y_g)**2
                    else
                        x_p = coord(3*(i_noe-1)+1)
                        i_y = i_y + coeno(i_no)*(x_p - x_g)**2
                    end if
                end do
            end if

            ! calcul de la rigidité
            do i_no = 1, nbno
                rigi = 0.d0
                if (i_val .eq. 4) then
                    rigi(3) = vale(i_val) * coeno(i_no) / i_x
                else if (i_val .eq. 5) then
                    rigi(3) = vale(i_val) * coeno(i_no) / i_y
                else
                    rigi(i_val) = vale(i_val) * coeno(i_no)
                end if
                ! on applique la rigidité sur le groupe de SEG_2
                iv = 1
                call affdis(ndim, 1, 0.d0, 'K_TR_D_L', rigi, &
                            jdc, jdv, infdonn%ivr, iv, kma, &
                            ncmp, ll, jdcinf, jdvinf, 1)
                call nocart(cartdi, 3, dimcar, mode='NUM', nma=1, limanu=[parcell(i_no)])
                call nocart(cart(ll), 3, ncmp, mode='NUM', nma=1, limanu=[parcell(i_no)])

!               affectation de matrice masse nulle
                iv = 1
                mass = 0.d0
                call affdis(ndim, 1, 0.d0, 'M_TR_D_L', mass, &
                            jdc, jdv, infdonn%ivr, iv, kma, &
                            ncmp, ll, jdcinf, jdvinf, 1)
                call nocart(cartdi, 3, dimcar, mode='NUM', nma=1, limanu=[parcell(i_no)])
                call nocart(cart(ll), 3, ncmp, mode='NUM', nma=1, limanu=[parcell(i_no)])
            end do

        end do

        AS_DEALLOCATE(vr=coeno)
        AS_DEALLOCATE(vi=parno)
        AS_DEALLOCATE(vi=parcell)

    end do

    call jedema()

end subroutine
