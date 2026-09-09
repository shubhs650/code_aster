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
subroutine tgverm(option, carcri, compor, nno1, nno2, &
                  nno3, geom, ndim, nddl, deplp, &
                  sdepl, vu, vg, vp, vectu, &
                  svect, ncont, contp, scont, nvari, &
                  varip, svari, matuu, smatr, matsym, &
                  epsilo, epsilp, epsilg, varia, iret)
!
    implicit none
!
#include "asterc/r8miem.h"
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mavec.h"
#include "asterfort/r8inir.h"
#include "asterfort/tgveri_use.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    aster_logical :: matsym
    character(len=16) :: option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8) :: iret, nno1, nno2, nno3, ndim
    integer(kind=8) :: vu(3, MT_NNOMAX), vg(MT_NNOMAX), vp(MT_NNOMAX)
    real(kind=8) :: sdepl(*), scont(*), svect(*)
    real(kind=8) :: geom(*), deplp(*), vectu(*), contp(*), matuu(*)
    real(kind=8) :: varip(*), svari(*), smatr(*), varia(*)
!
! ----------------------------------------------------------------------
! VAR OPTION NOM DE L'OPTION DE CALCUL
!             IN  : CELLE UTILISEE PAR LE TE
!             OUT : 'RAPH_MECA' SI BOUCLE, 'FULL_MECA' SI FIN DE BOUCLE
! IN  CARCRI  : CARCRI(1) = type de matrice tangente
!               0 : ANALYTIQUE, on ne passe pas ici
!               1 : PERTURBATION, on calcule Ktgte (FULL_MECA)
!               2 : VERIFICATION, on calcule Ktgte (FULL_MECA) + Kpertu
!               CARCRI(7) = valeur de la perturbation
! OUT IRET   SI IRET = 0 -> FIN, SINON -> BOUCLE
! ----------------------------------------------------------------------
!
! EXEMPLE D'INSERTION DANS UN TE DE L'OPTION FULL_MECA
!  1000 CONTINUE
!       CALL NIFINT(OPTION,...)
!       CALL TGVERM(OPTION,....., IRET)
!       IF (IRET.NE.0) GOTO 1000
!
! ----------------------------------------------------------------------
    character(len=24) :: matra, matrc
    integer(kind=8) :: ematra, ematrc, exi
    integer(kind=8) :: nddl
    integer(kind=8) :: nvari, ncont, iuse
    integer(kind=8) :: i, j, k, l, indi, nvar, init, pos
    real(kind=8) :: v, epsilo, fp, fm, pertu, maxdep, maxgeo, maxpre, maxgon
    real(kind=8) :: matper(nddl*nddl), epsilp, epsilg
    blas_int :: b_incx, b_incy, b_n
    save init, pos
    data matra/'PYTHON.TANGENT.MATA'/
    data matrc/'PYTHON.TANGENT.MATC'/
    data init, pos/1, 0/
! ----------------------------------------------------------------------
!
    call jemarq()
!
!     Calcul de la matrice TGTE par PERTURBATION
!
    iret = 0
    call tgveri_use(option, carcri, compor, iuse)
    if (iuse == 0) then
        goto 999
    end if
!
! --  INITIALISATION (PREMIER APPEL)
!
    if (init .eq. 1) then
!       PERTURBATION OU VERIFICATION => FULL_MECA
        if (option .ne. 'FULL_MECA') then
            goto 999
        end if
!
!       CALCUL de la valeur de la perturbation
!       Ici on est en mecanique seule, les DDL sont
!       seulement des deplacements
!
        maxdep = 0.d0
        maxpre = 0.d0
        maxgon = 0.d0
        maxgeo = 0.d0
        do i = 1, nno1
            do j = 1, ndim
                maxdep = max(maxdep, abs(deplp(vu(j, i))))
            end do
        end do
        do i = 1, nno2
            maxgon = max(maxgon, abs(deplp(vg(i))))
        end do
        do i = 1, nno3
            maxpre = max(maxpre, abs(deplp(vp(i))))
        end do
        do i = 1, nno1*ndim
            maxgeo = max(maxgeo, abs(geom(i)))
        end do
        pertu = carcri(7)
        if (maxdep .gt. pertu*maxgeo) then
            epsilo = pertu*maxdep
        else
            epsilo = pertu*maxgeo
        end if
        if (epsilo .lt. r8miem()) then
            call utmess('F', 'ALGORITH11_86')
        end if
        epsilp = pertu*maxpre
        if (epsilp .lt. r8miem()) then
            call utmess('F', 'ALGORITH11_86')
        end if
        epsilg = pertu*maxgon
        if (epsilg .lt. r8miem()) then
            call utmess('F', 'ALGORITH11_86')
        end if
!      ARCHIVAGE DES VALEURS DE REFERENCE
!
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deplp, b_incx, sdepl, b_incy)
        b_n = to_blas_int(ncont)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, contp, b_incx, scont, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vectu, b_incx, svect, b_incy)
        b_n = to_blas_int(nvari)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, varip, b_incx, svari, b_incy)
!
!       ARCHIVAGE DE LA MATRICE TANGENTE COHERENTE
        if (matsym) then
            k = 0
            do i = 1, nddl
                do j = 1, i
                    v = matuu(k+1)
                    k = k+1
                    smatr((i-1)*nddl+j) = v
                    smatr((j-1)*nddl+i) = v
                end do
            end do
        else
            b_n = to_blas_int(nddl*nddl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, matuu, b_incx, smatr, b_incy)
        end if
!
!      PREPARATION DES ITERATIONS
!
        option = 'RAPH_MECA'
        iret = 1
        init = 0
        pos = 0
    end if
!
! -- TRAITEMENT DES VARIATIONS
!
!
!    SAUVEGARDE DE LA FORCE INTERIEURE PERTURBEE
!
    nvar = int((pos+1)/2)
!
    if (nvar .gt. 0) then
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vectu, b_incx, varia(1+(pos-1)*nddl), b_incy)
    end if
!
    pos = pos+1
    nvar = int((pos+1)/2)
    indi = 1-2*mod(pos, 2)
!
    if (nvar .le. nddl) then
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sdepl, b_incx, deplp, b_incy)
        do i = 1, nno1
            do j = 1, ndim
                if (nvar .eq. vu(j, i)) then
                    deplp(nvar) = sdepl(nvar)+indi*epsilo
                    goto 800
                end if
            end do
        end do
!
        do i = 1, nno2
            if (nvar .eq. vg(i)) then
                deplp(nvar) = sdepl(nvar)+indi*epsilg
                goto 800
            end if
        end do
!
        do i = 1, nno3
            if (nvar .eq. vp(i)) then
                deplp(nvar) = sdepl(nvar)+indi*epsilp
                goto 800
            end if
        end do
!
800     continue
!      INITIALISATION DES CHAMPS 'E'
        call r8inir(ncont, 0.d0, contp, 1)
        call r8inir(nddl, 0.d0, vectu, 1)
        iret = 1
        goto 999
    end if
!
!    CALCUL DE LA MATRICE TANGENTE
!
    do i = 1, nddl
        do j = 1, nddl
            fm = varia((2*j-2)*nddl+i)
            fp = varia((2*j-1)*nddl+i)
            do k = 1, nno1
                do l = 1, ndim
                    if (j .eq. vu(l, k)) then
                        v = (fp-fm)/(2*epsilo)
                        goto 900
                    end if
                end do
            end do
            do k = 1, nno2
                if (j .eq. vg(k)) then
                    v = (fp-fm)/(2*epsilg)
                    goto 900
                end if
            end do
            do k = 1, nno3
                if (j .eq. vp(k)) then
                    v = (fp-fm)/(2*epsilp)
                    goto 900
                end if
            end do
900         continue
            matper((i-1)*nddl+j) = v
        end do
    end do
!
!    MENAGE POUR ARRET DE LA ROUTINE
!
    iret = 0
    init = 1
    option = 'FULL_MECA'
!
!    RETABLISSEMENT DE LA SOLUTION
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sdepl, b_incx, deplp, b_incy)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, svect, b_incx, vectu, b_incy)
    b_n = to_blas_int(ncont)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, scont, b_incx, contp, b_incy)
    b_n = to_blas_int(nvari)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, svari, b_incx, varip, b_incy)
!
!     PERTURBATION => SAUVEGARDE DE LA MATRICE CALCULEE PAR
!     DIFFERENCES FINIES COMME MATRICE TANGENTE
!
    if (abs(carcri(2)-1.d0) .lt. 0.1d0) then
        if (matsym) then
            call mavec(matper, nddl, matuu, nddl*(nddl+1)/2)
        else
            b_n = to_blas_int(nddl*nddl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, matper, b_incx, matuu, b_incy)
        end if
!
!     VERIFICATION
!
    else if (abs(carcri(2)-2.d0) .lt. 0.1d0) then
        if (matsym) then
            call mavec(smatr, nddl, matuu, nddl*(nddl+1)/2)
        else
            b_n = to_blas_int(nddl*nddl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, smatr, b_incx, matuu, b_incy)
        end if
!
!      CREATION DES OBJETS
!      CE N'EST PAS LA PREMIERE FOIS QU'ON CALCULE LA MATRICE TANGENTE
!      -> ON NE CONSERVE QUE LE DERNIER CALCUL (EN COURS)
        call jeexin(matra, exi)
        if (exi .ne. 0) then
            call jedetr(matra)
            call jedetr(matrc)
        end if
!        ON CONSERVE L'ALLOCATION DYNAMIQUE AU DETRIMENT DE L'ALLOCATION
!        STATIQUE, CAR MATRA ET MATRB SONT UTILIES A L'EXTERIEUR DES
!        ROUTINES ELEMENTAIRES
        call wkvect(matra, 'G V R', nddl*nddl, ematra)
        call wkvect(matrc, 'G V R', nddl*nddl, ematrc)
        b_n = to_blas_int(nddl*nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, smatr, b_incx, zr(ematra), b_incy)
        b_n = to_blas_int(nddl*nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, matper, b_incx, zr(ematrc), b_incy)
    end if
!
999 continue
!
    call jedema()
end subroutine
