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
subroutine te0442(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/cylrep.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxsiro.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/normev.h"
#include "asterfort/tecach.h"
#include "asterfort/tpsivp.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT/DKTG/DST/Q4G/Q4GG
!
! Options: REPE_TENS, REPE_GENE
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, npg, ivf, iret(4)
    integer(kind=8) :: jvGeom, jvFieldIn, jvFieldOut, jvAngrep, np, itab(7), iret1, iret2, nbsp
    integer(kind=8) ::  ncmp, vali(2)
    real(kind=8) :: repAlpha, repBeta
    real(kind=8), dimension(3, 3) :: pig, pgcyl, picyl
    real(kind=8), dimension(2, 2) :: t2iu2
    integer(kind=8), parameter :: nptmax = 9, nspmax = 162
    real(kind=8), dimension(3) :: axe_z, orig, x, xsp, xbary
    integer(kind=8) :: repType
    character(len=8) :: fami
    character(len=8) :: paraInName, paraOutName
    integer(kind=8) :: ipt, ino, joff, type_pt, ipaxe, ipaxe2
    real(kind=8) :: a, b, xnorm, epais, excen, zic, hLayer
    integer(kind=8), parameter :: pt_gauss = 1, pt_noeud = 2
    integer(kind=8) :: nbLayer, iLayer, isp
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'REPE_TENS' .or. option .eq. 'REPE_GENE')

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pig)

! - Compute coordinate system for plate
    call compCoorSystPlate(pig, plateCara, plateOrie)

    if (option .eq. 'REPE_TENS') then
        ncmp = 6
        call tecach('ONO', 'PCOGAIN', 'L', iret(1), nval=7, itab=itab)
        call tecach('ONO', 'PCONOIN', 'L', iret(2), nval=7, itab=itab)
        call tecach('ONO', 'PDEGAIN', 'L', iret(3), nval=7, itab=itab)
        call tecach('ONO', 'PDENOIN', 'L', iret(4), nval=7, itab=itab)
        iret1 = iret(1)+iret(2)+iret(3)+iret(4)
        ASSERT(iret1 .eq. 6)
        if (iret(1) .eq. 0) then
            paraInName = 'PCOGAIN'
            paraOutName = 'PCOGAOUT'
        else if (iret(2) .eq. 0) then
            paraInName = 'PCONOIN'
            paraOutName = 'PCONOOUT'
        else if (iret(3) .eq. 0) then
            paraInName = 'PDEGAIN'
            paraOutName = 'PDEGAOUT'
        else if (iret(4) .eq. 0) then
            paraInName = 'PDENOIN'
            paraOutName = 'PDENOOUT'
        end if
    else if (option .eq. 'REPE_GENE') then
        ncmp = 8
        call tecach('ONO', 'PEFGAIN', 'L', iret(1), nval=7, itab=itab)
        call tecach('ONO', 'PEFNOIN', 'L', iret(2), nval=7, itab=itab)
        call tecach('ONO', 'PDGGAIN', 'L', iret(3), nval=7, itab=itab)
        call tecach('ONO', 'PDGNOIN', 'L', iret(4), nval=7, itab=itab)
        iret1 = iret(1)+iret(2)+iret(3)+iret(4)
        ASSERT(iret1 .eq. 6)
        if (iret(1) .eq. 0) then
            paraInName = 'PEFGAIN'
            paraOutName = 'PEFGAOUT'
        else if (iret(2) .eq. 0) then
            paraInName = 'PEFNOIN'
            paraOutName = 'PEFNOOUT'
        else if (iret(3) .eq. 0) then
            paraInName = 'PDGGAIN'
            paraOutName = 'PDGGAOUT'
        else if (iret(4) .eq. 0) then
            paraInName = 'PDGNOIN'
            paraOutName = 'PDGNOOUT'
        end if
    end if

    if (paraInName(4:5) .eq. 'NO') then
        fami = 'NOEU'
        type_pt = pt_noeud
    else if (paraInName(4:5) .eq. 'GA') then
        fami = 'RIGI'
        type_pt = pt_gauss
    end if

! Infos sur les noeuds et points de Gauss de l'élément
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jvf=ivf)

! Nombre de points en fonction de la localisation des champs
    if (paraInName(4:5) .eq. 'NO') then
        np = nno
    else if (paraInName(4:5) .eq. 'GA') then
        np = npg
    end if
    ASSERT(np .le. nptmax)

! - Get input and output fields
    call jevech(paraInName, 'L', jvFieldIn)
    call tecach('OOO', paraInName, 'L', iret2, nval=7, itab=itab)
    nbsp = itab(7)
    if ((nbsp .ne. 1) .and. (mod(nbsp, 3) .ne. 0)) then
        call utmess('F', 'ELEMENTS5_54', si=nbsp)
    end if
    vali(1) = nspmax
    vali(2) = nbsp
    if (nbsp .gt. nspmax) then
        call utmess('F', 'ELEMENTS5_4', ni=2, vali=vali)
    end if
    call jevech(paraOutName, 'E', jvFieldOut)

! - Paramètres de la coque
    epais = plateCara%thick
    excen = plateCara%offset
    nbLayer = plateCara%nbLayer

! - Definition of new frame for output field
    call jevech('PANGREP', 'L', jvAngrep)
    repType = nint(zr(jvAngrep-1+3))

    if (repType == 0) then
        repAlpha = zr(jvAngrep-1+1)*r8dgrd()
        repBeta = zr(jvAngrep-1+2)*r8dgrd()
        call coqrep(pig, repAlpha, repBeta, t2iu2)
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, plateOrie%t2ui, zr(jvFieldIn), zr(jvFieldOut))
        else if (option .eq. 'REPE_GENE') then
            call dxefro(np, plateOrie%t2ui, zr(jvFieldIn), zr(jvFieldOut))
        end if
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, t2iu2, zr(jvFieldOut), zr(jvFieldOut))
        else if (option .eq. 'REPE_GENE') then
            call dxefro(np, t2iu2, zr(jvFieldOut), zr(jvFieldOut))
        end if

    else if (repType == 1) then
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, plateOrie%t2iu, zr(jvFieldIn), zr(jvFieldOut))
        else if (option .eq. 'REPE_GENE') then
            call dxefro(np, plateOrie%t2iu, zr(jvFieldIn), zr(jvFieldOut))
        end if

    else if (repType == 2) then
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, plateOrie%t2ui, zr(jvFieldIn), zr(jvFieldOut))
        else if (option .eq. 'REPE_GENE') then
            call dxefro(np, plateOrie%t2ui, zr(jvFieldIn), zr(jvFieldOut))
        end if

    else if (repType == 3) then
!   PASSAGE DES CONTRAINTES DU REPERE LOCAL 1
!   A L'ELEMENT AU REPERE CYLINDRIQUE
!   REPERE ='COQUE_UTIL_CYL'
        if (option .eq. 'REPE_TENS') then
            ASSERT(nbsp /= np*3*nbLayer)
            hLayer = epais/nbLayer
! Définition du repère
            axe_z = zr(jvAngrep-1+4:jvAngrep-1+6)
            call normev(axe_z, xnorm)
! Si l'axe z n'est pas initialisé dans la carte, on ne veut pas changer de repère
! => on ne fait rien
            if (xnorm > r8prem()) then
                orig = zr(jvAngrep-1+7:jvAngrep-1+9)
! Passage des contraintes du repère local 1 (repère utilisateur défini par
! cara_elem) au repère intrinsèque pour tous les points et sous-points
! matrice de passage plateOrie%t2ui
                call dxsiro(np*nbsp, plateOrie%t2ui, zr(jvFieldIn), zr(jvFieldOut))
! Passage des contraintes du repère intrinsèque au repère local cylindrique
                joff = 0
                do ipt = 1, np
! -- Calcul des coordonnées du point
                    x = 0.d0
                    if (type_pt == pt_gauss) then
!  Le point est un point de Gauss
                        do ino = 1, nno
                            x(1) = x(1)+zr(jvGeom+3*(ino-1)-1+1)*zr(ivf+(ipt-1)*nno+ino-1)
                            x(2) = x(2)+zr(jvGeom+3*(ino-1)-1+2)*zr(ivf+(ipt-1)*nno+ino-1)
                            x(3) = x(3)+zr(jvGeom+3*(ino-1)-1+3)*zr(ivf+(ipt-1)*nno+ino-1)
                        end do
                    else if (type_pt == pt_noeud) then
!  Le point est un noeud de l'élément
                        x(1) = zr(jvGeom+3*(ipt-1)-1+1)
                        x(2) = zr(jvGeom+3*(ipt-1)-1+2)
                        x(3) = zr(jvGeom+3*(ipt-1)-1+3)
                    end if

                    do iLayer = 1, nbLayer
                        do isp = 1, 3
                            zic = excen-epais/2.d0+(iLayer-1)*hLayer
                            if (isp .eq. 1) then
                                zic = zic
                            else if (isp .eq. 2) then
                                zic = zic+hLayer/2.d0
                            else
                                zic = zic+hLayer
                            end if
                            xsp(:) = x(:)+zic*pig(3, :)

!-- Calcul de la matrice de passage du repère global au repère local cylindrique
                            ipaxe = 0
                            call cylrep(ndim, xsp, axe_z, orig, pgcyl, &
                                        ipaxe)
                            if (ipaxe > 0) then
! le point est sur l'axe du repere cylindrique, on essaie de se placer au
! centre de gravité pour calculer le repère
                                xbary(:) = 0
                                do ino = 1, nno
                                    xbary(1) = xbary(1)+zr(jvGeom+3*(ino-1)-1+1)
                                    xbary(2) = xbary(2)+zr(jvGeom+3*(ino-1)-1+2)
                                    xbary(3) = xbary(3)+zr(jvGeom+3*(ino-1)-1+3)
                                end do
                                xbary(:) = xbary(:)/nno
                                ipaxe2 = 0
                                call cylrep(ndim, xbary, axe_z, orig, pgcyl, &
                                            ipaxe2)
                                if (ipaxe2 > 0) then
                                    call utmess('A', 'ALGORITH2_13')
                                end if
                            end if
! picyl: matrice de passage du repère intrinsèque au repère cylindrique
! picyl = pig*pgcyl
                            a = 1.d0
                            b = 0.d0
                            b_ldc = to_blas_int(3)
                            b_ldb = to_blas_int(3)
                            b_lda = to_blas_int(3)
                            b_m = to_blas_int(3)
                            b_n = to_blas_int(3)
                            b_k = to_blas_int(3)
                            call dgemm('N', 'N', b_m, b_n, b_k, &
                                       a, pig(1, 1), b_lda, pgcyl(1, 1), b_ldb, &
                                       b, picyl(1, 1), b_ldc)
! Appliquer le changement de base
                            call tpsivp(picyl, zr(jvFieldOut+joff-1+1:jvFieldOut+joff-1+ncmp))
! Mettre à jour l'offset (un tenseur 3D symétrique a 6 composantes)
                            joff = joff+ncmp
                        end do
                    end do
                end do
            else
! on ne veut pas changer de repère, on recopie le champ
                zr(jvFieldOut:jvFieldOut+np*nbsp*ncmp) = zr(jvFieldIn:jvFieldIn+np*nbsp*ncmp)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
end subroutine te0442
