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
subroutine te0516(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/bsigma.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jeexin.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/matela.h"
#include "asterfort/matrot.h"
#include "asterfort/mavec.h"
#include "asterfort/moytem.h"
#include "asterfort/pmfdef.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/pmfite.h"
#include "asterfort/pmfitg.h"
#include "asterfort/pmfits.h"
#include "asterfort/pmfmats.h"
#include "asterfort/pmfmcf.h"
#include "asterfort/porea1.h"
#include "asterfort/porea3.h"
#include "asterfort/pouex7.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "blas/dscal.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!           ELEMENTS DE POUTRE MULTI-FIBRES DE TIMOSHENKO AVEC GAUCHISSEMENT.
!
!       OPTION       RAPH_MECA FULL_MECA RIGI_MECA_TANG
!       NOMTE        MECA_POU_D_TGM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: nc, nno, dimklv, npg, iret, codrep
    parameter(nc=7, dimklv=2*nc*(2*nc+1)/2, nno=2, npg=3)
    real(kind=8) :: hoel(nc), fl(2*nc), hota(nc, nc), d1b(nc, 2*nc)
    real(kind=8) :: rg0(2*nc, 2*nc), eps(nc), deps(nc), u(2*nc), du(2*nc)
    real(kind=8) :: klv(dimklv), work(nc, 2*nc), co(npg)
    real(kind=8) :: rigge0(2*nc, 2*nc), ddu(2*nc), effgep(nc), d1bsig(4, 2*nc)
    integer(kind=8) :: ne, cara, idepla, ifgp
    integer(kind=8) :: ii, jcret, npge
    integer(kind=8) :: jvGeom, jvMaterc, icontm, iorien, ivarim, jvInstpr, ipoids
    integer(kind=8) :: jvCarcri, ideplm, ideplp, jvInstmr, ivectu, icontp, ivarip, imat
    integer(kind=8) :: jacf, ivarmp, codret
    integer(kind=8) :: ncomp, nbvalc, isdcom
    integer(kind=8) :: kpg, jj, kk, istrxm, istrxp, istrmp, ncomp2
    real(kind=8) :: ey, ez, gamma, xl, xls2, Nx, My, Mz
    real(kind=8) :: aa, xiy, xiz, alfay, alfaz, xjx, xjg
    real(kind=8) :: e, g, nu, temp, phiy, phiz
    real(kind=8) :: defam(6), defap(6), angp(3)
    real(kind=8) :: pgl(3, 3), ffp(3), matsct(6)
    real(kind=8) :: xiyr2, xizr2, hotage(4, 4), epsm
    real(kind=8) :: ksi1, d1b3(2, 3), sigfib, mflex(4)
    real(kind=8) :: carsec(6)
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=8) :: materSect
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: reactu, isgrot
    character(len=16) :: defoComp, multComp, rigiGeom, typeComp
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10)
    real(kind=8), pointer :: defmfib(:) => null()
    real(kind=8), pointer :: defpfib(:) => null()
    real(kind=8), pointer :: modufib(:) => null()
    real(kind=8), pointer :: vsigfib(:) => null()
    real(kind=8), pointer :: varfib(:) => null()
    integer(kind=8), parameter :: nb_cara = 8
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8), parameter :: noms_cara(nb_cara) = (/'AY1  ', 'AZ1  ', 'EY1  ', 'EZ1  ', &
                                                          'JX1  ', 'JG1  ', 'IYR21', 'IZR21'/)
    blas_int :: b_incx, b_n
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, &
                     npg=npge, jpoids=ipoids)
    ASSERT(npg .ge. npge)
    co(1:npg) = zr(ipoids:ipoids+npg-1)
!
    hoel(:) = 0.0d0
    fl(:) = 0.0d0
    rg0(:, :) = 0.0d0
    rigge0(:, :) = 0.0d0
    mflex(:) = 0.d0
    codret = 0
    codrep = 0

!   Récupération des caractéristiques des fibres
    call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, &
                 jacf=jacf)
!
!   Nombre de composantes du champs PSTRX?? par points de gauss
!   La 15eme composante ne concerne pas les POU_D_TGM
    ncomp = 18
!   Nombre de composantes d'efforts et de déformations généralisés
    ncomp2 = 7

!   Longueur de l'élément et pointeur sur le géométrie
    xl = lonele(igeom=jvGeom)

!  Paramètres en entrée
    call jevech('PCAORIE', 'L', iorien)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PSTRXMR', 'L', istrxm)
!   la presence du champ de deplacement a l instant t+ devrait etre conditionne  par l'option
!   (avec rigi_meca_tang ca n a pas de sens). Ce champ est initialise a 0 par la routine nmmatr.
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PDDEPLA', 'L', idepla)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)

!   Pour matric=(RIGI_MECA_TANG) : valeurs "+" égalent valeurs "-"
    icontp = icontm
    ivarip = ivarim
    istrxp = istrxm
    ivarmp = ivarim

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Get fields for non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Properties of behaviour
    defoComp = compor(DEFO)
    multComp = compor(MULTCOMP)
    rigiGeom = compor(RIGI_GEOM)
    typeComp = compor(INCRELAS)
    read (compor(NVAR), '(I16)') nbvalc
    reactu = defoComp .eq. 'GROT_GDEP' .or. rigiGeom .eq. 'OUI'
    isgrot = defoComp .eq. 'GROT_GDEP'

! - Some checks
    if (typeComp .eq. 'COMP_ELAS') then
        call utmess('F', 'POUTRE0_90')
    else if ((defoComp .ne. 'PETIT') .and. (defoComp .ne. 'GROT_GDEP')) then
        call utmess('F', 'POUTRE0_40', sk=defoComp)
    end if
    if ((defoComp .ne. 'PETIT') .and. (defoComp .ne. 'GROT_GDEP')) then
        call utmess('F', 'POUTRE0_40', sk=defoComp)
    end if
    call jeexin(multComp, iret)
    if (iret .eq. 0) then
        call utmess('F', 'POUTRE0_14', sk=nomte)
    end if

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)
    if (option .eq. 'RIGI_MECA_TANG') then
        lVect = ASTER_FALSE
        lSigm = ASTER_FALSE
        lVari = ASTER_FALSE
    end if

! - Get time
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)

!   Recuperation de la SD_COMPOR ou le comportement des groupes de fibres est stocke
!   pour chaque groupe : (nom, mater, loi, algo1d, deformation nbfig) dans
!   l'ordre croissant des numeros de groupes
    call jeveuo(multComp, 'L', isdcom)

!   deformations anelastiques
    defam = 0.d0
    defap = 0.d0

!  Paramètres en sortie
    if (lMatr) then
        call jevech('PMATUUR', 'E', imat)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIMP', 'L', ivarmp)
        call jevech('PSTRXMP', 'L', istrmp)
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PSTRXPR', 'E', istrxp)
    end if

! - Update geometry
    if (reactu) then
!       recuperation du 3eme angle nautique au temps t-
        gamma = zr(istrxm+18-1)
!       calcul de PGL,XL et ANGP
        if (isgrot) then
!           Rotation sur iteration non convergee
            call porea1(nno, nc, zr(ideplm), zr(ideplp), zr(jvGeom+1), &
                        gamma, lVect, pgl, xl, angp)
        else
!          Rotation sur dernier pas convergee
            call porea3(nno, nc, zr(ideplm), zr(ideplp), zr(jvGeom+1), &
                        gamma, pgl, xl, angp)
        end if
!       sauvegarde des angles nautiques
        if (lVect) then
            zr(istrxp+16-1) = angp(1)
            zr(istrxp+17-1) = angp(2)
            zr(istrxp+18-1) = angp(3)
        end if
    else
        call matrot(zr(iorien), pgl)
    end if
    xls2 = xl/2.d0

!   recuperation des caracteristiques de la section
    call pmfitg(tygrfi, nbfibr, nbcarm, zr(jacf), carsec)
    aa = carsec(1)
    xiy = carsec(5)
    xiz = carsec(4)
    call poutre_modloc('CAGNP2', noms_cara, nb_cara, lvaleur=vale_cara)
    alfay = vale_cara(1)
    alfaz = vale_cara(2)
    ey = vale_cara(3)
    ez = vale_cara(4)
    xjx = vale_cara(5)
    xjg = vale_cara(6)
    xiyr2 = vale_cara(7)
    xizr2 = vale_cara(8)

!   calcul des deplacements et de leurs increments passage dans le repere local
    call utpvgl(nno, nc, pgl, zr(ideplm), u)
    call utpvgl(nno, nc, pgl, zr(ideplp), du)
    call utpvgl(nno, nc, pgl, zr(idepla), ddu)
!   prise en compte de la position du centre de torsion
    do ii = 1, 2
        u(7*(ii-1)+2) = u(7*(ii-1)+2)-ez*u(7*(ii-1)+4)
        u(7*(ii-1)+3) = u(7*(ii-1)+3)+ey*u(7*(ii-1)+4)
        du(7*(ii-1)+2) = du(7*(ii-1)+2)-ez*du(7*(ii-1)+4)
        du(7*(ii-1)+3) = du(7*(ii-1)+3)+ey*du(7*(ii-1)+4)
        ddu(7*(ii-1)+2) = ddu(7*(ii-1)+2)-ez*ddu(7*(ii-1)+4)
        ddu(7*(ii-1)+3) = ddu(7*(ii-1)+3)+ey*ddu(7*(ii-1)+4)
    end do
!   coefficient dependant de la temperature moyenne
    call moytem(fami, npg, 1, '+', temp, iret)

!   caracteristiques elastiques (pas de temperature pour l'instant)
!   on prend le E et NU du materiau torsion (voir op0059)
    call pmfmats(materSect)
    ASSERT(materSect .ne. ' ')
    call matela(zi(jvMaterc), materSect, 1, temp, e, nu)
    g = e/(2.d0*(1.d0+nu))

!   matrice de raideur elastique : materiau integre sur la section
    hoel(1) = e*aa
    hoel(2) = g*aa/alfay
    hoel(3) = g*aa/alfaz
    hoel(4) = g*xjx
    hoel(5) = e*xiy
    hoel(6) = e*xiz
    hoel(7) = e*xjg
    phiy = e*xiz*12.d0*alfay/(xl*xl*g*aa)
    phiz = e*xiy*12.d0*alfaz/(xl*xl*g*aa)
!   deformatiions moins et increment de deformation pour chaque fibre
    AS_ALLOCATE(vr=defmfib, size=nbfibr)
    AS_ALLOCATE(vr=defpfib, size=nbfibr)
!   module et contraintes sur chaque fibre (comportement)
    AS_ALLOCATE(vr=modufib, size=nbfibr)
    AS_ALLOCATE(vr=vsigfib, size=nbfibr)
    AS_ALLOCATE(vr=varfib, size=nbfibr*nbvalc*npg)
!   boucle sur les points de gauss
    do kpg = 1, 3
!       calcul de EPS, DEPS, D1B ( EPSI = D1B * U ) :
        call jsd1ff(kpg, xl, phiy, phiz, d1b)
        eps(1:nc) = 0.0d0
        deps(1:nc) = 0.0d0
!       calcul de l'increment de deformation sur un pas en grands deplacements il est cumulatif.
!           - dans strxmp, on trouve l'increment de deformation jusqu'a l'iteration de newton
!             precedente (0 pour la première itération : strxmp n'est pas présent dans ce cas)
!           - dans ddu, on trouve l'increment de deplacement depuis l'iteration de newton
!             precedente (0 pour la première itération)
!           - apres calcul de deps, on le stocke dans strxpr
!       les deformations sont stockes a partir de la 8eme position
        kk = ncomp*(kpg-1)+ncomp2
        if (.not. reactu) then
!           calcul classique des deformations à partir de DU
            do ii = 1, nc
                do jj = 1, 2*nc
                    eps(ii) = eps(ii)+d1b(ii, jj)*u(jj)
                    deps(ii) = deps(ii)+d1b(ii, jj)*du(jj)
                end do
            end do
        else
!           calcul ameliore tenant compte de la reactualisation
!           on cumule les increments de def de chaque iteration
            if (.not. lVect) then
                do ii = 1, nc
                    do jj = 1, 2*nc
                        eps(ii) = eps(ii)+d1b(ii, jj)*u(jj)
                    end do
                    deps(ii) = 0.d0
                end do
            else
                do ii = 1, nc
                    deps(ii) = zr(istrmp+kk+ii)
                    do jj = 1, 2*nc
                        eps(ii) = eps(ii)+d1b(ii, jj)*u(jj)
                        deps(ii) = deps(ii)+d1b(ii, jj)*ddu(jj)
                    end do
                    zr(istrxp+kk+ii) = deps(ii)
                end do
            end if
        end if
!       calcul des deformations et des increments de def sur les fibres
        call pmfdef(tygrfi, nbfibr, nbcarm, zr(jacf), eps, &
                    defmfib)
        call pmfdef(tygrfi, nbfibr, nbcarm, zr(jacf), deps, &
                    defpfib)
        epsm = (u(8)-u(1))/xl
!       module et contraintes sur chaque fibre (comportement)
        call pmfmcf(materPara, &
                    option, zr(jvCarcri), &
                    kpg, nbgrfi, nbfibr, nug, zk24(isdcom), &
                    zr(jvInstmr), zr(jvInstpr), &
                    nbvalc, defam, defap, zr(ivarim), zr(ivarmp), &
                    zr(icontm), defmfib, defpfib, epsm, modufib, &
                    vsigfib, varfib, codrep)
!
        if (codrep .ne. 0) then
            codret = codrep
!           code 3: on continue et on le renvoie a la fin. autres codes: sortie immediate
            if (codrep .ne. 3) goto 900
        end if
!       calcul de BT*H*B :
        if (lMatr) then
            hota(1:nc, 1:nc) = 0.0d0
!           calcul de la matrice tangente au comportement global
!           seuls 3 efforts sont concernés, les autres ==> élastique
!              effort normal X : composante 1
!              moment autour Y : composante 5
!              moment autour Z : composante 6
!           calcul de la raideur tangente au comportement par fibre
            call pmfite(tygrfi, nbfibr, nbcarm, zr(jacf), modufib, &
                        matsct)
!           MATSCT(1:3) : INT(E.DS)     INT(E.Y.DS)   INT(E.Z.DS)
!           MATSCT(4:6) : INT(E.Y.Y.DS) INT(E.Z.Z.DS) INT(E.Y.Z.DS)
            hota(2, 2) = hoel(2)
            hota(3, 3) = hoel(3)
            hota(4, 4) = hoel(4)
            hota(7, 7) = hoel(7)
!
            hota(1, 1) = matsct(1)
            hota(1, 5) = matsct(3)
            hota(1, 6) = -matsct(2)
            hota(5, 1) = matsct(3)
            hota(5, 5) = matsct(5)
            hota(5, 6) = -matsct(6)
            hota(6, 1) = -matsct(2)
            hota(6, 5) = -matsct(6)
            hota(6, 6) = matsct(4)
            b_n = to_blas_int(nc*nc)
            b_incx = to_blas_int(1)
            call dscal(b_n, xls2, hota, b_incx)
            b_n = to_blas_int(nc*nc)
            b_incx = to_blas_int(1)
            call dscal(b_n, co(kpg), hota, b_incx)
            call utbtab('CUMU', nc, 2*nc, hota, d1b, &
                        work, rg0)
        end if
!       On stocke a "+" : contraintes, fl, vari
        if (lSigm) then
!           Contraintes
            do ii = 1, nbfibr
                zr(icontp-1+nbfibr*(kpg-1)+ii) = vsigfib(ii)
            end do
        end if
        if (lVari) then
!           Variables internes
            do ii = 1, nbfibr*nbvalc*npg
                zr(ivarip-1+ii) = varfib(ii)
            end do
        end if
        if (lVect) then
            ASSERT(lSigm)
!           Efforts généralisés à "+" :
!               ffp : < int(sig.ds) int(sig.y.ds)  int(sig.z.ds) >
!               ffp : <     Nx          -Mz             My       >
            call pmfits(tygrfi, nbfibr, nbcarm, zr(jacf), vsigfib, &
                        ffp)
            Nx = ffp(1)
            My = ffp(3)
            Mz = -ffp(2)
!           Calcul des efforts généralisés
            ifgp = ncomp*(kpg-1)-1
            zr(istrxp+ifgp+1) = Nx
            zr(istrxp+ifgp+2) = zr(istrxm+ifgp+2)+hoel(2)*deps(2)
            zr(istrxp+ifgp+3) = zr(istrxm+ifgp+3)+hoel(3)*deps(3)
!           On rajoute l'effet WAGNER dû au gauchissement
            zr(istrxp+ifgp+4) = zr(istrxm+ifgp+4)+hoel(4)*deps(4)+(Nx*((xiy+xiz)/aa+ey**2+ez**2)&
                                &)*deps(4)+(My*(xizr2/xiy-2.0d0*ez))*deps(4)-(Mz*(xiyr2/xiz-2.0&
                                &d0*ey))*deps(4)
!
            zr(istrxp+ifgp+5) = My
            zr(istrxp+ifgp+6) = Mz
            zr(istrxp+ifgp+7) = zr(istrxm+ifgp+7)+hoel(7)*deps(7)
!
            do kk = 1, 2*nc
                do ii = 1, nc
                    fl(kk) = fl(kk)+xls2*zr(istrxp+ifgp+ii)*d1b(ii, kk)*co(kpg)
                end do
            end do
        end if
!       Calcul de la matrice de rigidite geometrique
        if (lMatr .and. reactu) then
            hotage(:, :) = 0.0d0
            ifgp = ncomp*(kpg-1)-1
            do ii = 1, ncomp2
                effgep(ii) = zr(istrxp+ifgp+ii)
            end do
            hotage(1, 2) = -effgep(3)
            hotage(1, 3) = effgep(2)
            hotage(1, 4) = -( &
                           ey*effgep(2)+ez*effgep(3))+(0.5d0*(xiyr2/xiz)*effgep(2))+(0.5d0*(xizr&
                           &2/xiy)*effgep(3) &
                           )
!           Terme non calcule exactement (on fait l'hypothese d'une torsion de saint-venant)
            hotage(2, 1) = hotage(1, 2)
            hotage(2, 2) = effgep(1)
            hotage(2, 4) = (ez*effgep(1)-effgep(5))
!
            hotage(3, 1) = hotage(1, 3)
            hotage(3, 3) = effgep(1)
            hotage(3, 4) = -(ey*effgep(1)+effgep(6))
!
            hotage(4, 1) = hotage(1, 4)
            hotage(4, 2) = hotage(2, 4)
            hotage(4, 3) = hotage(3, 4)
!           Moment de WAGNER
            hotage(4, 4) = ( &
                           effgep(1)*((xiy+xiz)/aa+ey**2+ez**2))+(effgep(5)*(xizr2/xiy-2.0d0*ez)&
                           &)-(effgep(6)*(xiyr2/xiz-2.0d0*ey) &
                           )
!           Terme non calcule actuellement car xiwr2 n'est pas fourni par l'utilisateur
!               XIWR2 = INT(W*(Y*Y+Z*Z)*DS) + (EFFGEP(7)*(XIWR2/XJG))
            b_n = to_blas_int(4*4)
            b_incx = to_blas_int(1)
            call dscal(b_n, xls2, hotage, b_incx)
            b_n = to_blas_int(4*4)
            b_incx = to_blas_int(1)
            call dscal(b_n, co(kpg), hotage, b_incx)
!           Recuperation de la matrice des fonctions de forme d1bsig
!           Le dernier argument permet de choisir l'interpolation :
!               lineaire (0)
!               cubique flexion-torsion(1)
            call bsigma(kpg, xl, phiy, phiz, d1bsig, &
                        1)
            call utbtab('CUMU', 4, 2*nc, hotage, d1bsig, &
                        work, rigge0)
        end if
    end do
!
    if (lMatr) then
        if (reactu) then
            if (isgrot) then
!           Calcul de la matrice de correction des GR
!           rappel :
!               le calcul de la matrice de correction kc est fait a part, on tient compte à
!               posteriori des rotations modérées entre deux iterations
!
!           Les MFY et MFZ intervenant ici sont ceux aux extremites et on ne les connait
!           qu'aux points de gauss il faut donc utiliser des fonctions de forme pour les
!           transporter aux noeuds (on prend l'interpolation polynomiale d'ordre 2)
!
!           On projette avec des fcts de forme sur les noeuds debut et fin de l'élément
!           pour le point 1
                ksi1 = -sqrt(5.d0/3.d0)
                d1b3(1, 1) = ksi1*(ksi1-1.d0)/2.0d0
                d1b3(1, 2) = 1.d0-ksi1*ksi1
                d1b3(1, 3) = ksi1*(ksi1+1.d0)/2.0d0
!           pour le point 2
                ksi1 = sqrt(5.d0/3.d0)
                d1b3(2, 1) = ksi1*(ksi1-1.d0)/2.0d0
                d1b3(2, 2) = 1.d0-ksi1*ksi1
                d1b3(2, 3) = ksi1*(ksi1+1.d0)/2.0d0
!
!           pour les noeuds 1 et 2 :
!               calcul des contraintes
!               calcul des efforts generalises a partir des contraintes
                do ne = 1, 2
                    do ii = 1, nbfibr
                        sigfib = 0.d0
                        do kpg = 1, 3
                            kk = icontp+nbfibr*(kpg-1)+ii-1
                            sigfib = sigfib+zr(kk)*d1b3(ne, kpg)
                        end do
                        kk = 2*(ne-1)
                        cara = jacf+(ii-1)*nbcarm
                        mflex(1+kk) = mflex(1+kk)+sigfib*zr(cara+2)*zr(cara+1)
                        mflex(2+kk) = mflex(2+kk)-sigfib*zr(cara+2)*zr(cara)
                    end do
                end do
!           on calcule la matrice tangente en sommant les termes de :
!               rigidité matérielle RG0
!               rigidité géométrique RIGGE0
!               matrice de correction pour la prise en compte de rotations modérées
                rigge0(4, 5) = rigge0(4, 5)+mflex(2)*0.5d0
                rigge0(4, 6) = rigge0(4, 6)-mflex(1)*0.5d0
                rigge0(5, 4) = rigge0(5, 4)+mflex(2)*0.5d0
                rigge0(6, 4) = rigge0(6, 4)-mflex(1)*0.5d0
!
                rigge0(11, 12) = rigge0(11, 12)-mflex(4)*0.5d0
                rigge0(11, 13) = rigge0(11, 13)+mflex(3)*0.5d0
                rigge0(12, 11) = rigge0(12, 11)-mflex(4)*0.5d0
                rigge0(13, 11) = rigge0(13, 11)+mflex(3)*0.5d0
            end if
!       On remet tout dans rg0
            rg0 = rg0+rigge0
        end if
        call mavec(rg0, 2*nc, klv, dimklv)
    end if
!
!   On rend le FL dans le repere global
    if (lVect) then
! Prise en compte du centre de torsion
        do ii = 1, 2
            fl(7*(ii-1)+4) = fl(7*(ii-1)+4)-ez*fl(7*(ii-1)+2)+ey*fl(7*(ii-1)+3)
        end do
! passage local -> global
        call utpvlg(nno, nc, pgl, fl, zr(ivectu))
    end if
!
!   On sort la matrice tangente
    if (lMatr) then
!       Prise en compte du centre de torsion
        call pouex7(klv, ey, ez)
!       Passage local -> global
        call utpslg(nno, nc, pgl, klv, zr(imat))
    end if
!
900 continue
    if (lSigm) then
        zi(jcret) = codret
    end if
    AS_DEALLOCATE(vr=defmfib)
    AS_DEALLOCATE(vr=defpfib)
    AS_DEALLOCATE(vr=modufib)
    AS_DEALLOCATE(vr=vsigfib)
    AS_DEALLOCATE(vr=varfib)
end subroutine
