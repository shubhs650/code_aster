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
subroutine dktnli(plateCara, plateOrie, &
                  BEHInteg, option, typmod, &
                  instm, instp, &
                  xyzl, uml, dul, &
                  btsig, ktan, codret)
!
    use Behaviour_type
    use Behaviour_module
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dkqbf.h"
#include "asterfort/dkqtxy.h"
#include "asterfort/dktbf.h"
#include "asterfort/dkttxy.h"
#include "asterfort/dsxhft.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/dxqloc.h"
#include "asterfort/dxqloc2.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/dxtloc2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/nmcomp.h"
#include "asterfort/pmrvec.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: instm, instp
    real(kind=8), intent(in) :: xyzl(3, 4), uml(6, 4), dul(6, 4)
    real(kind=8), intent(out) :: ktan(576), btsig(6, 4)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Plate elements (DKT)
!
! Compute non-linear options
! Input/output fields are in local frame
!
! --------------------------------------------------------------------------------------------------
!
!     IN  OPT  : OPTION NON LINEAIRE A CALCULER
!                  'RAPH_MECA'
!                  'FULL_MECA'      OU 'FULL_MECA_ELAS'
!                  'RIGI_MECA_TANG' OU 'RIGI_MECA_ELAS'
!     IN  XYZL : COORDONNEES DES NOEUDS
!     IN  UL   : DEPLACEMENT A L'INSTANT T "-"
!     IN  DUL  : INCREMENT DE DEPLACEMENT
!     IN  PGL  : MATRICE DE PASSAGE GLOBAL - LOCAL ELEMENT
!     OUT KTAN : MATRICE DE RIGIDITE TANGENTE
!                    SI 'FULL_MECA' OU 'RIGI_MECA_TANG'
!     OUT BTSIG: DIV (SIGMA)
!                    SI 'FULL_MECA' OU 'RAPH_MECA'
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    character(len=16) :: defoComp, typeComp
    integer(kind=8), parameter :: nbNodeMaxi = 4
    real(kind=8) :: distn
!  CMPS DE DEPLACEMENT :
!   - MEMBRANE : DX(N1), DY(N1), DX(N2), ..., DY(NNO)
!   - FLEXION  : DZ(N1), BETAX(N1), BETAY(N1), DZ(N2), ..., BETAY(NNO)
!  CMPS DE DEFORMATION ET CONTRAINTE PLANE (DANS UNE COUCHE) :
!   -            EPSIXX,EPSIYY,2*EPSIXY
!   -            SIGMXX,SIGMYY,SIGMXY
!  CMPS DE DEFORMATION ET CONTRAINTE PLANE POUR APPEL NMCOMP :
!   -            EPSIXX,EPSIYY,EPSIZZ,SQRT(2)*EPSIXY
!   -            SIGMXX,SIGMYY,SIGMZZ,SQRT(2)*SIGMXY
!  CMPS DE DEFORMATION COQUE :
!   - MEMBRANE : EPSIXX,EPSIYY,2*EPSIXY
!   - FLEXION  : KHIXX,KHIYY,2*KHIXY
!  CMPS D' EFFORTS COQUE :
!   - MEMBRANE : NXX,NYY,NXY
!   - FLEXION  : MXX,MYY,MXY
! --------------------------------------------------------------------
    integer(kind=8) :: nbcou, jnbspi
!            NBCOU:  NOMBRE DE COUCHES (INTEGRATION DE LA PLASTICITE)
!            NPG:    NOMBRE DE POINTS DE GAUSS PAR ELEMENT
!            NC :    NOMBRE DE COTES DE L'ELEMENT
    integer(kind=8), parameter :: npgh = 3
!            NPGH:   NOMBRE DE POINT D'INTEGRATION PAR COUCHE
    integer(kind=8), parameter :: nbcon = 6
!            NBCON:  number of compoennts in stress tensor
    real(kind=8) :: poids, hic, h, zic, zmin, coef
!            POIDS:  POIDS DE GAUSS (Y COMPRIS LE JACOBIEN)
!            AIRE:   SURFACE DE L'ELEMENT
!            HIC:    EPAISSEUR D'UNE COUCHE
!            H :     EPAISSEUR TOTALE DE LA PLAQUE
!            ZIC:    COTE DU POINT D'INTEGRATION DANS L'EPAISSEUR
!            ZMIN:   COTE DU POINT D'INTEGRATION LE PLUS "BAS"
!            INSTM:  INSTANT "-"
!            INSTP:  INSTANT "+"
!            COEF:   POIDS D'INTEGRATION PAR COUCHE
    real(kind=8) :: um(2, nbNodeMaxi), uf(3, nbNodeMaxi), dum(2, nbNodeMaxi), duf(3, nbNodeMaxi)
!            UM:     DEPLACEMENT (MEMBRANE) "-"
!            UF:     DEPLACEMENT (FLEXION)  "-"
!           DUM:     INCREMENT DEPLACEMENT (MEMBRANE)
!           DUF:     INCREMENT DEPLACEMENT (FLEXION)
    real(kind=8) :: eps2d(6), deps2d(6), dsidep(6, 6)
!            EPS2D:  DEFORMATION DANS UNE COUCHE (2D C_PLAN)
!           DEPS2D:  INCREMENT DEFORMATION DANS UNE COUCHE (2D C_PLAN)
!            SIG2D:  CONTRAINTE DANS UNE COUCHE (2D C_PLAN)
!           DSIDEP:  MATRICE D(SIG2D)/D(EPS2D)
    real(kind=8) :: eps(3), khi(3), deps(3), dkhi(3), n(3), m(3), sigmPrep(4)
!    real(kind=8) :: q(2)
!            EPS:    DEFORMATION DE MEMBRANE "-"
!            DEPS:   INCREMENT DE DEFORMATION DE MEMBRANE
!            KHI:    DEFORMATION DE FLEXION  "-"
!            DKHI:   INCREMENT DE DEFORMATION DE FLEXION
!            N  :    EFFORT NORMAL "+"
!            M  :    MOMENT FLECHISSANT "+"
!            Q  :    EFFORT TRANCHANT
!            SIGM : CONTRAINTE "-"
    real(kind=8) :: df(9), dm(9), dmf(9), d2d(9)
!            D2D:    MATRICE DE RIGIDITE TANGENTE MATERIELLE (2D)
!            DF :    MATRICE DE RIGIDITE TANGENTE MATERIELLE (FLEXION)
!            DM :    MATRICE DE RIGIDITE TANGENTE MATERIELLE (MEMBRANE)
!            DMF:    MATRICE DE RIGIDITE TANGENTE MATERIELLE (COUPLAGE)
    real(kind=8) :: bf(3, 3*nbNodeMaxi), bm(3, 2*nbNodeMaxi), bmq(2, 3)
!            BF :    MATRICE "B" (FLEXION)
!            BM :    MATRICE "B" (MEMBRANE)
    real(kind=8) :: flex(3*nbNodeMaxi*3*nbNodeMaxi), memb(2*nbNodeMaxi*2*nbNodeMaxi)
    real(kind=8) :: mefl(2*nbNodeMaxi*3*nbNodeMaxi), work(3*nbNodeMaxi*3*nbNodeMaxi)
!           MEMB:    MATRICE DE RIGIDITE DE MEMBRANE
!           FLEX:    MATRICE DE RIGIDITE DE FLEXION
!           WORK:    TABLEAU DE TRAVAIL
!           MEFL:    MATRICE DE COUPLAGE MEMBRANE-FLEXION
!             LE MATERIAU EST SUPPOSE HOMOGENE
!             IL PEUT NEANMOINS Y AVOIR COUPLAGE PAR LA PLAinstm, instp, STICITE
!     ------------------ PARAMETRAGE ELEMENT ---------------------------
    integer(kind=8), parameter :: ndimLdc = 2
    real(kind=8), parameter :: deux = 2.d0, rac2 = sqrt(2.d0)
    integer(kind=8) :: ndim, nbNode, npg, ipoids, icoopg
    integer(kind=8) :: jtab(7), codkpg, i, ksp
    integer(kind=8) :: jvCarcri, icontm, icontp, icou, icpg, igauh
    integer(kind=8) :: ino, kpg, iret, isp, ivarim, ivarip, ivarix, ivpg
    integer(kind=8) :: j, k, nbsp, nbVari
    real(kind=8) :: qsi, eta, cara(25), jacob(5)
    real(kind=8) :: ctor, coehsd, zmax, quotient, a, b, c
    aster_logical :: lTria, lQuad, leul
    real(kind=8) :: dvt(2), vt(2)
    real(kind=8) :: dfel(3, 3), dmel(3, 3), dmfel(3, 3), dcel(2, 2), dciel(2, 2)
    real(kind=8) :: dmcel(3, 2), dfcel(3, 2)
    real(kind=8) :: d1iel(2, 2)
    real(kind=8) :: depfel(3*nbNodeMaxi)
    real(kind=8) :: hft2el(2, 6)
    aster_logical :: coupmfel, l_matr_symm
    integer(kind=8) :: multicel
    integer(kind=8) :: lg_varip
    real(kind=8), allocatable :: varip(:)
    character(len=8), parameter :: fami = 'RIGI'
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=16), pointer :: compor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    lSigm = L_SIGM(option)
    lVect = L_VECT(option)
    lMatr = L_MATR(option)
    lVari = L_VARI(option)
    flex = 0.d0
    memb = 0.d0
    mefl = 0.d0
    ktan = 0.d0
    btsig = 0.d0
    codret = 0

! - Get finite element parameters
    call elrefe_info(fami=fami, &
                     ndim=ndim, nno=nbNode, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg)
    lTria = ASTER_FALSE
    lQuad = ASTER_FALSE
    if (nbNode .eq. 3) then
        lTria = ASTER_TRUE
    else if (nbNode .eq. 4) then
        lQuad = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get plate parameters
    h = plateCara%thick
    distn = plateCara%offset
    ctor = plateCara%coefRigiDRZ
    if (lTria) then
        call gtria3(xyzl, cara)
    else if (lQuad) then
        call gquad4(xyzl, cara)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get input fields
    call tecach('OOO', 'PCONTMR', 'L', iret, nval=7, itab=jtab)
    nbsp = jtab(7)
    icontm = jtab(1)
    ASSERT(npg .eq. jtab(3))
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCARCRI', 'L', jvCarcri)
    l_matr_symm = .true.
    if (nint(zr(jvCarcri-1+CARCRI_MATRSYME)) .gt. 0) then
        l_matr_symm = .false.
    end if

! - Properties of behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    defoComp = compor(DEFO)
    typeComp = compor(INCRELAS)
    leul = defoComp .eq. 'GROT_GDEP'
    read (compor(NVAR), '(I16)') nbVari
    ASSERT((.not. defoComp .eq. 'GROT_GDEP') .or. (.not. typeComp .eq. 'COMP_ELAS'))

! - Output fields
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
!
    lg_varip = npg*nbsp*nbVari
    allocate (varip(lg_varip))
    if (lVari) then
        call jevech('PVARIMP', 'L', ivarix)
        varip(1:lg_varip) = zr(ivarix:ivarix+lg_varip-1)
    end if

! - Preparation of displacements
    do ino = 1, nbNode
        um(1, ino) = uml(1, ino)
        um(2, ino) = uml(2, ino)
        uf(1, ino) = uml(3, ino)
        uf(2, ino) = uml(5, ino)
        uf(3, ino) = -uml(4, ino)
        dum(1, ino) = dul(1, ino)
        dum(2, ino) = dul(2, ino)
        duf(1, ino) = dul(3, ino)
        duf(2, ino) = dul(5, ino)
        duf(3, ino) = -dul(4, ino)
    end do

! - Get layers
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    ASSERT(nbcou .gt. 0)
    hic = h/nbcou
    zmin = -h/deux+distn

! - Coefficients for shear stresses
    zmax = distn+h/2.d0
    quotient = 1.d0*zmax**3-3*zmax**2*zmin+3*zmax*zmin**2-1.d0*zmin**3
    a = -6.d0/quotient
    b = +6.d0*(zmin+zmax)/quotient
    c = -6.d0*zmax*zmin/quotient

! - Hooke matrix for shear
    call dxmate(plateCara, plateOrie, &
                fami, dfel, dmel, dmfel, dcel, &
                dciel, dmcel, dfcel, &
                multicel, coupmfel)
    do ino = 1, nbNode
        depfel(1+3*(ino-1)) = uf(1, ino)+duf(1, ino)
        depfel(2+3*(ino-1)) = uf(2, ino)+duf(2, ino)
        depfel(3+3*(ino-1)) = uf(3, ino)+duf(3, ino)
    end do

! - Loop on Gauss points
    do kpg = 1, npg
        n = 0.d0
        m = 0.d0
        df = 0.d0
        dm = 0.d0
        dmf = 0.d0
        dvt = 0.d0
        vt = 0.d0

! ----- Current Gauss point
        qsi = zr(icoopg-1+ndim*(kpg-1)+1)
        eta = zr(icoopg-1+ndim*(kpg-1)+2)

! ----- Prepare quantities for strain operators
        if (lQuad) then
            call jquad4(xyzl, qsi, eta, jacob)
            poids = zr(ipoids+kpg-1)*jacob(1)
            call dxqbm(qsi, eta, jacob(2), bm)
            call dkqbf(qsi, eta, jacob(2), cara, bf)
            call dsxhft(dfel, jacob(2), hft2el)
            call dkqtxy(qsi, eta, hft2el, depfel, cara(13), &
                        cara(9), vt)
        else if (lTria) then
            poids = zr(ipoids+kpg-1)*cara(7)
            call dxtbm(cara(9), bm)
            call dktbf(qsi, eta, cara, bf)
            call dsxhft(dfel, cara(9), hft2el)
            call dkttxy(cara(16), cara(13), hft2el, depfel, vt)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Compute strain and curvature
        call pmrvec('ZERO', 3, 2*nbNode, bm, um, eps)
        call pmrvec('ZERO', 3, 2*nbNode, bm, dum, deps)
        call pmrvec('ZERO', 3, 3*nbNode, bf, uf, khi)
        call pmrvec('ZERO', 3, 3*nbNode, bf, duf, dkhi)

! ----- Quadratic terms for GROT_GDEP
        if (leul) then
            bmq = 0.d0
            do i = 1, 2
                do k = 1, nbNode
                    do j = 1, 2
                        bmq(i, j) = bmq(i, j)+bm(i, 2*(k-1)+i)*dum(j, k)
                    end do
                    bmq(i, 3) = bmq(i, 3)+bm(i, 2*(k-1)+i)*duf(1, k)
                end do
            end do
            do k = 1, 3
                do i = 1, 2
                    deps(i) = deps(i)-0.5d0*bmq(i, k)*bmq(i, k)
                end do
                deps(3) = deps(3)-bmq(1, k)*bmq(2, k)
            end do
        end if

! ----- Loop on layers to integrate behaviour
        do icou = 1, nbcou
            do igauh = 1, npgh
! ------------- Current (sub)-point
                ksp = (icou-1)*npgh+igauh
                isp = (icou-1)*npgh+igauh
                ivpg = ((kpg-1)*nbsp+isp-1)*nbVari
                icpg = ((kpg-1)*nbsp+isp-1)*nbcon

! ------------- Value of height for integration point
                if (igauh .eq. 1) then
                    zic = zmin+(icou-1)*hic
                    coef = 1.d0/3.d0
                else if (igauh .eq. 2) then
                    zic = zmin+hic/deux+(icou-1)*hic
                    coef = 4.d0/3.d0
                else
                    zic = zmin+hic+(icou-1)*hic
                    coef = 1.d0/3.d0
                end if

! ------------- Compute 2D strains
                eps2d(1) = eps(1)+zic*khi(1)
                eps2d(2) = eps(2)+zic*khi(2)
                eps2d(3) = 0.0d0
                eps2d(4) = (eps(3)+zic*khi(3))/rac2
                eps2d(5) = 0.d0
                eps2d(6) = 0.d0
                deps2d(1) = deps(1)+zic*dkhi(1)
                deps2d(2) = deps(2)+zic*dkhi(2)
                deps2d(3) = 0.0d0
                deps2d(4) = (deps(3)+zic*dkhi(3))/rac2
                deps2d(5) = 0.d0
                deps2d(6) = 0.d0

! ------------- Elastic matrix for shear stresses
                d1iel(1, 1) = a*zic*zic+b*zic+c
                d1iel(2, 2) = d1iel(1, 1)
                d1iel(1, 2) = 0.d0
                d1iel(2, 1) = 0.d0

! ------------- Prepare stresses
                do j = 1, 4
                    sigmPrep(j) = zr(icontm+icpg-1+j)
                end do
                sigmPrep(4) = sigmPrep(4)*rac2
                if (lSigm) then
                    do j = 1, 6
                        zr(icontp+icpg-1+j) = 0.d0
                    end do
                end if

! ------------- Set main parameters for behaviour (on point)
                call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ------------- Integrator
                call nmcomp(BEHInteg, &
                            ndimLdc, option, typmod, &
                            instm, instp, &
                            compor, zr(jvCarcri), multComp, &
                            4, eps2d, deps2d, &
                            4, sigmPrep, &
                            zr(ivarim+ivpg), &
                            zr(icontp+icpg), varip(1+ivpg), &
                            36, dsidep, &
                            codkpg)

                if (codkpg .ne. 0) then
                    if (codret .ne. 1) then
                        codret = codkpg
                    end if
                end if

! ------------- Get stresses
                if (lSigm) then
                    zr(icontp+icpg+3) = zr(icontp+icpg+3)/rac2
                    zr(icontp+icpg+4) = d1iel(1, 1)*vt(1)+d1iel(1, 2)*vt(2)
                    zr(icontp+icpg+5) = d1iel(2, 1)*vt(1)+d1iel(2, 2)*vt(2)
                end if

! ------------- Compute vector
                if (lVect) then
                    coehsd = coef*hic/deux
                    n(1) = n(1)+coehsd*zr(icontp+icpg-1+1)
                    n(2) = n(2)+coehsd*zr(icontp+icpg-1+2)
                    n(3) = n(3)+coehsd*zr(icontp+icpg-1+4)
                    m(1) = m(1)+coehsd*zic*zr(icontp+icpg-1+1)
                    m(2) = m(2)+coehsd*zic*zr(icontp+icpg-1+2)
                    m(3) = m(3)+coehsd*zic*zr(icontp+icpg-1+4)
                end if
! ------------- Compute matrix
                if (lMatr) then
                    d2d(1) = dsidep(1, 1)
                    d2d(2) = dsidep(1, 2)
                    d2d(3) = dsidep(1, 4)/rac2
                    d2d(4) = dsidep(2, 1)
                    d2d(5) = dsidep(2, 2)
                    d2d(6) = dsidep(2, 4)/rac2
                    d2d(7) = dsidep(4, 1)/rac2
                    d2d(8) = dsidep(4, 2)/rac2
                    d2d(9) = dsidep(4, 4)/deux
                    do k = 1, 9
                        dm(k) = dm(k)+coef*hic/deux*poids*d2d(k)
                        dmf(k) = dmf(k)+coef*hic/deux*poids*zic*d2d(k)
                        df(k) = df(k)+coef*hic/deux*poids*zic*zic*d2d(k)
                    end do
                end if
            end do
        end do
! ----- BTSIG = BTSIG + BFT*M + BMT*N
        if (lVect) then
            do ino = 1, nbNode
                do k = 1, 3
                    btsig(1, ino) = btsig(1, ino)+bm(k, 2*(ino-1)+1)*n(k)*poids
                    btsig(2, ino) = btsig(2, ino)+bm(k, 2*(ino-1)+2)*n(k)*poids
                    btsig(3, ino) = btsig(3, ino)+bf(k, 3*(ino-1)+1)*m(k)*poids
                    btsig(5, ino) = btsig(5, ino)+bf(k, 3*(ino-1)+2)*m(k)*poids
                    btsig(4, ino) = btsig(4, ino)-bf(k, 3*(ino-1)+3)*m(k)*poids
                end do
            end do
        end if
! ----- Elementary matrices (membrane, bending, ...)
        if (lMatr) then
            call utbtab('CUMU', 3, 2*nbNode, dm, bm, work, memb)
            call utbtab('CUMU', 3, 3*nbNode, df, bf, work, flex)
            call utctab('CUMU', 3, 3*nbNode, 2*nbNode, dmf, bf, bm, work, mefl)
        end if
    end do
!
! - Add elementary matrices in global matrix
!
    if (lMatr) then
        if (lTria) then
            if (l_matr_symm) then
                call dxtloc(flex, memb, mefl, ctor, ktan)
            else
                call dxtloc2(flex, memb, mefl, ctor, ktan)
            end if
        else if (lQuad) then
            if (l_matr_symm) then
                call dxqloc(flex, memb, mefl, ctor, ktan)
            else
                call dxqloc2(flex, memb, mefl, ctor, ktan)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Get internal variables
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        zr(ivarip:ivarip+lg_varip-1) = varip(1:lg_varip)
    end if
!
    deallocate (varip)
end subroutine
