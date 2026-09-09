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

subroutine nmmaba(jvMaterCode, relaComp, e, dsde, sigy)
!
    implicit none
!
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: relaComp
    integer(kind=8) :: jvMaterCode
    real(kind=8) :: e, dsde, sigy
!     RECUPERATION DES CARACTERISTIQUES DU MATERIAU POUR LES ELEMENTS
!     MECA_BARRE DANS LE CAS DES COMPORTEMENTS NON LINEAIRES.
! ----------------------------------------------------------------------
!
! IN  : ICODMA : NUMERO DU MATERIAU CODE
!       COMPOR : NOM DE LA RELATION DE COMPORTEMENT
!
! OUT : E      : MODULE D'YOUNG
!       DSDE   : PENTE D'ECROUISSAGE
!       SIGY   : LIMITE ELASTIQUE POUR L,ECROUISSAGE LINEAIRE
!       NCSTPM : NOMBRE DE CONSTANTES DE MATERIAU POUR INTO
!       CSTPM  : CONSTANTES DE MATERIAU :
!           E      : MODULE D'YOUNG
!           SY     : LIMITE ELASTIQUE
!           EPSU   : DEFORMATION ULTIME
!           SU     : CONTRAINTE ULTIME
!           EPSH   : DEFORMATION A LA FIN DU PALIER PLASTIQUE PARFAIT
!           R0     : COEFFICIENT EXPERIMENTAL
!           B      : COEFFICIENT
!           A1     : COEFFICIENT EXPERIMENTAL
!           A2     : COEFFICIENT EXPERIMENTAL
!           ELAN   : RAPPORT LONGUEUR/DIAMETRE DE LA BARRE
!           A6     : COEFFICIENT EXPERIMENTAL FLAMMBAGE
!           C      : COEFFICIENT EXPERIMENTAL FLAMMBAGE
!           COA    : COEFFICIENT EXPERIMENTAL FLAMMBAGE
!
!
! *************** DECLARATION DES VARIABLES LOCALES ********************
!
!-----------------------------------------------------------------------omecl
    integer(kind=8) :: nbpar, nbres, nbval
!-----------------------------------------------------------------------
    parameter(nbval=12)
    real(kind=8) :: valpar, valres(nbval)
    integer(kind=8) :: codres(nbval)
    character(len=8) :: nompar, nomela(1)
    character(len=16) :: nomecl(2), nompim(12)
    character(len=8) :: fami
!
! *********** FIN DES DECLARATIONS DES VARIABLES LOCALES ***************
!
! ****************************** DATA  *********************************
!
    data nomela/'E'/
    data nomecl/'D_SIGM_EPSI', 'SY'/
    data nompim/'SY', 'EPSI_ULTM', 'SIGM_ULTM', 'ELAN', 'EPSP_HARD',&
     &            'R_PM', 'EP_SUR_E', 'A1_PM', 'A2_PM', 'A6_PM',&
     &            'C_PM', 'A_PM'/
!
! ********************* DEBUT DE LA SUBROUTINE *************************
!
! --- MESSAGE D'ERREUR SI COMPORTEMENT NON REPERTORIE POUR LES BARRES
!
    if ((relaComp .ne. 'ELAS') .and. (relaComp .ne. 'VMIS_ISOT_LINE') .and. &
        (relaComp .ne. 'VMIS_CINE_LINE') .and. (relaComp .ne. 'VMIS_ASYM_LINE') .and. &
        (relaComp .ne. 'GRILLE_CINE_LINE') .and. (relaComp .ne. 'GRILLE_ISOT_LINE')) then
        call utmess('F', 'COMPOR4_32', sk=relaComp)
    end if
!
! --- INITIALISATIONS
!
    fami = 'RIGI'
!
    call r8inir(nbval, 0.d0, valres, 1)
    nbpar = 0
    nompar = '  '
    valpar = 0.d0
!
! --- CARACTERISTIQUES ELASTIQUES
!
    nbres = 2
    call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                ' ', 'ELAS', nbpar, nompar, [valpar], &
                1, nomela, valres, codres, 1)
    e = valres(1)
!
! --- CARACTERISTIQUES ECROUISSAGE LINEAIRE
!
    if ((relaComp .eq. 'VMIS_ISOT_LINE') .or. (relaComp .eq. 'VMIS_CINE_LINE') .or. &
        (relaComp .eq. 'GRILLE_CINE_LINE') .or. (relaComp .eq. 'GRILLE_ISOT_LINE')) then
        nbres = 2
!
!
        call r8inir(nbval, 0.d0, valres, 1)
        nbpar = 0
        nompar = '  '
        valpar = 0.d0
        call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                    ' ', 'ECRO_LINE', nbpar, nompar, [valpar], &
                    1, nomecl, valres, codres, 1)
        call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                    ' ', 'ECRO_LINE', nbpar, nompar, [valpar], &
                    1, nomecl(2), valres(2), codres(2), 0)
        if (codres(2) .ne. 0) valres(2) = 0.d0
        dsde = valres(1)
        sigy = valres(2)
    end if
!
end subroutine
