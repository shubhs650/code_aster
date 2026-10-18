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

subroutine apmain(action, kptsc, rsolu, vcine, istop, &
                  iret)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
    use aster_petsc_module
    use petsc_data_module
    use saddle_point_data_module, only: clear_augmented_lagrangian_context, clear_saddle_point_context
    use saddle_point_module
    use lmp_module, only: lmp_update
!
    implicit none
!
    character(len=*) :: action
    integer(kind=8) :: kptsc
    real(kind=8) :: rsolu(*)
    character(len=19) :: vcine
    integer(kind=8) :: istop, iret
!--------------------------------------------------------------
!
! IN  : ACTION :
!     /'DETR_MAT': POUR DETRUIRE L'INSTANCE PETSC ASSOCIEE A UNE MATRICE
!     /'PRERES'  : POUR CONSTRUIRE LE PRECONDITIONNEUR
!                 (ATTENTION EN // LA CONSTRUCTION DE CERTAINS PC EST
!                  RETARDEE)
!     /'RESOUD'  : POUR RESOUDRE LE SYSTEME LINEAIRE
!
! IN  : KPTSC (I): INDICE DES INSTANCES PETSC DANS Ap,Kp
! I/O : RSOLU (R): EN ENTREE : VECTEUR SECOND MEMBRE (REEL)
!                  EN SORTIE : VECTEUR SOLUTION      (REEL)
!                 (SI ACTION=RESOUD)
! IN  : VCINE (K19): NOM DU CHAM_NO DE CHARGEMENT CINEMATIQUE
!                   (SI ACTION=RESOUD)
! IN  : ISTOP (I)  : COMPORTEMENT EN CAS D'ERREUR
! OUT : IRET  (I)  : CODE RETOUR
!---------------------------------------------------------------
#include "asterc/asmpi_comm.h"
#include "asterc/create_custom_ksp.h"
#include "asterc/matfpe.h"
#include "asterfort/ap2foi.h"
#include "asterfort/apalmc.h"
#include "asterfort/apalmd.h"
#include "asterfort/apalmh.h"
#include "asterfort/apksp.h"
#include "asterfort/apmamc.h"
#include "asterfort/apmamd.h"
#include "asterfort/apmamh.h"
#include "asterfort/appcpr.h"
#include "asterfort/appcrs.h"
#include "asterfort/apsolu.h"
#include "asterfort/apvsmb.h"
#include "asterfort/apvsmbh.h"
#include "asterfort/assert.h"
#include "asterfort/cpysol.h"
#include "asterfort/csmbgg.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/filter_rhs.h"
#include "asterfort/filter_smd.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrconl.h"
#include "asterfort/mtdscr.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "asterfort/isParallelMatrix.h"
!
#ifdef ASTER_HAVE_PETSC
!----------------------------------------------------------------
!
!     VARIABLES LOCALES
    integer(kind=8) :: ifm, niv, ierd, nmaxit, ptserr, iaux
    integer(kind=8) :: lmat, idvalc, icode
    integer(kind=8), dimension(:), pointer :: slvi => null()
    mpi_int :: mpicomm
!
    character(len=24) :: precon, algo
    character(len=24), dimension(:), pointer :: slvk => null()
    character(len=19) :: nomat, nosolv
    character(len=14) :: nonu
    character(len=3) :: matd
    character(len=1) :: rouc
!
    real(kind=8) :: divtol, resipc
    real(kind=8), dimension(:), pointer :: slvr => null()
    complex(kind=8) :: cbid
!
    aster_logical :: lmd, l_parallel_matrix, lap2foi
    aster_logical, parameter :: dbg = .false.
!
!----------------------------------------------------------------
!     Variables PETSc
!
    PetscInt :: its, maxits
    PetscErrorCode ::  ierr
    PetscInt :: low, high
    PetscReal :: rtol, atol, dtol
    Vec :: r
    PetscScalar :: ires, fres
    PetscScalar, pointer :: xx(:) => null()
    KSPConvergedReason :: indic
    Mat :: a
    KSP, pointer :: ksp => null()
    PC :: pc
!----------------------------------------------------------------
    cbid = dcmplx(0.d0, 0.d0)
    call jemarq()
!
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicomm)
!
!   -- ON DESACTIVE LA LEVEE D'EXCEPTION FPE DANS LES BIBLIOTHEQUES MATHEMATIQUES
    call matfpe(-1_8)
!
    call infniv(ifm, niv)
!
!     -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
!
    call exisd('MATR_ASSE', nomat, icode)
    if (icode == 0) then
!   si la matrice n'existe pas, on peut quand meme
!   vouloir la detruire, mais c'est la seule action
!   autorisee
        ASSERT(action == 'DETR_MAT')
    else
        if (action .ne. 'DETR_MAT') then
            call jeveuo(nosolv//'.SLVK', 'L', vk24=slvk)
            precon = slvk(2)
            algo = slvk(6)
            call dismoi('MATR_DISTRIBUEE', nomat, 'MATR_ASSE', repk=matd)
            lmd = matd .eq. 'OUI'
            l_parallel_matrix = isParallelMatrix(nomat)
            ASSERT(.not. (lmd .and. l_parallel_matrix))
        end if
!
    end if
!
!
    if (action .eq. 'PRERES') then
!     ----------------------------
!
!        1.1 CREATION ET PREALLOCATION DE LA MATRICE PETSc :
!        ---------------------------------------------------
!
        if (.not. l_parallel_matrix) then
            if (lmd) then
                call apalmd(kptsc)
            else
                call apalmc(kptsc)
            end if
        else
            call apalmh(kptsc)
        end if
!
!        1.2 COPIE DE LA MATRICE ASTER VERS LA MATRICE PETSc :
!        -----------------------------------------------------
!
        if (.not. l_parallel_matrix) then
            if (lmd) then
                call apmamd(kptsc)
            else
                call apmamc(kptsc)
            end if
        else
            call apmamh(kptsc)
        end if
!
!        1.3 ASSEMBLAGE DE LA MATRICE PETSc :
!        ------------------------------------
!
        call MatAssemblyBegin(ap(kptsc), MAT_FINAL_ASSEMBLY, ierr)
        ASSERT(ierr .eq. 0)
        call MatAssemblyEnd(ap(kptsc), MAT_FINAL_ASSEMBLY, ierr)
        ASSERT(ierr .eq. 0)

        ! if(dbg) then
        !     fres = 0.d0
        !     call MatNorm(ap(kptsc), NORM_FROBENIUS, fres, ierr)
        !     ASSERT( ierr == 0 )
        !     call MatGetSize( ap(kptsc), mm, nn, ierr )
        !     ASSERT( ierr == 0 )
        !     print*, "SIZE LHS PETSC: ", mm, nn
        !     print*, "NORME LHS PETSC: ", fres
        ! end if
!
        if (precon == 'BLOC_LAGR') then
            call convert_mat_to_saddle_point(nomat, ap(kptsc))
        end if
!
!        1.4 CREATION DU PRECONDITIONNEUR PETSc (EXTRAIT DU KSP) :
!        ---------------------------------------------------------
!
        ksp => kp(kptsc)
        if (precon == 'UTILISATEUR') then
            call create_custom_ksp(ksp, ap(kptsc), ierr)
            ASSERT(ierr .eq. 0)
            user_ksp(kptsc) = ASTER_TRUE

            if (niv >= 2) then
                call KSPView(ksp, PETSC_VIEWER_STDOUT_SELF, ierr)
                ASSERT(ierr .eq. 0)
            end if
        else
            if (PetscObjectIsNull(ksp)) then
                call KSPCreate(mpicomm, ksp, ierr)
                ASSERT(ierr .eq. 0)
            end if
            !
            call KSPSetOperators(ksp, ap(kptsc), ap(kptsc), ierr)
            ASSERT(ierr == 0)
        end if
        !
        !  Initialisation du préconditionneur
        !
        call appcpr(kptsc)
!
    else if (action .eq. 'RESOUD') then
!     ---------------------------------
!
!        2.0 RECUPERATION DES POINTEURS DANS LE COMMUN :
!        -----------------------------------------------
!
        a = ap(kptsc)
        ksp => kp(kptsc)
!
!        2.1 PRETRAITEMENT DU SECOND MEMBRE :
!        ------------------------------------

!        -- DANS LE CAS DES LAGRANGE EN HPC, CERTAINES CONTRIBUTIONS
!           AU SECOND MEMBRE PEUVENT VENIR D'UN AUTRE PROCESSEUR
!           CAR LES MAILLES TARDIVES NE SONT PAS PARTAGEES DONC
!           CALCULEES SEULEMENT SUR UN PROC
        call filter_rhs(rsolu, nonu//".NUME")
!
!        -- MISE A L'ECHELLE DES LAGRANGES DANS LE SECOND MEMBRE
        call mtdscr(nomat)
        call jeveuo(nomat//'.&INT', 'L', lmat)
        call mrconl('MULT', lmat, 0_8, 'R', rsolu, 1_8)
!
!        -- MISE A ZERO DES TERMES NON CINEMATIQUES DONT LE PROC
!           COURANT N'EST PAS SEUL PROPRIETAIRE
        call filter_smd(nomat, rsolu)
!        -- PRISE EN COMPTE DES CHARGES CINEMATIQUES :
        call jeexin(vcine//'.VALE', ierd)
        if (ierd .ne. 0) then
            call jeveuo(vcine//'.VALE', 'L', idvalc)
            call jelira(vcine//'.VALE', 'TYPE', cval=rouc)
            ASSERT(rouc .eq. 'R')
            call csmbgg(lmat, rsolu, zr(idvalc), [cbid], [cbid], &
                        'R')
        end if
!
!        2.2 CREATION DU VECTEUR SECOND MEMBRE PETSc :
!        ---------------------------------------------
!
        if (.not. l_parallel_matrix) then
            call apvsmb(kptsc, lmd, rsolu)
        else
            call apvsmbh(kptsc, rsolu)
        end if

!        2.3 PARAMETRES DU KSP :
!        -----------------------
!
        call apksp(kptsc)
!
!        2.3b CREATION DES PRECONDITIONNEURS RETARDES :
!        ----------------------------------------------
!
        call appcrs(kptsc, lmd)
!
!        2.4 RESOLUTION :
!        ----------------
!
        call VecDuplicate(b, x, ierr)
        ASSERT(ierr .eq. 0)
!
        if (dbg) then
            fres = 0.d0
            call VecNorm(b, norm_2, fres, ierr)
            ASSERT(ierr == 0)
            print *, "NORME RHS PETSC: ", fres
        end if
!
        call KSPSolve(ksp, b, x, ierr)

        if (dbg) then
            fres = 0.d0
            call VecNorm(x, norm_2, fres, ierr)
            ASSERT(ierr == 0)
            print *, "NORME SOL PETSC: ", fres
        end if
!
!        2.5 DIAGNOSTIC :
!        ----------------
!
!       ARRET ANORMAL DU KSP
        if (ierr .gt. 0) call utmess('F', 'PETSC_13')
!
!       ANALYSE DE LA CONVERGENCE DU KSP
        call KSPGetConvergedReason(ksp, indic, ierr)
        ASSERT(ierr .eq. 0)
        call KSPGetIterationNumber(ksp, its, ierr)
        ASSERT(ierr .eq. 0)

!
!       -- si LDLT_SP/DP et its > maxits, on essaye une 2eme fois
!       -- apres avoir actualise le preconditionneur :
        if ((indic .eq. KSP_DIVERGED_ITS) .and. ((precon .eq. 'LDLT_SP') &
                                                 .or. (precon .eq. 'LDLT_DP'))) then
            call ap2foi(kptsc, mpicomm, nosolv, lmd, indic, its)
!           -- ksp a ete modifie par ap2foi :
            ksp => kp(kptsc)
            lap2foi = .true.
        else
            lap2foi = .false.
        end if
!
!
!       ANALYSE DES CAUSES ET EMISSION EVENTUELLE D'UN MESSAGE
!       EN CAS DE DIVERGENCE
        if (indic%v .lt. 0) then
            call KSPGetTolerances(ksp, rtol, atol, dtol, maxits, &
                                  ierr)
            ASSERT(ierr .eq. 0)
!

            if (indic .eq. KSP_DIVERGED_ITS) then
!               -- NOMBRE MAX D'ITERATIONS
                if (istop == 0) then
!                  ERREUR <F>
                    nmaxit = maxits
                    call utmess('F', 'PETSC_5', si=nmaxit)
                else if (istop == 2) then
!                  ON CONTINUE ET ON REMONTE UN CODE D'ERREUR
                    iret = 1
                    goto 999
                else
                    ASSERT(.false.)
                end if
!
            else if (indic .eq. KSP_DIVERGED_DTOL) then
!               DIVERGENCE
                if (istop == 0) then
!                  ERREUR <F>
                    divtol = dtol
                    call utmess('F', 'PETSC_6', sr=divtol)
                else if (istop == 2) then
!                  ON CONTINUE ET ON REMONTE UN CODE D'ERREUR
                    iret = 1
                    goto 999
                else
                    ASSERT(.false.)
                end if
!
            else if (indic .eq. KSP_DIVERGED_BREAKDOWN) then
!               BREAKDOWN
                if (istop == 0) then
!                  ERREUR <F>
                    call utmess('F', 'PETSC_7')
                else if (istop == 2) then
!                  ON CONTINUE ET ON REMONTE UN CODE D'ERREUR
                    iret = 1
                    goto 999
                else
                    ASSERT(.false.)
                end if
!
            else if (indic .eq. KSP_DIVERGED_NONSYMMETRIC) then
!               MATRICE NON SYMETRIQUE
                call utmess('F', 'PETSC_9')

!
            else if (indic .eq. KSP_DIVERGED_INDEFINITE_PC) then
!              PRECONDITIONNEUR NON DEFINI
                call utmess('F', 'PETSC_10')
!
            else if (indic .eq. KSP_DIVERGED_NANORINF) then
!               NANORINF
                if (istop == 0) then
!                  ERREUR <F>
                    call utmess('F', 'PETSC_8')
                else if (istop == 2) then
!                  ON CONTINUE ET ON REMONTE UN CODE D'ERREUR
                    iret = 1
                    goto 999
                else
                    ASSERT(.false.)
                end if
!
            else if (indic .eq. KSP_DIVERGED_INDEFINITE_MAT) then
!               MATRICE NON DEFINIE
                call utmess('F', 'PETSC_11')
!
            else
!              AUTRE ERREUR
                ptserr = indic%v
                call utmess('F', 'PETSC_12', si=ptserr)
            end if
        end if
!
!        2.5b VERIFICATION DE LA SOLUTION :
!        ----------------------------------
!
!        -- DOIT-ON VERIFIER LE CRITERE EN NORME NON PRECONDITIONNEE ?
        call jeveuo(nosolv//'.SLVR', 'L', vr=slvr)
        resipc = slvr(5)
!
        if (resipc .ge. 0.d0) then
            call VecDuplicate(x, r, ierr)
            ASSERT(ierr .eq. 0)
!           r = Ax
            call MatMult(ap(kptsc), x, r, ierr)
            ASSERT(ierr .eq. 0)
!           r = b - Ax
            call VecAYPX(r, -1.d0, b, ierr)
            ASSERT(ierr .eq. 0)
!           fres = ||r||_2
            call VecNorm(r, norm_2, fres, ierr)
            ASSERT(ierr .eq. 0)
!           ires = ||b||_2
            call VecNorm(b, norm_2, ires, ierr)
            ASSERT(ierr .eq. 0)
!
            call VecDestroy(r, ierr)
            ASSERT(ierr .eq. 0)
!
            call KSPGetTolerances(ksp, rtol, atol, dtol, maxits, &
                                  ierr)
            ASSERT(ierr .eq. 0)
!
            if (fres .gt. sqrt(rtol)*ires) then
                call utmess('F', 'PETSC_16', sr=fres)
            end if
        end if
!
!        2.6 RECOPIE DE LA SOLUTION :
!        ----------------------------
        if (.not. l_parallel_matrix) then
            call apsolu(kptsc, lmd, rsolu)
        else
            call VecGetOwnershipRange(x, low, high, ierr)
            ASSERT(ierr .eq. 0)
!
!        -- RECOPIE DE DANS RSOLU
            call VecGetArray(x, xx, ierr)
            ASSERT(ierr .eq. 0)
!
            call cpysol(nomat, nonu, rsolu, low, xx)
!
            call VecRestoreArray(x, xx, ierr)
            ASSERT(ierr .eq. 0)
        end if

!        2.7 UTILISATION DU LMP EN 2ND NIVEAU
!        -------------------------------------

        call jeveuo(nosolv//'.SLVK', 'L', vk24=slvk)
        algo = slvk(6)
        if (algo == 'GMRES_LMP') then
            call KSPGetPC(ksp, pc, ierr)
            ASSERT(ierr == 0)
            call lmp_update(pc, ksp, ierr)
            ASSERT(ierr == 0)
        end if
!
!         2.8 NETTOYAGE PETSc (VECTEURS) :
!         --------------------------------
!
!        -- EN CAS D'ERREUR DANS LES ITERATIONS DE KRYLOV ON SAUTE ICI
999     continue
        call VecDestroy(b, ierr)
        ASSERT(ierr .eq. 0)
        call VecDestroy(x, ierr)
        ASSERT(ierr .eq. 0)
!
!        -- PRECONDITIONNEUR UTILISE
!
!        -- TRAITEMENT PARTICULIER DU PRECONDITIONNEUR LAGRANGIEN AUGMENTE
        if (precon .eq. 'BLOC_LAGR') then
!
!           ON STOCKE LE NOMBRE D'ITERATIONS DU KSP
            call KSPGetIterationNumber(ksp, maxits, ierr)
            ASSERT(ierr .eq. 0)
            nmaxit = maxits
            call jeveuo(nosolv//'.SLVI', 'E', vi=slvi)
            slvi(5) = nmaxit
        end if
!
!        -- TRAITEMENT PARTICULIER DU PRECONDITIONNEUR LDLT_SP et LDLT_DP
        if (precon .eq. 'LDLT_SP' .or. precon .eq. 'LDLT_DP') then
!           MENAGE
            spsomu = slvk(3) (1:19)
            call detrsd('SOLVEUR', spsomu)
            spsomu = ' '
!
            call VecDestroy(xlocal, ierr)
            ASSERT(ierr .eq. 0)
            PetscObjectNullify(xlocal)
!
            call VecDestroy(xglobal, ierr)
            ASSERT(ierr == 0)
            PetscObjectNullify(xglobal)
!
            call VecScatterDestroy(xscatt, ierr)
            ASSERT(ierr == 0)
            PetscObjectNullify(xscatt)
!           ON STOCKE LE NOMBRE D'ITERATIONS DU KSP
            call KSPGetIterationNumber(ksp, maxits, ierr)
            ASSERT(ierr .eq. 0)
            iaux = maxits
            if (lap2foi) then
                if (niv .ge. 2) call utmess('I', 'PETSC_21', si=iaux)
!             ON REMET A ZERO POUR RELANCER LE CALCUL DU PRECONDITIONNEUR
!             SUR LA NOUVELLE MATRICE ET AVEC LE PARAMETRAGE INITIAL AU
!             PROCHAIN PAS DE NEWTON
                maxits = 0
            else
                if (niv .ge. 2) call utmess('I', 'PETSC_22', si=iaux)
            end if
            nmaxit = maxits
            call jeveuo(nosolv//'.SLVI', 'E', vi=slvi)
            slvi(5) = nmaxit
        end if
!
    else if (action .eq. 'DETR_MAT') then
!     -----------------------------------
!
!        3.0 RECUPERATION DES POINTEURS :
!        --------------------------------
!
        a = ap(kptsc)
        ksp => kp(kptsc)
!
!        3.1 NETTOYAGE PETSc :
!        ---------------------
!
!        -- DESTRUCTION DES OBJETS PETSC GENERAUX
        call MatDestroy(a, ierr)
        ASSERT(ierr .eq. 0)
!
!       user KSP will be removed on Python object deletion
        if (.not. user_ksp(kptsc)) then
            call KSPDestroy(ksp, ierr)
            ASSERT(ierr .eq. 0)
        end if
        call clear_augmented_lagrangian_context(nomat)
        call clear_saddle_point_context(nomat)
!
!        -- SUPRESSION DE L'INSTANCE PETSC
        nomats(kptsc) = ' '
        nosols(kptsc) = ' '
        nonus(kptsc) = ' '
        options(kptsc) = ' '
        PetscObjectNullify(ap(kptsc))
        PetscObjectNullify(kp(kptsc))
        tblocs(kptsc) = -1
!
!        -- PRECONDITIONNEUR UTILISE
!
        spmat = ' '
        spsolv = ' '
!
    else
        ASSERT(.false.)
    end if
!
!     -- ON REACTIVE LA LEVEE D'EXCEPTION
    call matfpe(1_8)
!
    call jedema()
!
#else
    character(len=1) :: kdummy
    integer(kind=8) :: idummy
    real(kind=8) :: rdummy
    kdummy = action(1:1)
    idummy = kptsc
    rdummy = rsolu(1)
    kdummy = vcine(1:1)
    idummy = istop
    idummy = iret
    ASSERT(.false.)
#endif
!
end subroutine
