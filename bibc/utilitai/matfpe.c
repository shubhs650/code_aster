/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

/*
   Cette fonction permet de désactiver temporairement la levée d'exceptions
   autour des appels des fonctions des librairies mathématiques blas et
   lapack.
      CALL MATFPE(-1) : on désactive la levée d'exceptions
      CALL MATFPE(1) : on active la levée d'exceptions

   Problème rencontré sur Linux ia64 avec MKL 8.0.
*/
#include "aster.h"

#include "aster_sigfpe.h"

#include <stdio.h>

#if defined ASTER_HAVE_SUPPORT_FPE
#include <float.h>
#include <signal.h>

#ifdef ASTER_PLATFORM_POSIX
#define _GNU_SOURCE 1
#include <fenv.h>
#endif

void hanfpe( int sig );
void aster_fe_invalid( const int );

static int compteur_fpe = 1;
static int compteur_invalid = 0;
#endif

void DEFP( MATFPE, matfpe, ASTERINTEGER *enable ) {
#if defined ASTER_HAVE_SUPPORT_FPE

    /* permet juste de vérifier où on en est si besoin ! */
    if ( *enable == 0 ) {
        printf( "#MATFPE var = %ld (compteur %d)\n", *enable, compteur_fpe );
        return;
    }
    compteur_fpe = compteur_fpe + *enable;
    if ( compteur_fpe < 1 ) {
#if defined ASTER_PLATFORM_MINGW
        _controlfp( _MCW_EM, _MCW_EM );
#else
        fedisableexcept( ASTER_SIGFPE );
#endif
        /* définition du handler : hanfpe appelle UTMFPE qui fait UTMESS('F') */
        signal( SIGFPE, hanfpe );
    } else if ( compteur_fpe >= 1 ) {
        /* avant de reactiver le controle des FPE, on abaisse les flags */
#if defined ASTER_PLATFORM_MINGW
        _clearfp();
        _controlfp( _EM_INVALID | _EM_UNDERFLOW | _EM_DENORMAL | _EM_INEXACT, _MCW_EM );
#else
        feclearexcept( ASTER_SIGFPE );
        feenableexcept( ASTER_SIGFPE );
#endif
        /* définition du handler : hanfpe appelle UTMFPE qui fait UTMESS('F') */
        signal( SIGFPE, hanfpe );
    }
    aster_fe_invalid( (int)*enable );

#endif
}

void aster_fe_invalid( const int enable ) {

#if defined ASTER_HAVE_SUPPORT_FPE && defined ASTER_PLATFORM_LINUX
    // FE_INVALID was not enabled on mingw

    /* permet juste de vérifier où on en est si besoin ! */
    if ( enable == 0 ) {
        printf( "#ASTER_FE_INVALID var = %ld (compteur %d)\n", enable, compteur_invalid );
        return;
    }
    compteur_invalid = compteur_invalid + enable;
    if ( compteur_invalid < 1 ) {
        fedisableexcept( FE_INVALID );
        /* définition du handler : hanfpe appelle UTMFPE qui fait UTMESS('F') */
        signal( SIGFPE, hanfpe );
    } else if ( compteur_invalid >= 1 ) {
        // printf( "<i> FE_INVALID exception enabled\n" );
        /* avant de reactiver le controle des FPE, on abaisse les flags */
        feclearexcept( FE_INVALID );
        feenableexcept( FE_INVALID );
        /* définition du handler : hanfpe appelle UTMFPE qui fait UTMESS('F') */
        signal( SIGFPE, hanfpe );
    }
#endif
}
