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

#include "aster.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef ASTER_HAVE_MPI
#include "mpi.h"
#ifdef OPEN_MPI
#include <dlfcn.h>
#endif
#endif

#include "aster_fort_mpi.h"
#include "aster_fort_utils.h"
#include "aster_mpi.h"
#include "aster_utils.h"

/*! Global object that store the entire tree */
static aster_comm_t aster_mpi_world;

/*! and a pointer to the current node */
static aster_comm_t *aster_mpi_current = NULL;

#ifdef ASTER_HAVE_MPI
static MPI_Errhandler errhdlr;
#endif

#ifdef ASTER_HAVE_PETSC
#include "petsc.h"
#endif

/*! This module defines functions:
 * - to manage the MPI communicators
 * - to properly interrupt a MPI execution.
 *
 * The communicators are managed in C. Fortran subroutines call these functions.
 * But all the communications are initiated from the Fortran subroutines.
 *
 * Communicators are store in fortran as Code_Aster mpi_int (== MPI_Fint).
 * They are converted to MPI_Comm with MPI_Comm_f2c((MPI_Fint)fortran_comm)
 *
 * Naming convention:
 *      aster_mpi_xxx : C functions and global variable
 *      asmpi_xxx : Fortran functions
 *
 * @todo this text should comment the module file, not the first next
 * subroutine.
 */
/*
 *   PUBLIC FUNCTIONS
 *
 */
void aster_mpi_init( const MPI_Fint init_comm_world ) {
    /*! MPI initialization */
    /* aster_world is COMM_WORLD for code_aster */
    MPI_Comm aster_world;
#ifdef ASTER_HAVE_MPI
    int isdone;
#ifdef OPEN_MPI
    void *handle = 0;
    int mode = RTLD_NOW | RTLD_GLOBAL;
    mode |= RTLD_NOLOAD;
    if ( !handle )
        handle = dlopen( "libmpi.so.15", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.14", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.13", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.12", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.11", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.10", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.1", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so.0", mode );
    if ( !handle )
        handle = dlopen( "libmpi.so", mode );
#endif
    printf( "checking MPI initialization...\n" );
    if ( init_comm_world == 0 ) {
        aster_world = MPI_COMM_WORLD;
        printf( "using COMM_WORLD.\n" );
    } else {
        printf( "using MPI communicator passed in argument.\n" );
        aster_world = MPI_Comm_f2c( init_comm_world );
    }
    // MPI should have been initialized during mpi4py import
    MPI_Initialized( &isdone );
    if ( !isdone ) {
        printf( "calling MPI_Init...\n" );
        AS_ASSERT( MPI_Init( 0, NULL ) == MPI_SUCCESS );
    } else {
        printf( "MPI is initialized.\n" );
    }
#ifdef ASTER_HAVE_PETSC
    // Ensure it has not been already initialized
    PetscBool initializedByPetsc4Py;
    PetscInitialized( &initializedByPetsc4Py );
    AS_ASSERT( initializedByPetsc4Py == PETSC_FALSE );
    // Pass communicator to PETSc
    PETSC_COMM_WORLD = aster_world;
#endif
    AS_ASSERT( atexit( terminate ) == 0 );
    /* set the error handler */
    AS_ASSERT( MPI_Comm_create_errhandler( errhdlr_func, &errhdlr ) == MPI_SUCCESS );
    AS_ASSERT( MPI_Comm_set_errhandler( aster_world, errhdlr ) == MPI_SUCCESS );
#else
    aster_world = 0;
#endif
    aster_mpi_world.id = aster_world;
    aster_mpi_world.parent = NULL;
    aster_mpi_world.level = 0;
    strncpy( aster_mpi_world.name, "WORLD", NAME_LENGTH );
    aster_mpi_current = &aster_mpi_world;
#ifdef ASTER_DEBUG_UNITTEST_MPI
    _unittest_aster_mpi();
#endif
    return;
}

/* API that works on aster_comm_t */
aster_comm_t *aster_get_comm_world() {
    /*! Return the original "COMM_APP" node */
    return &aster_mpi_world;
}

aster_comm_t *aster_get_current_comm() {
    /*! Return the current node */
    return aster_mpi_current;
}

void aster_set_current_comm( aster_comm_t *node ) {
    /*! Assign the current communicator */
    aster_mpi_current = node;
}

void aster_get_mpi_info( aster_comm_t *node, int *rank, int *size ) {
    /*! Return the rank of the process in `node` and its size
     * @param[in]   node    communicator
     * @param[out]  rank    rank of the current processor
     * @param[out]  size    number of processors in this communicator
     */
    *rank = 0;
    *size = 1;
    COMM_DEBUG( *node );
#ifdef ASTER_HAVE_MPI
    MPI_Comm_rank( node->id, rank );
    MPI_Comm_size( node->id, size );
#endif
    return;
}

aster_comm_t *aster_split_comm( aster_comm_t *parent, int color, int key, char *name ) {
    /*! Split the given communicator using color/key args,
     * return the sub-communicator */
    aster_comm_t *new;
#ifdef ASTER_HAVE_MPI
    MPI_Errhandler hdlr;
    int ierr;

    new = (aster_comm_t *)malloc( sizeof( aster_comm_t ) );
    ierr = MPI_Comm_split( parent->id, color, key, &( new->id ) );
    AS_ASSERT( ierr == MPI_SUCCESS );
    /* the parent has a new child */
    AS_ASSERT( parent->nbchild < MAX_CHILDS );
    parent->childs[parent->nbchild] = new;
    parent->nbchild++;
    /* fill the new node */
    new->parent = parent;
    new->level = parent->level + 1;
    new->nbchild = 0;
    strncpy( new->name, name, NAME_LENGTH );
    /* transfert the error handler - maybe already done by MPI_Comm_split */
    AS_ASSERT( MPI_Comm_get_errhandler( parent->id, &hdlr ) == MPI_SUCCESS );
    AS_ASSERT( MPI_Comm_set_errhandler( new->id, hdlr ) == MPI_SUCCESS );
    COMM_DEBUG( *new );
#else
    new = NULL;
#endif
    return new;
}

void aster_free_comm( aster_comm_t *node ) {
/*! delete this node */
#ifdef ASTER_HAVE_MPI
    aster_comm_t *parent;
    int i = 0;
    int nb, j, ierr;

    AS_ASSERT( node->nbchild == 0 );
    /* remove node from its parent childs list*/
    parent = node->parent;
    nb = parent->nbchild;
    while ( i < nb && parent->childs[i] != node ) {
        i++;
    }
    AS_ASSERT( i < nb );
    for ( j = i + 1; j < nb; j++ ) {
        parent->childs[j - 1] = parent->childs[j];
    }
    parent->childs[nb - 1] = NULL;
    parent->nbchild = nb - 1;
    /* delete the MPI_Comm */
    ierr = MPI_Comm_free( &( node->id ) );
    AS_ASSERT( ierr == MPI_SUCCESS );
    free( node );
#endif
    return;
}

/*
 * Wrapper around MPI_Barrier (because the communicator is optional in
 * asmpi_barrier)
 * Do not check returncode because all errors raise
 */
void DEFP( ASMPI_BARRIER_WRAP, asmpi_barrier_wrap, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    aster_comm_t *node;
    mpicom = MPI_Comm_f2c( *comm );
    node = get_node_by_id( &mpicom );
    aster_set_mpi_barrier( node );
#endif
    return;
}

int aster_set_mpi_barrier( aster_comm_t *node ) {
/*! Set a MPI barrier */
#ifdef ASTER_HAVE_MPI
    ASTERINTEGER iret, n0 = 0, n1 = 1, ibid = 0;
    ASTERDOUBLE rbid = 0.;
    char *valk;

    DEBUG_MPI( "mpi_barrier: %s is %d\n", "communicator", (int)MPI_Comm_c2f( node->id ) )
    CALL_ASMPI_CHECK( &iret );
    if ( iret != 0 ) {
        // valk = MakeCStrFromFStr("MPI_Barrier", VALK_SIZE);
        valk = MakeTabFStr( 1, VALK_SIZE );
        SetTabFStr( valk, 0, "MPI_Barrier", VALK_SIZE );
        CALL_UTMESS_CORE( "I", "APPELMPI_83", &n1, valk, &n0, &ibid, &n0, &rbid, &n0, " " );
        FreeStr( valk );
        return 1;
    }

    AS_MPICHECK( MPI_Barrier( node->id ) );
#endif
    return 0;
}

/* Tools allowing collective communications between processes */
int aster_mpi_bcast( void *buffer, int count, MPI_Datatype datatype, int root,
                     aster_comm_t *node ) {
/*! Broadcasts a message from one process to all other processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_Bcast: send %d values from proc #%d ...\n", count, root );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Bcast( buffer, count, datatype, root, node->id ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Bcast: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_allreduce( void *sendbuf, void *recvbuf, int count, MPI_Datatype sendtype, MPI_Op op,
                         aster_comm_t *node ) {
/*! Reduces a message and distributes the result to all other processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_Allreduce: send %d values to all %s...\n", count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allreduce( sendbuf, recvbuf, count, sendtype, op, node->id ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allreduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_gather( void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt,
                      MPI_Datatype recvtype, int root, aster_comm_t *node ) {
/*! Gathers together values from a group of processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_Gather: %d gathered values by proc #%d ...\n", sendcnt, root );
    double start = MPI_Wtime();
    AS_MPICHECK(
        MPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, node->id ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Gather: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_gatherv( void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
                       int *recvcnt, int *displ, MPI_Datatype recvtype, int root,
                       aster_comm_t *node ) {
/*! Gathers into specified locations from all processes in a group */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_Gatherv: %d gathered values by proc #%d ...\n", sendcnt, root );
    double start = MPI_Wtime();
    AS_ASSERT( MPI_Gatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, displ, recvtype, root,
                            node->id ) == MPI_SUCCESS );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Gatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_allgather( void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
                         int recvcnt, MPI_Datatype recvtype, aster_comm_t *node ) {
/*! Gathers together values from a group of processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_AllGather: %d gathered values by all %s...\n", sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK(
        MPI_Allgather( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, node->id ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_AllGather: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_allgatherv( void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
                          int *recvcnt, int *displs, MPI_Datatype recvtype, aster_comm_t *node ) {
/*! Gathers together values from a group of processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_AllGatherv: %d gathered values by all %s...\n", sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allgatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, displs, recvtype,
                                 node->id ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_AllGatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_send( void *sendbuf, int sendcnt, MPI_Datatype sendtype, int dest, int tag,
                    aster_comm_t *node ) {
/*! Gathers together values from a group of processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_Send: send %d values to %s...\n", sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Send( sendbuf, sendcnt, sendtype, dest, tag, node->id ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Send: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_recv( void *recvbuf, int recvcnt, MPI_Datatype recvtype, int source, int tag,
                    aster_comm_t *node ) {
/*! Gathers together values from a group of processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_Recv: receive %d values from %s...\n", recvcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Recv( recvbuf, recvcnt, recvtype, source, tag, node->id, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Recv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

int aster_mpi_sendrecv( void *sendbuf, int sendcnt, MPI_Datatype sendtype, int recv, int sendtag,
                        void *recvbuf, int recvcnt, MPI_Datatype recvtype, int source, int recvtag,
                        aster_comm_t *node ) {
/*! Gathers together values from a group of processes */
#ifdef ASTER_HAVE_MPI
    DEBUG_MPI( "MPI_SendRecv: receive %d values from %s...\n", recvcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Sendrecv( sendbuf, sendcnt, sendtype, recv, sendtag, recvbuf, recvcnt,
                               recvtype, source, recvtag, node->id, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_SendRecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return 0;
}

/* Access functions */
aster_comm_t *_search_id( aster_comm_t *node, MPI_Comm *id ) {
    /*! Search for 'id' in 'node' and its childs
     * Return NULL if not found.
     */
    aster_comm_t *found;
    int i;

    if ( node->id == *id ) {
        return node;
    } else {
        for ( i = 0; i < node->nbchild; i++ ) {
            found = _search_id( node->childs[i], id );
            if ( found ) {
                return found;
            }
        }
    }
    return NULL;
}

aster_comm_t *get_node_by_id( MPI_Comm *id ) {
    /*! Return the node that has the given 'id' */
    aster_comm_t *node;

    node = _search_id( &aster_mpi_world, id );
    AS_ASSERT( node );
    return node;
}

/*
 *  Fortran interfaces - wrappers of the C functions
 *
 */
void DEFSP( ASMPI_COMM, asmpi_comm, _IN const char *action, STRING_SIZE lact,
            _INOUT MPI_Fint *comm ) {
    /*! Wrapper around:
     *  aster_get_comm_world:   action = 'GET_WORLD', comm is OUT
     *  aster_get_current_comm: action = 'GET',       comm is OUT
     *  aster_set_current_comm: action = 'SET',       comm is IN
     *  aster_free_comm:        action = 'FREE',      comm is IN
     */
    MPI_Comm mpicom;
    aster_comm_t *node;
    char *act;

    act = MakeCStrFromFStr( action, lact );
    if ( strcmp( act, "GET_WORLD" ) == 0 ) {
        *comm = MPI_Comm_c2f( aster_get_comm_world()->id );
    } else if ( strcmp( act, "GET" ) == 0 ) {
        *comm = MPI_Comm_c2f( aster_get_current_comm()->id );
    } else if ( strcmp( act, "SET" ) == 0 ) {
        mpicom = MPI_Comm_f2c( *comm );
        node = get_node_by_id( &mpicom );
        aster_set_current_comm( node );
    } else if ( strcmp( act, "FREE" ) == 0 ) {
        mpicom = MPI_Comm_f2c( *comm );
        node = get_node_by_id( &mpicom );
        aster_free_comm( node );
    } else {
        AS_ASSERT( 0 );
    }
    FreeStr( act );
    return;
}

void DEFPPPSP( ASMPI_SPLIT_COMM, asmpi_split_comm, _IN MPI_Fint *parent, _IN MPI_Fint *color,
               MPI_Fint *key, _IN const char *name, STRING_SIZE lname, _OUT MPI_Fint *newcomm ) {
    /*! Wrapper around aster_split_comm */
    MPI_Comm mpicom;
    aster_comm_t *new;
    char *newname;

    newname = MakeCStrFromFStr( name, lname );
    mpicom = MPI_Comm_f2c( *parent );
    new = aster_split_comm( get_node_by_id( &mpicom ), (int)*color, (int)*key, newname );
    *newcomm = MPI_Comm_c2f( new->id );
    FreeStr( newname );
    return;
}

void DEFPPP( ASMPI_INFO_WRAP, asmpi_info_wrap, MPI_Fint *comm, MPI_Fint *rank, MPI_Fint *size ) {
    /*! Wrapper around aster_get_mpi_info
     * Called by the fortran subroutine asmpi_info where all arguments are
     * optional.
     */
    MPI_Comm mpicom;
    aster_comm_t *node;
    int irank = 0, isize = 1;

    AS_ASSERT( sizeof( MPI_Fint ) == 4 );

    mpicom = MPI_Comm_f2c( *comm );
    node = get_node_by_id( &mpicom );
    COMM_DEBUG( *node );
    aster_get_mpi_info( node, &irank, &isize );
    *rank = (MPI_Fint)irank;
    *size = (MPI_Fint)isize;
    return;
}

/*
 * Wrappers around MPI_Send
 * Do not check returncode because all errors raise
 */
void DEFPPPPP( ASMPI_SEND_R, asmpi_send_r, ASTERDOUBLE *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *dest, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Send: send %d double values to proc #%d ...\n", *count, *dest );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Send( (void *)buf, *count, MPI_DOUBLE_PRECISION, *dest, *tag, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Send: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_SEND_I, asmpi_send_i, ASTERINTEGER *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *dest, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Send: send %d integer8 values to proc #%d ...\n", *count, *dest );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Send( (void *)buf, *count, MPI_INTEGER8, *dest, *tag, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Send: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_SEND_I4, asmpi_send_i4, ASTERINTEGER4 *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *dest, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Send: send %d integer4 values to proc #%d ...\n", *count, *dest );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Send( (void *)buf, *count, MPI_INTEGER4, *dest, *tag, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Send:  ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_SEND_C, asmpi_send_c, ASTERDOUBLE *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *dest, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Send: send %d complex values to proc #%d ...\n", *count, *dest );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Send( (void *)buf, *count, MPI_DOUBLE_COMPLEX, *dest, *tag, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Send: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Recv
 * Do not check returncode because all errors raise
 */
void DEFPPPPP( ASMPI_RECV_R, asmpi_recv_r, ASTERDOUBLE *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *source, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Recv: recieve %d double values from proc #%d ...\n", *count, *source );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Recv( (void *)buf, *count, MPI_DOUBLE_PRECISION, *source, *tag, mpicom,
                           MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Recv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_RECV_I, asmpi_recv_i, ASTERINTEGER *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *source, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Recv: recieve %d integer8 values from proc #%d ...\n", *count, *source );
    double start = MPI_Wtime();
    AS_MPICHECK(
        MPI_Recv( (void *)buf, *count, MPI_INTEGER8, *source, *tag, mpicom, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Recv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_RECV_I4, asmpi_recv_i4, ASTERINTEGER4 *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *source, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Recv: recieve %d integer4 values from proc #%d ...\n", *count, *source );
    double start = MPI_Wtime();
    AS_MPICHECK(
        MPI_Recv( (void *)buf, *count, MPI_INTEGER4, *source, *tag, mpicom, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Recv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_RECV_C, asmpi_recv_c, ASTERDOUBLE *buf, ASTERINTEGER4 *count,
               ASTERINTEGER4 *source, ASTERINTEGER4 *tag, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Recv: recieve %d complex values from proc #%d ...\n", *count, *source );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Recv( (void *)buf, *count, MPI_DOUBLE_COMPLEX, *source, *tag, mpicom,
                           MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Recv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_SendRev
 * Do not check returncode because all errors raise
 */
void DEFPPPPPPPPP( ASMPI_SENDRECV_R, asmpi_sendrecv_r, ASTERDOUBLE *buffer_send,
                   ASTERINTEGER4 *count_send, ASTERINTEGER4 *recipient, ASTERINTEGER4 *tag_send,
                   ASTERDOUBLE *buffer_recv, ASTERINTEGER4 *count_recv, ASTERINTEGER4 *sender,
                   ASTERINTEGER4 *tag_recv, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_SendRecv: send/recv %d double values to proc #%d ...\n",
               *count_send + *count_recv, *recipient );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Sendrecv( (void *)buffer_send, *count_send, MPI_DOUBLE_PRECISION, *recipient,
                               *tag_send, (void *)buffer_recv, *count_recv, MPI_DOUBLE_PRECISION,
                               *sender, *tag_recv, mpicom, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Sendrecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPPPPP( ASMPI_SENDRECV_I, asmpi_sendrecv_i, ASTERINTEGER *buffer_send,
                   ASTERINTEGER4 *count_send, ASTERINTEGER4 *recipient, ASTERINTEGER4 *tag_send,
                   ASTERINTEGER *buffer_recv, ASTERINTEGER4 *count_recv, ASTERINTEGER4 *sender,
                   ASTERINTEGER4 *tag_recv, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_SendRecv: send/recv %d integer values to proc #%d ...\n",
               *count_send + *count_recv, *recipient );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Sendrecv( (void *)buffer_send, *count_send, MPI_INTEGER8, *recipient,
                               *tag_send, (void *)buffer_recv, *count_recv, MPI_INTEGER8, *sender,
                               *tag_recv, mpicom, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Sendrecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPPPPP( ASMPI_SENDRECV_I4, asmpi_sendrecv_i4, ASTERINTEGER4 *buffer_send,
                   ASTERINTEGER4 *count_send, ASTERINTEGER4 *recipient, ASTERINTEGER4 *tag_send,
                   ASTERINTEGER4 *buffer_recv, ASTERINTEGER4 *count_recv, ASTERINTEGER4 *sender,
                   ASTERINTEGER4 *tag_recv, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_SendRecv: send/recv %d short integer values to proc #%d ...\n",
               *count_send + *count_recv, *recipient );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Sendrecv( (void *)buffer_send, *count_send, MPI_INTEGER4, *recipient,
                               *tag_send, (void *)buffer_recv, *count_recv, MPI_INTEGER4, *sender,
                               *tag_recv, mpicom, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Sendrecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPPPPP( ASMPI_SENDRECV_C, asmpi_sendrecv_c, ASTERDOUBLE *buffer_send,
                   ASTERINTEGER4 *count_send, ASTERINTEGER4 *recipient, ASTERINTEGER4 *tag_send,
                   ASTERDOUBLE *buffer_recv, ASTERINTEGER4 *count_recv, ASTERINTEGER4 *sender,
                   ASTERINTEGER4 *tag_recv, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_SendRecv: send/recv %d complex values to proc #%d ...\n",
               *count_send + *count_recv, *recipient );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Sendrecv( (void *)buffer_send, *count_send, MPI_DOUBLE_COMPLEX, *recipient,
                               *tag_send, (void *)buffer_recv, *count_recv, MPI_DOUBLE_COMPLEX,
                               *sender, *tag_recv, mpicom, MPI_STATUS_IGNORE ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Sendrecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrapper around MPI_ISend
 * Do not check returncode because all errors raise
 */
void DEFPPPPPP( ASMPI_ISEND_I, asmpi_isend_i, ASTERINTEGER *buf, ASTERINTEGER4 *count,
                ASTERINTEGER4 *dest, ASTERINTEGER4 *tag, MPI_Fint *comm, MPI_Fint *request ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Request mpireq;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Isend: isend %d integer8 values to proc #%d ...\n", *count, *dest );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Isend( (void *)buf, *count, MPI_INTEGER8, *dest, *tag, mpicom, &mpireq ) );
    *request = MPI_Request_c2f( mpireq );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Isend: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrapper around MPI_IRecv
 * Do not check returncode because all errors raise
 */
void DEFPPPPPP( ASMPI_IRECV_I, asmpi_irecv_i, ASTERINTEGER *buf, ASTERINTEGER4 *count,
                ASTERINTEGER4 *source, ASTERINTEGER4 *tag, MPI_Fint *comm, MPI_Fint *request ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Request mpireq;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Irecv: irecieve %d integer8 values from proc #%d ...\n", *count, *source );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Irecv( (void *)buf, *count, MPI_INTEGER8, *source, *tag, mpicom, &mpireq ) );
    *request = MPI_Request_c2f( mpireq );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Irecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrapper around MPI_ISend
 * Do not check returncode because all errors raise
 */
void DEFPPPPPP( ASMPI_ISEND_I4, asmpi_isend_i4, ASTERDOUBLE *buf, ASTERINTEGER4 *count,
                ASTERINTEGER4 *dest, ASTERINTEGER4 *tag, MPI_Fint *comm, MPI_Fint *request ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Request mpireq;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Isend: isend %d integer4 values to proc #%d ...\n", *count, *dest );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Isend( (void *)buf, *count, MPI_INTEGER4, *dest, *tag, mpicom, &mpireq ) );
    *request = MPI_Request_c2f( mpireq );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Isend: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrapper around MPI_IRecv
 * Do not check returncode because all errors raise
 */
void DEFPPPPPP( ASMPI_IRECV_I4, asmpi_irecv_i4, ASTERDOUBLE *buf, ASTERINTEGER4 *count,
                ASTERINTEGER4 *source, ASTERINTEGER4 *tag, MPI_Fint *comm, MPI_Fint *request ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Request mpireq;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Irecv: irecieve %d integer4 values from proc #%d ...\n", *count, *source );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Irecv( (void *)buf, *count, MPI_INTEGER4, *source, *tag, mpicom, &mpireq ) );
    *request = MPI_Request_c2f( mpireq );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Irecv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*!
 * Wrapper around MPI_Test
 * Do not check returncode because all errors raise
 */
void DEFPP( ASMPI_TEST, asmpi_test, MPI_Fint *request, ASTERINTEGER4 *flag ) {
#ifdef ASTER_HAVE_MPI
    MPI_Request mpireq;
    int iflag;

    mpireq = MPI_Request_f2c( *request );
    DEBUG_MPI( "MPI_Test: Test %d communication request %s \n", mpireq, " " );
    AS_MPICHECK( MPI_Test( &mpireq, &iflag, MPI_STATUS_IGNORE ) );
    /* true=1, false=0 */
    *flag = (ASTERINTEGER4)iflag;
#endif
    return;
}

/*!
 * Wrapper around MPI_Cancel
 * Do not check returncode because all errors raise
 */
void DEFP( ASMPI_CANCEL, asmpi_cancel, MPI_Fint *request ) {
#ifdef ASTER_HAVE_MPI
    MPI_Request mpireq;

    mpireq = MPI_Request_f2c( *request );
    DEBUG_MPI( "MPI_Cancel: cancel %d communication request %s \n", mpireq, " " );
    AS_MPICHECK( MPI_Cancel( &mpireq ) );
#endif
    return;
}

/*!
 * Wrapper around MPI_Wtime
 * Do not check returncode because all errors raise
 */
ASTERDOUBLE DEF0( ASMPI_WTIME, asmpi_wtime ) {
#ifdef ASTER_HAVE_MPI
    return (ASTERDOUBLE)MPI_Wtime();
#else
    return (ASTERDOUBLE)0.0;
#endif
}

/*!
 * Wrappers around MPI_Reduce
 * Do not check returncode because all errors raise
 */
void DEFPPPPPP( ASMPI_REDUCE_R, asmpi_reduce_r, ASTERDOUBLE *sendbuf, ASTERDOUBLE *recvbuf,
                ASTERINTEGER4 *count, MPI_Fint *op, ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Reduce: %d double reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Reduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_DOUBLE_PRECISION, mpiop,
                             *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Reduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPP( ASMPI_REDUCE_C, asmpi_reduce_c, ASTERDOUBLE *sendbuf, ASTERDOUBLE *recvbuf,
                ASTERINTEGER4 *count, MPI_Fint *op, ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Reduce: %d double complex reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Reduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_DOUBLE_COMPLEX, mpiop,
                             *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Reduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPP( ASMPI_REDUCE_I, asmpi_reduce_i, ASTERINTEGER *sendbuf, ASTERINTEGER *recvbuf,
                ASTERINTEGER4 *count, MPI_Fint *op, ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Reduce: %d integer8 reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Reduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_INTEGER8, mpiop, *root,
                             mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Reduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPP( ASMPI_REDUCE_I4, asmpi_reduce_i4, ASTERINTEGER4 *sendbuf, ASTERINTEGER4 *recvbuf,
                ASTERINTEGER4 *count, MPI_Fint *op, ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Reduce: %d integer4 reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Reduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_INTEGER4, mpiop, *root,
                             mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Reduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Allreduce
 * Do not check returncode because all errors raise
 */
void DEFPPPPP( ASMPI_ALLREDUCE_R, asmpi_allreduce_r, ASTERDOUBLE *sendbuf, ASTERDOUBLE *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Allreduce: %d double reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allreduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_DOUBLE_PRECISION,
                                mpiop, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allreduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_ALLREDUCE_C, asmpi_allreduce_c, ASTERDOUBLE *sendbuf, ASTERDOUBLE *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Allreduce: %d double complex reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allreduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_DOUBLE_COMPLEX, mpiop,
                                mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allreduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_ALLREDUCE_I, asmpi_allreduce_i, ASTERINTEGER *sendbuf, ASTERINTEGER *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Allreduce: %d integer8 reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK(
        MPI_Allreduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_INTEGER8, mpiop, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allreduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPP( ASMPI_ALLREDUCE_I4, asmpi_allreduce_i4, ASTERINTEGER4 *sendbuf,
               ASTERINTEGER4 *recvbuf, ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    DEBUG_MPI( "MPI_Allreduce: %d integer4 reduced values by all ...%s\n", *count, " " );
    double start = MPI_Wtime();
    AS_MPICHECK(
        MPI_Allreduce( (void *)sendbuf, (void *)recvbuf, *count, MPI_INTEGER4, mpiop, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allreduce: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Bcast
 * Do not check returncode because all errors raise
 */
void DEFPPPP( ASMPI_BCAST_R, asmpi_bcast_r, ASTERDOUBLE *buffer, ASTERINTEGER4 *count,
              ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Bcast: send %d double values from proc #%d ...\n", *count, *root );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Bcast( (void *)buffer, *count, MPI_DOUBLE_PRECISION, *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Bcast: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPP( ASMPI_BCAST_C, asmpi_bcast_c, ASTERDOUBLE *buffer, ASTERINTEGER4 *count,
              ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Bcast: send %d double complex values from proc #%d ...\n", *count, *root );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Bcast( (void *)buffer, *count, MPI_DOUBLE_COMPLEX, *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Bcast: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPP( ASMPI_BCAST_I, asmpi_bcast_i, ASTERINTEGER *buffer, ASTERINTEGER4 *count,
              ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Bcast: send %d integer8 values from proc #%d ...\n", *count, *root );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Bcast( (void *)buffer, *count, MPI_INTEGER8, *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Bcast: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPP( ASMPI_BCAST_I4, asmpi_bcast_i4, ASTERINTEGER4 *buffer, ASTERINTEGER4 *count,
              ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Bcast: send %d integer4 values from proc #%d ...\n", *count, *root );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Bcast( (void *)buffer, *count, MPI_INTEGER4, *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Bcast: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFSPPP( ASMPI_BCAST_CHAR80, asmpi_bcast_char80, char *buffer, STRING_SIZE lbuff,
              ASTERINTEGER4 *count, ASTERINTEGER4 *root, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Bcast: send %d char80 values from proc #%d ...\n", *count, *root );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Bcast( (void *)buffer, ( *count ) * 80, MPI_CHAR, *root, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Bcast: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Allgather
 * Do not check returncode because all errors raise
 */
void DEFPPPPP( ASMPI_ALLGATHER_I, asmpi_allgather_i, ASTERINTEGER *sendbuf, ASTERINTEGER4 *sendcnt,
               ASTERINTEGER *recvbuf, ASTERINTEGER4 *recvcnt, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgather: %d gather integer values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allgather( (void *)sendbuf, *sendcnt, MPI_INTEGER8, (void *)recvbuf, *recvcnt,
                                MPI_INTEGER8, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgather: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Allgather
 * Do not check returncode because all errors raise
 */
void DEFPPPPP( ASMPI_ALLGATHER_R, asmpi_allgather_r, ASTERDOUBLE *sendbuf, ASTERINTEGER4 *sendcnt,
               ASTERDOUBLE *recvbuf, ASTERINTEGER4 *recvcnt, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgather: %d gather integer values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allgather( (void *)sendbuf, *sendcnt, MPI_DOUBLE_PRECISION, (void *)recvbuf,
                                *recvcnt, MPI_DOUBLE_PRECISION, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgather: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Allgather
 * Do not check returncode because all errors raise
 */
void DEFSPSPP( ASMPI_ALLGATHER_CHAR8, asmpi_allgather_char8, char *sendbuf, STRING_SIZE sbuff,
               ASTERINTEGER4 *sendcnt, char *recvbuf, STRING_SIZE rbuff, ASTERINTEGER4 *recvcnt,
               MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgather: %d gather integer values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allgather( (void *)sendbuf, ( *sendcnt ) * 8, MPI_CHAR, (void *)recvbuf,
                                ( *recvcnt ) * 8, MPI_CHAR, mpicom ) );
    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgather: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Allgatherv
 * Do not check returncode because all errors raise
 */
void DEFPPPPPP( ASMPI_ALLGATHERV_I, asmpi_allgatherv_i, ASTERINTEGER *sendbuf,
                ASTERINTEGER4 *sendcnt, ASTERINTEGER *recvbuf, ASTERINTEGER4 *recvcnt,
                ASTERINTEGER4 *displs, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgatherv: %d gather integer values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allgatherv( (void *)sendbuf, *sendcnt, MPI_INTEGER8, (void *)recvbuf, recvcnt,
                                 displs, MPI_INTEGER8, mpicom ) );

    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFPPPPPP( ASMPI_ALLGATHERV_R, asmpi_allgatherv_r, ASTERDOUBLE *sendbuf,
                ASTERINTEGER4 *sendcnt, ASTERDOUBLE *recvbuf, ASTERINTEGER4 *recvcnt,
                ASTERINTEGER4 *displs, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgatherv: %d gather double values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    AS_MPICHECK( MPI_Allgatherv( (void *)sendbuf, *sendcnt, MPI_DOUBLE_PRECISION, (void *)recvbuf,
                                 recvcnt, displs, MPI_DOUBLE_PRECISION, mpicom ) );

    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFSPSPPP( ASMPI_ALLGATHERV_CHAR16, asmpi_allgatherv_char16, char *sendbuf, STRING_SIZE sbuff,
                ASTERINTEGER4 *sendcnt, char *recvbuf, STRING_SIZE rbuff, ASTERINTEGER4 *recvcnt,
                ASTERINTEGER4 *displs, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgatherv: %d gather char16 values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    // We have to change size because of size of K16
    ASTERINTEGER4 size = 0;
    AS_MPICHECK( MPI_Comm_size( mpicom, &size ) );
    ASTERINTEGER4 *recvcnt_k16, *displs_k16;
    recvcnt_k16 = (ASTERINTEGER4 *)malloc( sizeof( ASTERINTEGER4 ) * size );
    displs_k16 = (ASTERINTEGER4 *)malloc( sizeof( ASTERINTEGER4 ) * size );

    ASTERINTEGER4 i;
    for ( i = 0; i < size; i++ ) {
        recvcnt_k16[i] = 16 * recvcnt[i];
        displs_k16[i] = 16 * displs[i];
    }

    AS_MPICHECK( MPI_Allgatherv( (void *)sendbuf, ( *sendcnt ) * 16, MPI_CHAR, (void *)recvbuf,
                                 recvcnt_k16, displs_k16, MPI_CHAR, mpicom ) );

    free( recvcnt_k16 );
    free( displs_k16 );

    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFSPSPPP( ASMPI_ALLGATHERV_CHAR24, asmpi_allgatherv_char24, char *sendbuf, STRING_SIZE sbuff,
                ASTERINTEGER4 *sendcnt, char *recvbuf, STRING_SIZE rbuff, ASTERINTEGER4 *recvcnt,
                ASTERINTEGER4 *displs, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgatherv: %d gather char24 values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    // We have to change size because of size of K24
    ASTERINTEGER4 size = 0;
    AS_MPICHECK( MPI_Comm_size( mpicom, &size ) );
    ASTERINTEGER4 *recvcnt_k24, *displs_k24;
    recvcnt_k24 = (ASTERINTEGER4 *)malloc( sizeof( ASTERINTEGER4 ) * size );
    displs_k24 = (ASTERINTEGER4 *)malloc( sizeof( ASTERINTEGER4 ) * size );

    ASTERINTEGER4 i;
    for ( i = 0; i < size; i++ ) {
        recvcnt_k24[i] = 24 * recvcnt[i];
        displs_k24[i] = 24 * displs[i];
    }

    AS_MPICHECK( MPI_Allgatherv( (void *)sendbuf, ( *sendcnt ) * 24, MPI_CHAR, (void *)recvbuf,
                                 recvcnt_k24, displs_k24, MPI_CHAR, mpicom ) );

    free( recvcnt_k24 );
    free( displs_k24 );

    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

void DEFSPSPPP( ASMPI_ALLGATHERV_CHAR80, asmpi_allgatherv_char80, char *sendbuf, STRING_SIZE sbuff,
                ASTERINTEGER4 *sendcnt, char *recvbuf, STRING_SIZE rbuff, ASTERINTEGER4 *recvcnt,
                ASTERINTEGER4 *displs, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;

    mpicom = MPI_Comm_f2c( *comm );
    DEBUG_MPI( "MPI_Allgatherv: %d gather char80 values by all ...%s\n", *sendcnt, " " );
    double start = MPI_Wtime();
    // We have to change size because of size of K80
    ASTERINTEGER4 size = 0;
    AS_MPICHECK( MPI_Comm_size( mpicom, &size ) );
    ASTERINTEGER4 *recvcnt_k80, *displs_k80;
    recvcnt_k80 = (ASTERINTEGER4 *)malloc( sizeof( ASTERINTEGER4 ) * size );
    displs_k80 = (ASTERINTEGER4 *)malloc( sizeof( ASTERINTEGER4 ) * size );

    ASTERINTEGER4 i;
    for ( i = 0; i < size; i++ ) {
        recvcnt_k80[i] = 80 * recvcnt[i];
        displs_k80[i] = 80 * displs[i];
    }

    AS_MPICHECK( MPI_Allgatherv( (void *)sendbuf, ( *sendcnt ) * 80, MPI_CHAR, (void *)recvbuf,
                                 recvcnt_k80, displs_k80, MPI_CHAR, mpicom ) );

    free( recvcnt_k80 );
    free( displs_k80 );

    double end = MPI_Wtime();
    DEBUG_MPI( "MPI_Allgatherv: ... in %f sec %s\n", ( end - start ), " " );
#endif
    return;
}

/*
 * Wrappers around MPI_Scan
 * Do not check returncode because all errors raise
 */
void DEFPPPPP( ASMPI_SCAN_R, asmpi_scan_r, ASTERDOUBLE *sendbuf, ASTERDOUBLE *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    AS_MPICHECK(
        MPI_Scan( (void *)sendbuf, (void *)recvbuf, *count, MPI_DOUBLE_PRECISION, mpiop, mpicom ) );
#endif
    return;
}

void DEFPPPPP( ASMPI_SCAN_C, asmpi_scan_c, ASTERDOUBLE *sendbuf, ASTERDOUBLE *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    AS_MPICHECK(
        MPI_Scan( (void *)sendbuf, (void *)recvbuf, *count, MPI_DOUBLE_COMPLEX, mpiop, mpicom ) );
#endif
    return;
}

void DEFPPPPP( ASMPI_SCAN_I, asmpi_scan_i, ASTERINTEGER *sendbuf, ASTERINTEGER *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    AS_MPICHECK(
        MPI_Scan( (void *)sendbuf, (void *)recvbuf, *count, MPI_INTEGER8, mpiop, mpicom ) );
#endif
    return;
}

void DEFPPPPP( ASMPI_SCAN_I4, asmpi_scan_i4, ASTERINTEGER4 *sendbuf, ASTERINTEGER4 *recvbuf,
               ASTERINTEGER4 *count, MPI_Fint *op, MPI_Fint *comm ) {
#ifdef ASTER_HAVE_MPI
    MPI_Comm mpicom;
    MPI_Op mpiop;

    mpicom = MPI_Comm_f2c( *comm );
    mpiop = MPI_Op_f2c( *op );
    AS_MPICHECK(
        MPI_Scan( (void *)sendbuf, (void *)recvbuf, *count, MPI_INTEGER4, mpiop, mpicom ) );
#endif
    return;
}

/*
 * Define a dedicated function to abort a Code_Aster execution.
 */
int gErrFlg = 0;

void DEFP( ASABRT, asabrt, _IN ASTERINTEGER *iret ) {
    /*! \brief Define a dedicated function to abort a Code_Aster execution.
     *
     * Function to interrupt the execution.
     * - In a sequential version, it just calls abort().
     * - In a MPI execution, it set a global flag and calls `MPI_Abort`.
     *
     * The usage of atexit seems required in a Python interpreter
     * certainly because `sys.exit()` probably calls the `exit` system
     * function (so we can't add a `MPI_Finalize` call before exiting).
     * But if `MPI_Finalize` is executed after a `MPI_Abort` call, all
     * the processes are not interrupted.
     * That's why a global flag is used to by-pass `MPI_Finalize` in
     * case of error.
     * But the same problem appears if `MPI_Finalize` is called before
     * a `MPI_Abort`. That's why a call to asmpi_check has been added before
     * calling MPI_Finalize.
     *
     * to test MPI_Abort : http://www.netlib.org/blacs/blacs_errata.html
     */
    gErrFlg = 1;
#ifdef ASTER_HAVE_MPI
    printf( "calling MPI_Abort with errorcode %d...\n", (int)*iret );
    fflush( stdout );
    MPI_Abort( aster_mpi_world.id, (int)( *iret ) );
#else
    CALL_ABORTF();
#endif
    return;
}

void terminate( void ) {
    /*! Function registered using atexit() in main.
     */
    printf( "End of the Code_Aster execution\n" );
#ifdef ASTER_HAVE_MPI
    int isdone;
    ASTERINTEGER dummy = 0;
    if ( gErrFlg == 0 ) {
        /* see help of asabrt */
        printf( "Code_Aster MPI exits normally\n" );
        MPI_Finalized( &isdone );
        if ( !isdone ) {
            CALL_ASMPI_CHECK( &dummy );
            MPI_Errhandler_free( &errhdlr );
            MPI_Finalize();
        }
    } else {
        printf( "Code_Aster MPI exits with errors\n" );
    }
#endif
    printf( "Exited\n" );
    return;
}

/*
 *   PRIVATE FUNCTIONS
 *
 */
#ifdef ASTER_HAVE_MPI
void errhdlr_func( MPI_Comm *comm, int *err, ... ) {
    /*! Error handler for calls to MPI functions */
    char errstr[MPI_MAX_ERROR_STRING];
    int len;

    AS_ASSERT( *err != MPI_SUCCESS );
    MPI_Error_string( *err, errstr, &len );
    printf( "\n<F> MPI Error code %d:\n    %s\n\n", *err, errstr );
    fflush( stdout );
    CALL_UTMESS( "F", "APPELMPI_5" );
    return;
}
#endif

/* ASTER_DEBUG_UNITTEST_MPI */
#ifdef ASTER_DEBUG_UNITTEST_MPI
void _unittest_aster_mpi() {
    /*! unittest of the functions on aster_comm_t tree */
    int size, rank, npband, npsolv;
    int color;
    aster_comm_t *node, *world, *scband, *sccross, *scsolv;

    COMM_DEBUG( aster_mpi_world );
    aster_get_mpi_info( &aster_mpi_world, &rank, &size );

    world = aster_get_comm_world();
    node = aster_get_current_comm();
    AS_ASSERT( world == node );

    npband = size / 2;
    npsolv = size / 4;
    if ( npsolv < 1 ) {
        printf( "this test requires at least 4 procs, 8 to be relevant\n" );
        return;
    }
    color = rank < npband;

    fprintf( stderr, "band color : %d\n", color );
    scband = aster_split_comm( world, color, rank, "band" );
    aster_set_current_comm( scband );
    COMM_DEBUG( *aster_mpi_current );

    color = rank % npband == 0;
    fprintf( stderr, "cross color : %d\n", color );
    sccross = aster_split_comm( world, color, rank, "cross" );
    COMM_DEBUG( *sccross );

    color = rank % npsolv == 0;
    fprintf( stderr, "solv color : %d\n", color );
    scsolv = aster_split_comm( aster_get_current_comm(), color, rank, "solver" );

    aster_get_mpi_info( world, &rank, &size );
    fprintf( stderr, "%-8s: rank=%d  size=%d\n", world->name, rank, size );

    aster_get_mpi_info( scband, &rank, &size );
    fprintf( stderr, "%-8s: rank=%d  size=%d\n", scband->name, rank, size );

    aster_get_mpi_info( sccross, &rank, &size );
    fprintf( stderr, "%-8s: rank=%d  size=%d\n", sccross->name, rank, size );

    aster_get_mpi_info( scsolv, &rank, &size );
    fprintf( stderr, "%-8s: rank=%d  size=%d\n", scsolv->name, rank, size );

    aster_set_current_comm( world );
    return;
}
#endif
