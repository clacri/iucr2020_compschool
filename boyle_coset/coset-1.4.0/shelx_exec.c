/* contains implementation for running SHELXL jobs to
 * test the twin laws derived from the 
 * Flack left coset decomposition program which uses
 * uses alogorithms outlined in Acta Cryst. (1987), A43,
 * 564-568, by H. D. Flack.
 *
 * Support for UNIX fork()/exec() and symbolic link semantics.
 * For other systems the ANSI C system() and copying the .hkl
 * file is used. Compile with
 * -DUSE_SYSTEM_FUNCTION if fork()/exec() are not available.
 *
 * Copyright (C) 2008, 2012, 2013 Paul D. Boyle
 *
 * Written by:
 *      Paul D. Boyle
 *      Department of Chemistry
 *      University of Western Ontario
 *      London, Ontario CANADA N6A 5B7
 *
 *      February 2008
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version
 * 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#ifdef NEED_C89_COMPATIBILITY
#undef _ISOC99_SOURCE
#else
#define _ISOC99_SOURCE
#endif

#ifndef USE_SYSTEM_FUNCTION
#define _XOPEN_SOURCE 500
#define FORKEXEC
#else 
#define USE_NONSTANDARD_FOPEN  /* for nonconforming C implementations of fopen() */
#endif

#ifdef FORKEXEC
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>


#include "sll.h"
#include "dynamic_sll.h"
#include "dupstr.h"
#include "shelx.h"    /* for get_basename */
#include "shelx_exec.h"

/* for compilers which don't have C99 functions available */
#ifdef NEED_C89_COMPATIBILITY
#include "c89_util.h"
#endif


#define HKLF_NAME_FORMAT "%s.hkl"


static char *build_command_buffer( char *exe, char *arg )
{
    int ret;
    size_t sz;
    char *buf;
  
#ifdef TRAP_SHELX_STDOUT
    sz = strlen(exe) + 2 * strlen(arg) + strlen(TRAP_SUFFIX) + 5;
#else
    sz = strlen(exe) + strlen(arg) + 2;
#endif
    errno = 0;
    buf = malloc( sz );
    if( NULL == buf )
        return NULL;

#ifdef TRAP_SHELX_STDOUT
    ret = snprintf( buf, sz, COMMAND_BUF_FORMAT, exe, arg, arg, TRAP_SUFFIX );
#else
    ret = snprintf( buf, sz, COMMAND_BUF_FORMAT, exe, arg );
#endif

    if( ret >= (int)sz ) {
        fprintf( stdout, "%s:%d: Command buffer for SHELXL command was truncated, expect unusual output or problems.\n",
                         __FILE__,__LINE__ );
    }

    return buf;

}

#ifdef FORKEXEC
static int create_symlink( const char *from, const char *to )
{
    int ret;

    errno = 0;
    ret = unlink( to );  /* we don't care if this fails, just be sure
                          * to clear old links before making new ones.
                          */
    ret = symlink( from, to );
    return ret;
}
#else /* for systems which don't provide symbolic links, copy the .hkl files */
static int create_symlink( const char *from, const char *to )
{   FILE *in, *out;
    int c;
#ifdef USE_NONSTANDARD_FOPEN
    const char *wmode = "wt";
    const char *rmode = "rt";
#else  /* use only ANSI C Standard flags */
    const char *wmode = "w";
    const char *rmode = "r";
#endif

    in = fopen( from, rmode );
    if( NULL == in )
        return -1;

    out = fopen( to, wmode );
    if( NULL == out ) {
        fclose( in );
        return -1;
    }
        
    while( EOF != (c = fgetc( in ) ) ) {
           if( ferror(in) ) {
               return -1;
           }
           if( ferror(out) ) {
               return -1;
           }
           fputc( c, out );
    }
    fclose( in );
    fclose( out );
    return 0;
}
#endif


char *real_hklf_filename( char *ins_filename )
{
    char *base,
         *hklf_name;
    size_t len;

/* take advantage of usual SHELX conventions between the names of 
 * the .ins file and the .hkl file to derive the name of the 
 * .hkl file.
 */
    base = get_basename( ins_filename, '.' );
    if( NULL == base )
        return NULL;

    len = strlen(base);
    len += 5; /* add some to make room for the ".hkl" */
    hklf_name = malloc( len );
    if( NULL == hklf_name ) {
        free( base );
        return NULL;
    }

    snprintf( hklf_name, len, HKLF_NAME_FORMAT, base );

    free( base );
    return hklf_name;
}

int create_hklf_file_links( const char *real_hklf, SLinkedList *link_names )
{
    char *hkl_name;
    int ret;

    while( -1 != sll_remove_next( link_names, NULL, (void **)&hkl_name ) ) {
         errno = 0;
         ret = create_symlink( real_hklf, hkl_name );
         if( -1 == ret ) {
             fprintf( stderr, "Could not create symbolic link: %s: %s\n",
                      hkl_name, 0 != errno ? strerror(errno) : "error not cataloged by errno" );
         }
         free( hkl_name );
    }

    return ret;
}


SLinkedList *make_job_name_list( SLinkedList *ins_file_names )
{
    SLinkedList *jobs = NULL;
    char *base, 
         *tmp;

    errno = 0;
    jobs = alloc_list_init( jobs, free );
    if( NULL == jobs ) {
        return NULL;
    }

    while( -1 != sll_remove_next(ins_file_names, NULL, (void **)&tmp ) ) {
           errno = 0;
           base = get_basename( tmp, '.' );
           if( NULL == base ) {
               dealloc_list( jobs );
               return NULL;
           }
           free( tmp );
           sll_insert_next(jobs, sll_list_tail(jobs), base );
    }

/* ins_file_names should be empty now, but leave it to caller to dealloc_list() the list */

    return jobs;
}

SLinkedList *create_hkl_filenames( SLinkedList *job_names )
{
    SLinkedList *hkl_names = NULL;
    SLinkedListElem *el;
    size_t len;
    int ret;
    char *hkl_name;

    hkl_names = alloc_list_init( hkl_names, free );
    if( NULL == hkl_names ) {
        return NULL;
    }

    for( el = sll_list_head(job_names); NULL != el; el = sll_list_next(el) ) {
         len = strlen(sll_list_data(el) ) + 5;
         errno = 0;
         hkl_name = malloc( len );
         if( NULL == hkl_name ) {
             dealloc_list( hkl_names );
             return NULL;
         }
         ret = snprintf( hkl_name, len, HKLF_NAME_FORMAT, (char *)sll_list_data(el) );
         if( ret >= (int)len ) {
             fprintf( stderr, "%s:%d:  WARNING: snprintf() reports %d bytes written to buffer of length %u\n",
                               __FILE__,__LINE__, ret, (unsigned)len ); 
         }
         sll_insert_next( hkl_names, sll_list_tail(hkl_names), hkl_name );
    } 

    return hkl_names;
}

SLinkedList *setup_shelx_jobs( SLinkedList *ins_file_names, char *orig_ins_filename )
{
    SLinkedList *job_list = NULL,
                *hklf_file_list = NULL;
    char *real_hkl_filename = NULL;
    int ret;

    errno = 0;
    job_list = make_job_name_list( ins_file_names );
    if( NULL == job_list ) {
        fprintf( stderr, "%s:%d: make_job_name_list() returned NULL %s",
                 __FILE__,__LINE__, errno != 0 ? strerror(errno) : "unspecified error" );
        return NULL;
    }
    hklf_file_list = create_hkl_filenames( job_list );
    if( NULL == hklf_file_list ) {
        dealloc_list( job_list );
        return NULL;
    }
    real_hkl_filename = real_hklf_filename( orig_ins_filename );
    if( NULL == real_hkl_filename ) {
        dealloc_list( job_list );
        dealloc_list( hklf_file_list );
        return NULL;
    }
    ret = create_hklf_file_links( real_hkl_filename, hklf_file_list );
    if( -1 == ret ) {
        dealloc_list( job_list );
        dealloc_list( hklf_file_list );
        free( real_hkl_filename );
        return NULL;
    }

    dealloc_list( hklf_file_list );
    free( real_hkl_filename );
    return job_list;
}

#ifdef FORKEXEC /* use normal UNIX fork()/exec() semantics */
int spawn_shelx_jobs( SLinkedList *jobs, char *shelx_exe_path )
{
    pid_t pid;
    int status = 0,
        ret = 0,
        n_jobs = 0,
        i_jobs = 0;
    char *arg1,
         *cmd_buffer;

    while( -1 != sll_remove_next(jobs, NULL, (void **)&arg1 ) ) {
          switch( pid = fork() ) {
              case -1:  /* whoops, fork() croaked */
                 fprintf( stderr, "%s:%d: %s\n", __FILE__, __LINE__, strerror(errno) );
                 break;
              case  0:  /* child does this */
                 cmd_buffer = build_command_buffer( shelx_exe_path, arg1 );
                 if( NULL == cmd_buffer ) {
                     fprintf( stderr, "%s:%d: %s\n", __FILE__, __LINE__, 
                              errno != 0 ? strerror(errno) : "memory allocation error in build_command_buffer()" ); 
                 }
                 fprintf( stdout, "Executing trial refinement for %s ...\n", arg1 );
                 ret = execl("/bin/sh", "sh", "-c", cmd_buffer, (char *)0 );
                 if( -1 == ret ) {
                     fprintf( stderr, "executing %s %s failed: %s\n",
                              shelx_exe_path, arg1, strerror(errno) );
                 }
                 break;
              default:  /* parent does this */
                 n_jobs++;
                 free( arg1 );
                 break;
          }
    } 
    i_jobs = n_jobs;
    while( n_jobs > 0 ) {
        pid = wait( &status );
        if( WIFEXITED(status) ) {
            fprintf( stdout, "Process %d exited normally with exit code %d\n",
                      pid, WEXITSTATUS(status) );
            fflush(stdout);
        }
        else if(WIFSIGNALED(status)) {
            fprintf( stderr, "Process %d terminated abnormally. Caught signal: %d\n",
                     pid, WTERMSIG(status) );
            fflush(stderr);
        }
        n_jobs--;
    }

    return i_jobs;
}   
#else /* use ANSI C system() function */
int spawn_shelx_jobs( SLinkedList *jobs, char *shelx_exe_path )
{
    int status = 0,
        n_jobs = 0;

    char *arg1,
         *cmd_buffer;
      
    while( -1 != sll_remove_next(jobs, NULL, (void **)&arg1 ) ) {
          errno = 0;
          cmd_buffer = build_command_buffer( shelx_exe_path, arg1 );
          if( NULL == cmd_buffer ) {
              return -1;
          } 
          status = system( cmd_buffer );
          if( -1 != status ) {
              n_jobs++;
          }
          free( arg1 );
          free( cmd_buffer );
    }
    return n_jobs;

}   
#endif
