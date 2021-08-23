/* contains public interface for running SHELXL jobs to
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

#define _ISOC99_SOURCE
#ifndef SHELX_EXEC_H
#define SHELX_EXEC_H

/* Change for coset-1.2.0 -- we now trap the output which SHELXL sends to
 * stdout.  We can have the old spew to the screen behaviour by adding
 * -DDONT_TRAP_SHELX_STDOUT to the CFLAGS in the Makefile.  Otherwise,
 * we go with the new behaviour by default. Set the macros below to configure.
 */
#ifndef DONT_TRAP_SHELX_STDOUT
#define TRAP_SHELX_STDOUT
#define COMMAND_BUF_FORMAT "%s %s > %s%s"

#ifdef UNIXY_SUFFIX  /* for people who enjoy C and UNIX programming conventions */
#define TRAP_SUFFIX ".stdout"
#else                /* for  the unwashed masses of other platforms */
#define TRAP_SUFFIX ".screen"
#endif

#else
#undef  TRAP_SHELX_STDOUT
#define COMMAND_BUF_FORMAT "%s %s"
#endif 


/* interface protypes */
char *real_hklf_filename(char *ins_filename);
int create_hklf_file_links(const char *real_hklf, SLinkedList *link_names);
SLinkedList *make_job_name_list(SLinkedList *ins_file_names);
SLinkedList *create_hkl_filenames(SLinkedList *job_names);
SLinkedList *setup_shelx_jobs(SLinkedList *ins_file_names, char *orig_ins_filename );
int spawn_shelx_jobs(SLinkedList *jobs, char *shelx_exe_path);
#endif

