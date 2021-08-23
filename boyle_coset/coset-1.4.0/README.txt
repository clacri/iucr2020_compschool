COSET is a program which derives potential merohedral and pseudomerohedral
twin laws.  In addition, the program can produce SHELX .ins file
containing TWIN instructions for the twinning possibilities and can test
them by executing SHELXL refinements automatically (assuming that a locally
installed SHELXL binary is available).  This allows for fast and easy
screening of potential twin laws.

The program performs left coset decompositions of a metrically available
supergroup point group symmetry by the point symmetry of a crystal in
which (pseudo)merohedral twinning is present.  The algorithms used in
this program are described in Howard Flack's paper: Acta Cryst. (1987),
A43, 564-568.

The input is simple and flexible, allowing the user to determine how
extensive of an analysis and testing of potential twin laws is to be
performed.  The program is released under the GNU Public License (GPL).

INSTALLATION:
First, edit the EIGEN_CODE macro in 'version.h' file to choose
which eigen code to use.  Allowed values are "lapacke" and "gsl".
 
The program is written in ANSI C99 (although it uses only a few
C99 features) so getting to work with a C89 compliant compiler, 
just add the following to the CFLAGS variable in the Makefile:

-DNEED_C89_COMPATIBILITY

When so desired, the program makes use of Unix system calls fork()
and execl() to spawn off SHELXL refinements to test trial twin laws.
If your system does not have fork(), symlink(), wait(), and execl(),
then add the following to the CFLAGS variable in the Makefile:

-DUSE_SYSTEM_FUNCTION

This will use the ANSI C system() function to execute the SHELXL binary
executable.

To compile this on a MS Windows machine you will need to add both
preprocessor directives, -DNEED_C89_COMPATIBILITY and -DUSE_SYSTEM_FUNCTION
to the CFLAGS variable in the Makefile.

If you want to use an extended version of the B algorithm which prints out
the centrically related twin laws for acentric structures, add 
-DUSE_EXTENDED_B_ALGORITHM to the CFLAGS variable in the Makefile.  See 
comments in the coset.c source file.

While editing the Makefile, consider where you want the COSET executable
to be stored. By default, it gets installed to /use/local/bin/.  Edit the
Makefile to suit your own preferences.

To compile the program simply type:
make
followed by:
make install

Python Extension Module of COSET 
As of version 1.1.0, the program can also be compiled as a Python extension module:
make FlackCoset.so
make install_py_module

For the Python extension module to work, the PYTHONPATH variable should be set.  The
Makefile assumes this is set to /usr/local/python_code.  Once the module is installed,
invoke the python interpreter, and type:
>>> import FlackCoset

The function is called "decomp" and is called like so:
>>> FlackCoset.decomp( '<input_filename>' ) 
where <input_filename> is the name of the COSET input file described below.

PROGRAM EXECUTION:
The usage is as follows:

On the command line, type:

coset <input_filename>

Typing the program name without an input file displays a help screen
on the terminal.

<input_filename> is the name of a plaintext file which contains a number
of directives and parameters used to govern the execution of the program.
The program can process a number of coset analyses in a given execution,
which are designated as "tasks".  The directive and parameters formats
are:

TITLE <character string data>  [mandatory]
ALGORITHM   <single character> [mandatory]
SUPERGROUP  <character string data> [mandatory]
SUBGROUP  <character string> < integer> [mandatory]
RMAT  r11 r12 r13 r21 r22 r23 r31 f32 r33 [mandatory]
TRANS t11 t12 t13 t21 t22 t23 t31 t32 t33 [optional]
INSFILE  <character string data> [optional but needs TRANS]
OUTFILE  <character string data> [optional]
NEWINS   <character string data> [optional but needs INSFILE and TRANS]
EXEC     <character string data> [optional but needs TRANS and NEWINS]
END     

The '#' character at the beginning of a line designates a comment and
is ignored by the program.

Explanation of input directives and parameters:
* TITLE  must be the first directive for a given task.  The parameter
  for the this directive is a short description of less than 75 characters.

* ALGORITHM allowed parameters are 'A' or 'B' (without quotes).

* SUPERGROUP parameter is one of the following character strings:
  -1
  2/m
  mmm
  4/mmm
  -3m
  6/mmm
  m-3m

*SUBGROUP takes two parameters.  The first being a character string designation
 for the crystal's point group, e.g. -3 or mm2.  The second parameter is
 an integer which is equal to the number of symmetry operators for the
 crystal's point group.  Including the identity operator.

*RMAT takes 9 numeric parameters, which are the matrix elements for the
 symmetry operators for the point group.  Each symmetry operator takes a
 separate RMAT statment.  The number of RMAT statements should equal the 
 numeric value given in the SUBGROUP statement.  When constructing the
 RMAT statements, simply translate the International Table's equivalent
 positions (in direct space) to the matrix representation.  For example,
 if the crystal's spacegroup is I 4(1)/a, then drop the lattice centering
 symbol and convert the translation symmety elements to the non-translational
 equivalents, so with I 4(1)/a one would use the equivalent positions for
 P 4/m. The identity operator must be included in this list.

*TRANS takes 9 numeric elements which transform the crystal's unit cell 
 parameters to the metrically available supergroup cell.  These elements
 are normally obtained from a cell reduction program.  If TRANS is omitted
 the identity matrix is used.

*INSFILE takes a single character string which is filename of the SHELX
 .ins file for the structure.  This file is not altered by the program but
 provides basis for new .ins file(s) used to do twinned refinements.

*OUTFILE takes a single character string which is filename for the general
 output from COSET.  If this directive is not specified, the program writes
 these results to stdout ("standard output", i.e. the terminal).

*NEWINS takes a single character string which is the basename for the
 new set of .ins file which incorporate the SHELX BASF and TWIN instructions
 for twinned refinement.

*EXEC takes a single character string which is the full pathname of the
 local system's SHELX(T)L executable.  With UNIX systems, symbol links
 are created to the actual .hkl file.  In systems which use the system()
 call to spawn SHELXL jobs, the .hkl file is copied to each
 <new_basename>.hkl file.  This wastes disk space, but is portable.

*END takes no paramters and should be the last line of the file.

Example inputs:
Below is the input file for COSET for Regine's Herbst-Irmer's
pseudo-merohedral aniline case.  The results are written to 
stdout (i.e. "the terminal").

Here is a minimal case:
TITLE anilin test case
ALGORITHM A
SUPERGROUP mmm
SUBGROUP 2/m  4
RMAT 1 0 0 0 1 0 0 0 1
RMAT -1 0 0 0 1 0 0 0 -1
RMAT -1 0 0 0 -1 0 0 0 -1
RMAT 1 0 0 0 -1 0 0 0 1
TRANS 0 0 1 2 0 1 0 1 0
END

Below is the same input file, but with more of the directives used.
In this case, the user specifies that the results of the coset
decomposition will be written to a file called 'anilin.cosets', the
existing file 'anilin.ins' will be used to generate new .ins files which
have the basename 'twin', and finally, the new .ins files will be used
as input for the SHELXL executable binary located in /usr/local/bin/xl.
This, of course, assumes that the .hkl file is in the same directory as
the .ins file.

TITLE anilin test case
ALGORITHM A
SUPERGROUP mmm
SUBGROUP 2/m  4
RMAT 1 0 0 0 1 0 0 0 1
RMAT -1 0 0 0 1 0 0 0 -1
RMAT -1 0 0 0 -1 0 0 0 -1
RMAT 1 0 0 0 -1 0 0 0 1
TRANS 0 0 1 2 0 1 0 1 0
OUTFILE anilin.cosets
INSFILE anilin.ins
NEWINS twin
EXEC /usr/local/bin/xl
END


I hope you find this software useful.  Please direct questions,
comments, suggestions, and bug reports to the author, Paul D. Boyle,
at pboyle@uwo.ca
