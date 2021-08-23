# Makefile for COSET
# To install type the following commands:
# make
# make install
# make clean

#Linux make variables:


#if fork() and execl() system calls are not available on your
#system, then compile with -DUSE_SYSTEM_FUNCTION added to the 
#CFLAGS variable.  To compile on systems where the C99 functions,
#isblank(), snprintf(), and vsnprint() are not available, add
#the following to the CFLAGS variable:

#-DNEED_C89_COMPATIBILITY

#if you want to use the extended Flack B algorithm which lists
#centrically related twin laws, add the following to the CFLAGS
#variable:

#-DUSE_EXTENDED_B_ALGORITHM

#if you want to compile the code as a Python extension module add this
#to the CFLAGS variable:


#  Development flags (for gcc only)
CFLAGS = -g -std=c99 -pedantic  -Wextra -Wall \
         -Wshadow -Wbad-function-cast \
         -Wstrict-prototypes -DUNIXY_SUFFIX -DUSE_EXTENDED_B_ALGORITHM

PYTHON_INCLUDE = /usr/include/python2.7
PYFLAGS = -fPIC -DPYTHON_EXTENSION_MODULE 
PYMODULE_NAME = FlackCoset.so 
PYTHON_CODE_DIR = /usr/local/python_code

PYTHON_CFLAGS = -c -I$(PYTHON_INCLUDE) $(PYFLAGS)  $(CFLAGS) 

# install the executable here.  We choose /usr/local/bin on
# our local setup, you may want this different.
INSTALL_BIN_DIR =/usr/local/bin

#Production flags for Linux using gcc.
#CFLAGS = -std=c99 -pedantic  -Wextra -Wall \
#         -Wshadow -Wbad-function-cast \
#         -Wstrict-prototypes -DUNIXY_SUFFIX

#Other platforms

#SGI CFLAGS.
#CFLAGS = -ansi -mips2 -O3

#DEC OSF/1 FLAGS
#CFLAGS = std1 -O3

#version of the program
VERSION := $(shell grep '\#define VERSION ' version.h | awk '{print $$3}' | tr -d \")

# get which numerical library for the eigen code is to be used.  Should be either 'lapacke' or 'gsl'
EIGEN := $(shell grep '\#define EIGEN_CODE' version.h | awk '{print $$3}' | tr -d \")


#uncomment to give executable file different names depending on
#eigen code used.
USE_SPECIAL_EXE_NAME = true

ifeq ($(EIGEN),lapacke)
LIBS   = -llapack -llapacke -lm
EIGEN_SRC = eigen_lapacke.c
else
LIBS   = -lgslcblas -lgsl -lm
EIGEN_SRC = eigen_gsl.c
endif

ifdef USE_SPECIAL_EXE_NAME
EXE_NAME = coset-$(EIGEN)
else
EXE_NAME = coset
endif

#General definitions.
 
EXE    = $(EXE_NAME)
SRCS =   coset.c \
         c89_util.c \
         dupstr.c \
         dynamic_sll.c \
         $(EIGEN_SRC)  \
         float_util.c \
         input.c \
         main.c \
         matrix.c \
         queue.c \
         shelx.c \
         shelx_exec.c \
         sll.c \
         symm_mat.c \
         task.c \
         usage.c

OBJS = $(SRCS:.c=.o)

$(EXE): $(OBJS)
	$(CC)  -o $@ $(OBJS) $(LIBS)

$(PYMODULE_NAME): $(SRCS)
	rm -f *.o
	$(CC) $(PYTHON_CFLAGS) $?
	$(CC) -shared -fPIC $(OBJS) -o $@ $(LIBS)

.c.o:
	$(CC) -c $(CFLAGS) $<


install:
	mv $(EXE) $(INSTALL_BIN_DIR)

install_py_module:
	mv $(PYMODULE_NAME) $(PYTHON_CODE_DIR)

clean:
	rm -f *.o $(EXE) $(PYMODULE_NAME)

archive:
	cd ../; tar -zcvf coset-$(VERSION).tar.gz --exclude=.svn --exclude='*.o'  coset-$(VERSION)

zip:
	cd ../; zip -r coset-$(VERSION).zip coset-$(VERSION) --exclude '*.svn*' \*.o
