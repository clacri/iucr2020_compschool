#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "dupstr.h"

#define LEN  64


int main( void )
{
    char tmp[LEN];
    char *filler = "the dog jumps high.";
    char *p;

    strcpy( tmp, filler );
    p = dupstr( tmp );

    printf( "Duplicate string: %s\n", p );
    free( p );
    exit( EXIT_SUCCESS );
}

