/* Implementation of a singly linked lists.  Based closely on Kyle Loudon's
 * singly linked lists as specified and implemented in his book,
 * "Mastering Algorithms in C"
 */

#ifndef SLL_H
#define SLL_H

#include <stdlib.h>

/* Here are the data types */
/* define the data structures for the list element and the list as whole
 * itself.
 */

typedef struct sll_elem SLinkedListElem;
typedef struct sll SLinkedList; 

struct sll_elem {
       struct sll_elem *next;
       void *data;
       };

struct sll {
       void (*destroy)( void *data );
       int (*match)(const void *key1, const void *key2 );
       int size;
       SLinkedListElem  *head;
       SLinkedListElem  *tail;
       };


/* user callable macros and functions */
void sll_init( SLinkedList *list, void (*destroy)(void *data) );
void sll_destroy( SLinkedList *list );
int sll_insert_next( SLinkedList *list, SLinkedListElem *elem, void *data );
int sll_remove_next( SLinkedList *list, SLinkedListElem *elem, void **data );

/* here are some utility macros */
#define sll_list_size(list) ((list)->size)
#define sll_list_head(list) ((list)->head)
#define sll_list_tail(list) ((list)->tail)
#define sll_is_head(list, element) ((element) == (list)->head ? 1 : 0 )
#define sll_is_tail(element) ((element)->next == NULL ? 1 : 0 )
#define sll_list_data(element) ((element)->data)
#define sll_list_next(element) ((element)->next)
#endif
