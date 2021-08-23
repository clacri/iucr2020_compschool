/* implementation of singly linked list routines.  Based very closely
 * on Kyle Loudon's  implementation as given in "Mastering Algorithms
 * in C".
 */

#include <stdlib.h>
#include <stdio.h>

#include "sll.h"


/* now here are the function implementations */
/* sll_init() -- initialize the linked list */
void sll_init( SLinkedList *list, void (*destroy)(void *data) )
{
    list->size = 0;
    list->destroy = destroy;
    list->match = NULL;
    list->head = NULL;
    list->tail = NULL;

    return;
}

/* sll_destroy() -- does deallocation and garbage cleanup  */

void sll_destroy( SLinkedList *list )
{
    void *data;

    while( sll_list_size(list) > 0 ) {
        if( (sll_remove_next(list, NULL, (void **)&data) == 0) && (list->destroy != NULL) ) {
            list->destroy(data);
        }
    }

/* clear the structure */
    list->size = 0;
    list->destroy = NULL;
    list->head = NULL;
    list->tail = NULL;
    list->match = NULL;
    
    return;
}

/* sll_insert_next() -- insert an element into the list after the specified
 * element.
 */

int sll_insert_next( SLinkedList *list, SLinkedListElem *elem, void *data )
{
    SLinkedListElem *new_elem;

    new_elem = malloc( sizeof( *new_elem ) );
    if( !new_elem )
        return -1;

    new_elem->data = data;

    if( NULL == elem ) {

/* handle insertion at the head */
       if( 0 == sll_list_size(list) ) 
           list->tail = new_elem;

       new_elem->next = list->head;
       list->head = new_elem;
       
    }
    else { /* handle insertion elsewhere in the list */
       if( NULL == elem->next )
           list->tail = new_elem;

       new_elem->next = elem->next;
       elem->next = new_elem;
    }

    list->size++;

    return 0;
}
        
/* sll_remove_next() -- removes the element in the list which is after 
 * the specified element.
 */

int sll_remove_next( SLinkedList *list, SLinkedListElem *elem, void **data )
{
    SLinkedListElem *old_elem;

    if( 0 == sll_list_size(list) )
        return -1;

    if( NULL == elem ) {  /* handle removal from the head of the list */
       *data = list->head->data;
       old_elem = list->head;
       list->head = list->head->next;
       if( 0 == sll_list_size(list) )
            list->tail = NULL;
    }
    else {  /* handle removal from elsewhere in the list */
       if( NULL == elem->next )
           return -1;

       *data = elem->next->data;
       old_elem = elem->next;
       elem->next = elem->next->next;
       if( NULL == elem->next )
           list->tail = elem;
    }

    free( old_elem );
    list->size--;

    return 0;
}
