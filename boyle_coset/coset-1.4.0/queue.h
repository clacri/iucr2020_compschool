/* this is an implementation of a queue data structure.  Based very closely
 * on Kyle Loudon's "Mastering Algorithms in C", pp. 105-110.  The 
 * implementation is based on a singly linked list.
 */
#ifndef QUEUE_H
#define QUEUE_H

#include <stdlib.h>

#include "sll.h"  /* my own singly linked list stuff */

typedef SLinkedList Queue;

int queue_enqueue( Queue *queue, void *data );
int queue_dequeue( Queue *queue, void **data );

#define queue_init  sll_init
#define queue_destroy sll_destroy
#define queue_peek(queue) ((queue)->head == NULL ? NULL : (queue)->head->data)
#define queue_size sll_list_size
#endif
