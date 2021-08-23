/* this is an implementation of a queue data structure.  Based very closely
 * on Kyle Loudon's "Mastering Algorithms in C", pp. 105-110.  The 
 * implementation is based on a singly linked list.
 */

#include <stdlib.h>

#include "sll.h"
#include "queue.h"

int queue_enqueue( Queue *queue, void *data )
{
    return sll_insert_next(queue, sll_list_tail(queue), data );
}

int queue_dequeue( Queue *queue, void **data )
{
    return sll_remove_next( queue, NULL, data );
}
