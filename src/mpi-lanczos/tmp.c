/* This is nowhere near done - please don't look at it! :)
*/

#include <stdio.h>
#include "mpi.h"

#include <time.h>


#define MAX_NODES 64

#define WORKTAG 1
#define DIETAG 2


/* These are for the status field. */
#define ST_MASTER  0
#define ST_WORKING 1
#define ST_ONCALL  2
#define ST_DEAD    3

typedef struct {
  int   nodeID;
  long  c0;
  long  c1;
  int   status;
} node_t;

/* Globals for the master node. */
node_t nodeList[MAX_NODES];
int    numNodes, numSlaves, numOnCall;


/*********************************************************************/
int initializeNode(node_t *N)
/*********************************************************************/
/* Read the columns N->c0,..., (N->c1-1) of the matrix from file and */
/* send them to node N->nodeID.                                      */
/*********************************************************************/
{
#ifdef _NO
    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
#endif
}



/********************************************************************/
static void master(void)
/********************************************************************/
/* One of the MPI processes will be distinguished as a master node. */
/* That master node will call this function only, so this is the    */
/* function that controls all the slave nodes.                      */
/********************************************************************/
{ node_t           nodeList[MAX_NODES];
  int              numNodes, j, numSlaves, numOnCall;
  long             colsPerNode, c0, c1;
  nfs_sparse_mat_t M;
  double           blstart, blstop;
  long            *deps;

  /* Load the matrix. */

  /* This is not yet done, as it could be a bit tricky.
     It might not fit into RAM on a single node, which means 
     we should really preprocess it first on a machine with
     plenty of RAM, save the preprocessed version to disk,
     and use that one here. Then, instead of loading the
     whole matrix at once, the master node should just load
     columns as they are needed to send to the workers.
       All this function really needs to know is the matrix
     stats - dimensions and such. So we should here fill those
     fields in.  Then make the initializeNode() function above
     read the needed columns from disk and send them out.
  */
  M.numCols = M.numRows = 0; /* Change this to the real deal. */



  /* Find out how many nodes there are. */
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
  nodeList[0].status = ST_MASTER;

  /****************************************************************/
  /* If we have more than a handful of nodes, say more than 8,    */
  /* we should leave at least one `on call' in case another one   */
  /* dies. Why? Consider the alternative: if a node dies and we   */
  /* don't have a backup available, then we either have to        */
  /* give one node twice the work, or re-distribute and re-assign */
  /* the matrix columns. This would:                              */
  /*  (1) Hugely complicate the code.                             */
  /*  (2) Waste bandwidth.                                        */
  /*  (3) Waste time.                                             */
  /* Of course, we could assume that no node ever fails, but that */
  /* is a very bad idea. The is the best compromise between       */
  /* robustness and efficiency that I can think of.               */
  /* If you don't like it, you can change this single line of     */
  /* code, but beware that the resulting code will not be fault   */
  /* tolerant at all!                                             */
  /****************************************************************/
  numOnCall = ((numNodes >= 8) ? (MIN(1, 0.05*numNodes)) : 0);

  for (i=1; i<numNodes; i++)
    nodeList[i].status=ST_WORKING;
  for (i=0; i<numOnCall; i--)
    nodeList[numNodes-1-i].status = ST_ONCALL;

  /* How many columns does each node get? */
  for (i=1; i<numNodes; i++)
    if (nodeList[i].status == ST_WORKING)
      numSlaves++;
  colsPerNode = M->numCols/numSlaves;

  /* Initialize the slave nodes with their needed columns. */
  c0 = 0; c1 = c0 + colsPerNode;
  for (i=0,j=0; i<numSlaves; i++) {
    do {
      j++;
    } while (nodeList[j].status != ST_WORKING);
    nodeList[j].c0 = c0;
    nodeList[j].c1 = c1;
    initializeNode(&nodeList[j]);
    c0 = c1;
    c1 = MIN(c0+colsPerNode, M->numCols);
  }

  /* We should now be able to just do: */
  srand(seed);
  seedBlockLanczos(seed);
  startTime = sTime();
  msgLog("", "GGNFS-%s : mpi-matsolve", GGNFS_VERSION);
  msgLog("", "Running on %d (of %d) nodes.", numSlaves, numNodes);
  

  blstart = sTime();
  if (!(deps = (long *)malloc(M.numCols*sizeof(long)))) {
    fprintf(stderr, "Severe memory allocation error!\n");
    res = -1;
  } else {
    res = blockLanczos32(deps,  mpi-multB32, mpi-multB_T32, NULL, M.numCols);
    blstop = sTime();
    printf("Returned %d. Block Lanczos took %1.2lf seconds.\n", res, blstop-blstart);
    msgLog("", "BLanczosTime: %1.1lf", blstop-blstart);
  }

  /* Shutdown all of the nodes except this one. */
  for (i=1; i<numNodes; i++) {
    MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
  }

  /* Now retrieve the actual dependencies (i.e., if there was a matrix
     pruning step, then we need to get back to the original matrix).
     Then store them to file.
  */


}

/***********************************/
static void slave(void)
/***********************************/
{ unit_of_work_t work;
  unit_result_t results;
  MPI_Status status;
  nfs_sparse_mat_t M;

  /* The first message should be an initialization,
     containing the columns that this node is responsible for.
     Upon receiving this message, we should call
        setupPartialMat()
     To setup our local copy of the matrix.
  */
  MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
           MPI_COMM_WORLD, &status);


  while (1) {
    /* Receive a message from the master. This should
       contain a vector and instructions to multiply
       it with B or B^T.
    */
    MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    /* Check the tag of the received message. */
    if (status.MPI_TAG == DIETAG) {
      return;
    }

    /* If we did not get a shutdown message, there is
       work to do. Call either MultB_partial32() or
       MultB_T_partial32() as appropriate. 
    */
    result = do_work(work);

    /* And send the result back */
    MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}


/*************************************************/
int mpi-multB32(long *Product, long *x, void *P)
/*************************************************/
/* This function should send out the vector 'x'  */
/* to the slave nodes and have them compute the  */
/* partial product of B|_{their columns} by 'x'. */
/* Then, gather up the results and combine them  */
/* to form the actual product.                   */
/*************************************************/
/* There must be a layer of fault tolerance here.*/
/* If a node dies during the computation, try to */
/* reassign it's task to an onCall node.         */
/*************************************************/
{
/*
  Use the functions MPI_Send() and MPI_recv() to
  send and receive work units to the slave nodes
  until the matrix multiplication is complete.
  If a node dies, re-assign it's columns to an
  onCall node if one is available. If not,
  simply fail.
*/
}

/*************************************************/
int mpi-multB_T32(long *Product, long *x, void *P)
/*************************************************/
/* Same as above, but with B^T instead of B.     */
/*************************************************/
{
}




/*************************************************/
int main(int argc, char *args[])
/*************************************************/
/* Entry point for master and slave nodes alike. */
/*************************************************/
{ int rank, size;
  time_t now, start;

 
  MPI_Init(&argc, &args);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("Hello world! I am %d of %d.\n", rank, size);

  if (rank == 0) {
    /* I am the controlling process. */
    master();    
  } else {
    slave();
  }

  printf("This is %d of %d signing off...\n", rank, size);
  MPI_Finalize();
  return 0;
}

#ifdef _NO

  /* This is the sample code I'm using as a template to
     build the MPI code. I'm keeping it here for handy reference.
  */

  /* Send a work unit to each node. */
  for (j=1; j<numNodes; j++) {
    work = get_next_work_item();
    /* Send it to each rank */
    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
  }

  /* Loop over getting new work requests until there is no more work
     to be done */
  work = get_next_work_item();
  while (work != NULL) {
    /* Receive results from a slave */
    MPI_Recv(&result,           /* message buffer */
             1,                 /* one data item */
             MPI_DOUBLE,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

    /* Send the slave a new work unit */
    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    /* Get the next unit of work to be done */
    work = get_next_work_item();
  }

  /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */
  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

#endif
