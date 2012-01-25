#include "graphtypes.h"
#include <stdio.h>
#include <math.h>


#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif


AddEdge (Graph g, int n, int m, int label)
{
  Edge edge1 = (Edge) NULL;
  Edge edge2 = (Edge) NULL;
  long bytes = 2 * sizeof(struct edge_ent);


#ifdef MEX_DEBUG
  mexPrintf( "Debug: Entering AddEdge().\n"  );
  mexPrintf( "Debug:   n     = %d.\n", n     );
  mexPrintf( "Debug:   m     = %d.\n", m     );
  mexPrintf( "Debug:   label = %d.\n", label );
#endif

#ifdef MATLAB_MEX_FILE
	edge1 = (Edge) mxMalloc( bytes );
#else
	edge1 = (Edge) malloc  ( bytes );
#endif

	edge2 = edge1 + 1;

	edge1->label     = label;
	edge1->endpoint  = m;
	edge1->otheredge = edge2;
	edge1->prevedge  = NULL;
	edge1->nextedge  = g[n].adj_list;

	if (edge1->nextedge != NULL)
  {
		edge1->nextedge->prevedge = edge1;
  }

	g[n].adj_list = edge1;
	g[n].degree++;

	edge2->label     = label;
	edge2->endpoint  = n;
	edge2->otheredge = edge1;
	edge2->prevedge  = NULL;
	edge2->nextedge  = g[m].adj_list;

	if (edge2->nextedge != NULL)
  {
		edge2->nextedge->prevedge = edge2;
  }

	g[m].adj_list = edge2;
	g[m].degree++;

#ifdef MEX_DEBUG
  mexPrintf("Debug: Leaving  AddEdge().\n");
#endif
}


Edge FindEdge(graph,i,j)
Graph graph;
int i,j;

{	Edge edge;

	edge = graph[i].adj_list;
	while (edge!=NULL && edge->endpoint!=j)
		edge = edge->nextedge;
	if (edge==NULL) return(NULL);
	else return(edge);
}

int RemoveEdge(graph,edge)
Graph graph;
Edge edge;

{	Edge other;
	int i,j;

	if (edge==NULL) return(0);
	other = edge->otheredge;
	i = other->endpoint;
	j = edge->endpoint;
	graph[i].degree--; graph[j].degree--;
	if (edge->prevedge == NULL) {
		graph[i].adj_list = edge->nextedge;
		if (edge->nextedge != NULL)
			edge->nextedge->prevedge = NULL;
		}
	else if (edge->nextedge == NULL)
        	(edge->prevedge)->nextedge = NULL;
	else {
		(edge->nextedge)->prevedge = edge->prevedge;
		(edge->prevedge)->nextedge = edge->nextedge;
		}
	if (other->prevedge == NULL) {
		graph[j].adj_list = other->nextedge;
		if (other->nextedge != NULL)
			other->nextedge->prevedge = NULL;
		}
	else if (other->nextedge == NULL)
		(other->prevedge)->nextedge = NULL;
	else {
		(other->nextedge)->prevedge = other->prevedge;
		(other->prevedge)->nextedge = other->nextedge;
		}
	/* changed, but appears to be unused by this program anyway */
#ifdef MATLAB_MEX_FILE
	if (edge < other)
	  mxFree(edge);
 	else
	  mxFree(other);
#else
	free((edge < other) ? edge : other);
#endif
	return(1);
}

int NumEdges(g)
Graph g;
{	int i,size,edges;

	edges = 0;
	size = Degree(g,0);
	for (i=1; i<=size; i++)
		edges += Degree(g,i);
	edges /= 2;
	return(edges);
}


Graph NewGraph (int size)
{
  Graph tmp   = (Graph) NULL;
	int   i     = 0;
  long  bytes = (size + 1) * sizeof(struct node_entry);


#ifdef MEX_DEBUG
  mexPrintf("Debug: Entering NewGraph().\n");
#endif


#ifdef MATLAB_MEX_FILE
	tmp = (Graph) mxMalloc( bytes );
#else
	tmp = (Graph) malloc  ( bytes );
#endif

	for (i = 1; i <= size; i++)
  {
		Degree   (tmp, i) = 0;
		FirstEdge(tmp, i) = NULL;
		NLabel   (tmp, i) = i;
  }

	Degree(tmp, 0) = size;


#ifdef MEX_DEBUG
  mexPrintf("Debug: Leaving  NewGraph().\n");
#endif

	return tmp;
}


EuclidGraph NewEuclid(size)
int size;
{
	EuclidGraph xy;

	xy = (EuclidGraph) malloc((size+1)*2*sizeof(int));
	xy[0][0] = size;
	return(xy);
}

MatrixGraph NewMatrix(size)
int size;
{
	MatrixGraph graph;
	int i;

	graph = (MatrixGraph) malloc((size*(size+1)+1)*sizeof(int));
	graph[0] = size;

	for (i=1; i<=size; i++)		/* zero the diagonal */
		graph[i*(size+1)] = 0;

	return(graph);
}

Graph CopyGraph(g)
Graph g;
{	int i,j,size;
	Edge edge;
	Graph cp;

	size = Degree(g,0);
	cp = NewGraph(size);
	for (i=1; i<=size; i++) {
		Xcoord(cp,i) = Xcoord(g,i);
		Ycoord(cp,i) = Ycoord(g,i);
		edge = FirstEdge(g,i);
		for (j=1; j<=Degree(g,i); j++) {
			if (i < EndPoint(edge))
				AddEdge(cp,i,EndPoint(edge),ELabel(edge));
			edge = NextEdge(edge);
			}
		}
	return (cp);
}

/* Euclidean distance routines */

int eucdist (graph,i,j) /* Find the distance between two points */
			/* 10K x 10K unit square */
EuclidGraph graph;
int i,j;
{	int dv,dh;
	register int k, l;

	dv = graph[i][0]-graph[j][0];
	dh = graph[i][1]-graph[j][1];
	k = dv*dv + dh*dh;
	if (k==0) return(0);
	if (dv<0) dv = -dv;
	if (dh<0) dh = -dh;
	l = dv + dh;
	l = (l + k/l)>>1;
	l = (l + k/l)>>1;
	l = (l + k/l)>>1;
	l = (l + k/l)>>1;
	return ((l*l<k) ? ++l : l);
}


int eucdist2 (graph,i,j) /* Find the distance between two points */
			/* 1M x 1M unit square */
EuclidGraph graph;
int i,j;
{	double dv,dh,d;
	int l;

	dv = (double) graph[i][0]-graph[j][0];
	dh = (double) graph[i][1]-graph[j][1];
	d  = sqrt(dv*dv + dh*dh);
	l  = (int) d;
	return((d-l > .000000001) ? l+1 : l);
}


int eucdistsq(graph,i,j) /* Find the square of the dist between two points */
EuclidGraph graph;
int i,j;
{
	int dv,dh;

	dv = graph[i][0]-graph[j][0];
	dh = graph[i][1]-graph[j][1];
	return(dv*dv+dh*dh);
}

