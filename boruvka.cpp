#include <cstdlib>
#include <assert.h>
#include <bits/stdc++.h>

using namespace std;

typedef struct Edge {
    int i;
    int j;
    int w;
} Edge;

typedef struct Node {
    int parent;
    int rank;
} Node;

void init_edge(Edge* edge, int i, int j, int w) {
    edge->i = min(i, j);
    edge->j = max(i, j);
	edge->w = w;
}

int cmp_edges(const void* left_ptr, const void* right_ptr) {
    const Edge* a = (const Edge*)left_ptr;
    const Edge* b = (const Edge*)right_ptr;

    // If the weights are unequal, prefer smaller weight
    if (a->w != b->w)
        return a->w - b->w;

    // If the weights are equal, prefer smaller first vertex
    if (a->i != b->i)
        return a->i - b->i;

    // Prefer smaller second vertex
    return a->j - b->j;
}

Node* create_nodes(int N) {
	Node* nodes = (Node*)calloc(N, sizeof(Node));
	for (int node = 0; node < N; ++node) {
		nodes[node].parent = node;
		nodes[node].rank = 1;
	}
	return nodes;
}

void print_tree(Edge* tree, int numEdges) {
	int sum = 0;
	for (int edge = 0; edge < numEdges; ++edge) {
        cout << tree[edge].i << " " << tree[edge].j << endl;
		sum += tree[edge].w;
	}    
    cout << "SUM : " << sum << endl;
}

// Find root of the tree
int find(Node* nodes, int node) {
	if (nodes[node].parent != node) {
		nodes[node].parent = find(nodes, nodes[node].parent);
	}
	return nodes[node].parent;
}

// Merge trees
void merge_trees(Node* nodes, int root1, int root2) {

    // Merge smaller rank tree into larger rank tree.
	if (nodes[root1].rank > nodes[root2].rank) {
		merge_trees(nodes, root2, root1);
		return ;
	}

    // Make root1 tree as child of root2 tree.
	nodes[root1].parent = root2;
	if (nodes[root1].rank == nodes[root2].rank) {
		nodes[root2].rank += 1;
	}
}

// Serial Boruvka's algorithm
int union_find(Edge* edges, int N, int M, Edge *tree) {
    int numEdges = 0;

    // Initially all nodes in a tree of its own.
	Node* nodes = create_nodes(N);

    // An array to store index of the cheapest edge of
    // node set. The stored index for indexing array 'edges[]'
    int *cheapest = new int[N];

    bool canMerge = true;
 
    // Keep combining components (or sets) until we can
    while (canMerge)
    {
        canMerge = false;
        for(int i=0; i<N; i++) cheapest[i] = -1;

        // Traverse through all edges and update
        // cheapest of every component
        for (int edge=0; edge<M; edge++)
        {
            // Find components (or node sets) of two corners
            // of current edge
            int root1 = find(nodes, edges[edge].i);
            int root2 = find(nodes, edges[edge].j);
 
            // If two corners of current edge belong to
            // same set, ignore current edge
            if (root1 == root2)
                continue;
 
            // Else check if current edge is closer to previous
            // cheapest edges of set1 and set2
            else
            {
               if (cheapest[root1] == -1 ||
                   edges[cheapest[root1]].w > edges[edge].w)
                 cheapest[root1] = edge;
 
               if (cheapest[root2] == -1 ||
                   edges[cheapest[root2]].w > edges[edge].w)
                 cheapest[root2] = edge;
            }
        }
 
        // Consider the above picked cheapest edges and add them
        // to MST
        for (int node=0; node<N; node++)
        {
            // Check if cheapest for current set exists
            if (cheapest[node] != -1)
            {
                int root1 = find(nodes, edges[cheapest[node]].i);
                int root2 = find(nodes, edges[cheapest[node]].j);
 
                if (root1 == root2)
                    continue;

                canMerge = true;

                // Select this edge
                // copy the forest edges into the tree variable
                memcpy(tree + numEdges++, edges + cheapest[node], sizeof(Edge));
 
                // Do a union of root1 and root2 and decrease number
                // of trees
                merge_trees(nodes, root1, root2);
            }
        }
    }

	free(nodes);

	return numEdges;
}

// Sets of edges received are already sorted. Just merge them.
void merge_lists(Edge* li, int ni, Edge* lj, int nj, Edge* output) {
	int i = 0;
	int j = 0;
	while (i < ni) {
		memcpy(output+i, li+i, sizeof(Edge));
		i += 1;
    }
    while(j < nj) {
		memcpy(output+i+j, lj+j, sizeof(Edge));
		j += 1;
    }
}

// Wrapper for merging sorted lists of edges.
void merge_into(Edge* oldBuf, int n1, Edge* newBuf, int n2) {
	Edge* output = (Edge* )calloc(n1+n2, sizeof(Edge));
	merge_lists(oldBuf, n1, newBuf, n2, output);
	memcpy(newBuf, output, (n1+n2)*sizeof(Edge));
	free(output);
}

int receive_edges_from(int stepSize, int targetRank, Edge* edges) {
    int numEdges;
 
    // Receive number of edges.
    MPI_Recv(&numEdges, 1, MPI_INT, targetRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Edge buffer is int buffer times 3.
    int* buffer = (int* )calloc(numEdges * 3, sizeof(int));
    MPI_Recv(buffer, numEdges * 3, MPI_INT, targetRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    for (int i = 0; i < numEdges; ++i) {
        int* edge = &buffer[i*3];
        edges[i].i = edge[0];
        edges[i].j = edge[1];
        edges[i].w = edge[2];
    }
    
    return numEdges;
}

int receive_edges(int procRank, int numProcs, int stepSize, int squareSize, Edge* edges,
    int numEdges, int numEdgesMax) {
    
    Edge* tmpBuf = (Edge* )calloc(numEdgesMax, sizeof(Edge));
    Edge* buf1 = edges;
    Edge* buf2 = tmpBuf;

    Edge* newEdges = (Edge* )calloc(stepSize * squareSize, sizeof(Edge));
    int parity = 0;
    for (int i = 0; i < stepSize; ++i) {
        if (procRank + stepSize + i >= numProcs) break;
        int numNewEdges = receive_edges_from(stepSize, procRank + stepSize + i, newEdges);
        merge_lists(buf1, numEdges, newEdges, numNewEdges, buf2);
        numEdges += numNewEdges;
        
        // Swap buffers buf1 and buf2
        Edge* tmp = buf1;
        buf1 = buf2;
        buf2 = tmp;        

        parity = 1 - parity;
    }

    // If last loop was run odd number of times, 
    if (parity != 0) {
        memcpy(edges, tmpBuf, numEdges*sizeof(Edge));
    }

    free(newEdges);
    free(tmpBuf);
    
    return numEdges;
}

int receive_forest(int procRank, int numProcs, int stepSize, Edge* edges) {
	int target = procRank + stepSize;
	if (target >= numProcs)
		return 0;
	return receive_edges_from(stepSize, target, edges);
}

int receive_new_edges(
	int procRank, int numProcs, int stepSize, int squareSize, int numRows,
	Edge* edges, int numMaxEdges, Edge* forest, int forestSize) {
	
    // Receive forest of i+2^k.
    int newForestSize = receive_forest(procRank, numProcs, stepSize, edges);
	
    // Receive edges with processes in range i+1 to i+2^k
    int numNewEdges = receive_edges(procRank, numProcs, stepSize, squareSize,
        edges, newForestSize, numMaxEdges);
	
    // Merge the two type of sorted edge sets
    merge_into(forest, forestSize, edges, numNewEdges);
	
    // Total edges received are forest edges + other edges
    numNewEdges += forestSize;
	
    return numNewEdges;
}

/* Add the edges which are connecting to target process*/
int add_bipartite_edges(int procRank, int stepSize, int* adj, int numRows, int N, Edge* edges) {
    // At this point target process' rows start from start and it is responsible for
    // numRows*stepSize number of rows.
	int start = ((procRank - procRank%stepSize) - stepSize)*numRows;
    int end = start + numRows*stepSize;

    int numEdges = 0;
    for (int i = 0; i < numRows; ++i) {

        // My vertices' id's
        int procIdx = procRank*numRows + i;
    
        if (procIdx >= N) break;

        for (int j = start; j < end; ++j) {
            if (adj[i*N + j] == 0) continue ;
            init_edge(edges+numEdges, procIdx, j, adj[i*N + j]);
            numEdges += 1;
        }
    }
    return numEdges;
}

/* Add the edges which are incident on vertices allotted to me. */
int add_local_edges(int procRank, int* adj, int numRows, int N, Edge* edges) {
    int numEdges = 0;

    for (int i = 0; i < numRows; ++i) {

        // My vertices' id's
        int procIdx = procRank*numRows + i;
    
        if (procIdx >= N) break;

        for (int j = procRank*numRows; j <= procIdx; ++j) {
            if (adj[i*N + j] == 0) continue ;
            init_edge(edges+numEdges, procIdx, j, adj[i*N + j]);
            numEdges += 1;
        }
    }
    return numEdges;
}

// Send edges to process with rank target.
void send_edges_to(int target, Edge* edges, int numEdges) {

    // Convert edges to integer buffer 
    int* buf = (int* )calloc(numEdges*3, sizeof(int));
    for (int i = 0; i < numEdges; ++i) {
        int* edge = buf + i*3;
        edge[0] = edges[i].i;
        edge[1] = edges[i].j;
        edge[2] = edges[i].w;
    }

    // Send number of edges
    MPI_Send(&numEdges, 1, MPI_INT, target, 0, MPI_COMM_WORLD);

    // Send edges
    MPI_Send(buf, numEdges * 3, MPI_INT, target, 0, MPI_COMM_WORLD);
    
    free(buf);
}

void send_bipartite_forest(int procRank, int stepSize, int squareSize, int numRows,
    int* adj, int N) {
	
    int target = (procRank - procRank%stepSize) - stepSize;

    // Max number of edges possible are squareSize*stepSize*stepSize, because vertex size
    // size doubles each time.
	Edge* edges = (Edge* )calloc(squareSize*stepSize*stepSize, sizeof(Edge));
	
    int numEdges = add_bipartite_edges(procRank, stepSize, adj, numRows, N, edges);
    
    // Sort the edges by smaller weight, then smaller vertex indices.
    // qsort(edges, numEdges, sizeof(Edge), cmp_edges);

    // Take a spanning forest of these bipartite edges and send.
    Edge* forest = (Edge* )calloc((stepSize+1)*numRows-1, sizeof(Edge));
	
    int forestSize = union_find(edges, N, numEdges, forest);
	send_edges_to(target, forest, forestSize);
	
    free(forest);
	free(edges);
}

void send_forest(int procRank, int stepSize, Edge* forest, int forestSize) {
	int target = procRank - stepSize;

    // Sending its forest means, sending all the forest edges.
	send_edges_to(target, forest, forestSize);
}

int create_forest(int procRank, int *adj, int numRows, int N, Edge* forest) {
    // Local edges
    Edge* edges = (Edge* )calloc(numRows*numRows, sizeof(Edge));

    // Add the edges which are incident on vertices allotted to me.
    int numEdges = add_local_edges(procRank, adj, numRows, N, edges);

    // Sort the edges by smaller weight, then smaller vertex indices.
    // qsort(edges, numEdges, sizeof(Edge), cmp_edges);

    // This creates the MSF from the sorted edges as in serial boruvka's algorithm
    int forestSize = union_find(edges, N, numEdges, forest);

    free(edges);
    
    return forestSize;
}

void parallel_boruvka(int procRank, int numProcs, int *adj, int N, int M, Edge* tree) {
	// Number of vertices per process, initially.
	int numRows = ceil((float)N / (float)numProcs);
	int squareSize = numRows * numRows;

	// Get the MSF for vertices allotted to me. This step can be done in parallel.
	Edge* forest = (Edge* )calloc(numRows - 1, sizeof(Edge));
	int forestSize = create_forest(procRank, adj, numRows, N, forest);

    // Each process starts in a receiver state at the beginning
	int receiver = 1;

    // A given processor is a receiver until it sends its first forest. After
    // that, it's not a receiver anymore. It becomes a sender. It will send
    // edges between two other forests for the rest of the algorithm. A given
    // processor is receiver for as many steps as the number of 0s before the
    // first 1 in the binary expansion of its rank.

    // At the kth step, (k starts from 0), stepsize is 2^k. processor i + 2^k
    // sends forest F(k,i+2^k) to processor i. Moreover the processors between
    // i + 1 and i + 2^k send their edges between F(k,i+2^k) and F(k,i).
    // Processor i merges these edges to form F(k+1,i).

    // At the kth step, a processor is responsible for stepSize*numRows number
    // of rows.

	for (int stepSize = 1, rank = procRank; stepSize*numRows < N; stepSize <<= 1, rank >>= 1) {
		
        // This is when process stops being a receiver. A processor sends edges
        // as many times as it has 1 in its rank's binary representation.
        if (rank & 1) {
			receiver = 0;

            // Processor i + 2^k sends its forest to processor i if i is a receiver.
			if (procRank % stepSize == 0) {
				send_forest(procRank, stepSize, forest, forestSize);
			}

            // the processors between i + 1 and i + 2^k send their edges between
            // F(k,i+2^k) and F(k,i). First they compute the spanning forest of
            // their respective edge set. This corresponds to a bipartite forest.
			send_bipartite_forest(procRank, stepSize, squareSize, numRows, adj, N);

		} else if (receiver) {
			int subnumEdges = stepSize*(numRows*(stepSize+1)-1)  + 2*(stepSize*numRows-1);
			int numMaxEdges = (M < subnumEdges) ? M : subnumEdges;
			Edge* edges = (Edge* )calloc(numMaxEdges, sizeof(Edge));
			int numNewEdges = receive_new_edges(procRank, numProcs, stepSize, squareSize,
                numRows, edges, numMaxEdges, forest, forestSize);
			forest = (Edge* )realloc(forest, (stepSize*numRows*2 - 1) * sizeof(Edge));

            // Create a spanning forest of all the edges received.
			forestSize = union_find(edges, N, numNewEdges, forest);
			free(edges);
		}
	}
	if (procRank == 0) {
		memcpy(tree, forest, sizeof(Edge) * forestSize);
	}
	free(forest);
}

/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 */
void computeMST(
		int N,
		int M,
		int *adj)
{
	int procRank, numProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // calloc additionally initialises space with zeros
	Edge* tree = (Edge* )calloc(N-1, sizeof(Edge));

	parallel_boruvka(procRank, numProcs, adj, N, M, tree);

	if (procRank == 0) {
		print_tree(tree, N-1);
	}
}
