#include <stdio.h>
#include <stdlib.h>
 
 
// Ogni nodo del grafo possiede una lista di adiacenza
// Servono V liste di adiacenza dove V è il numero di nodi
struct AdjListNode
{
    int dest;
    struct AdjListNode* next;
};
 
struct AdjList
{
    struct AdjListNode *head;
};
 
struct Graph
{
    int V;
    struct AdjList* array;
};
 
typedef struct AdjListNode* AdjListNode;

typedef struct AdjList* AdjList;

typedef struct Graph* Graph;

AdjListNode newAdjListNode(int dest) {
    AdjListNode newNode = (AdjListNode) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}
 
Graph createGraph(int V) {
    Graph graph = (Graph) malloc(sizeof(struct Graph));
    graph->V = V;
 
    graph->array = (AdjList) malloc(V * sizeof(struct AdjList));
 
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;
 
    return graph;
}
 
void addEdge(Graph graph, int src, int dest) {

    AdjListNode newNode = newAdjListNode(dest);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
 
}
 
void printGraph(Graph graph);

int main(int argc, char * argv[]) {

    int V = 5;
    Graph graph = createGraph(V);
	
	
	
	int d = 0;
	int k = 0;
	
	//AggiungiGrafo

	
//3,7,42   0
//0,7,2    1
//7,4,3    2

	int dest = 0;
	for(int i = 0; i<d; i++) {
		for(int j = 0; j<d; j++) {
			fscanf("%d,",&dest)
			addEdge(graph, i, dest);
		}
	}
 
    printGraph(graph);
 
    return 0;
}



void printGraph(Graph graph)
{
    int v;
    for (v = 0; v < graph->V; ++v)
    {
        AdjListNode pCrawl = graph->array[v].head;
        printf("\n Adjacency list of vertex %d\n head ", v);
        while (pCrawl)
        {
            printf("-> %d", pCrawl->dest);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}
 