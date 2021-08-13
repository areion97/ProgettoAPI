
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
 
 	int i = 0;
    for (; i < V; ++i)
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
	
	FILE *fp;
	int d = 0;
	int k = 0;
	char line[100];
	char input[20];
	fp = fopen("b.txt", "r");
	fscanf(fp, "%d,%d", &d, &k);
	printf("%d %d\n", d, k );
	
	fclose(fp);
	
	Graph graph = createGraph(d);
	while((fscanf(fp,"%s\n",&line)) != EOF) {
	
		if(strcmp(line,"AggiungiGrafo") == 0) {
			int dest = 0;

			fscanf(fp,"%d\n",&dest);
			for(int i = 0; i<d; i++) {
				for(int j = 0; j<d; j++) {
				    addEdge(graph, i, dest);
				}
			}
		
		}
		
		else if(strcmp(line,"TopK") == 0) {
			
		}
	}
	
	fclose(fp);
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
 


