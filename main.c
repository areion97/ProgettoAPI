
#include <stdio.h>
#include <stdlib.h>
 
 
// Ogni nodo del grafo possiede una lista di adiacenza
// Servono V liste di adiacenza dove V è il numero di nodi
// Classifica dei k migliori grafi

struct AdjListNode
{
    int dest;
    int weight;
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

typedef struct Graph* Graph;
typedef struct AdjListNode* AdjListNode;
typedef struct AdjList* AdjList;


AdjListNode newAdjListNode(int dest, int weight) {
    AdjListNode newNode = (AdjListNode) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
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
 
void addEdge(Graph graph, int src, int dest, int weight) {

    AdjListNode newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
 
}

 
void printGraph(Graph graph, int graphId);



// Top k stampa l'indice dei grafi che sono in classifica separati da uno spazio (non necessariamente in ordine)
void topK() {
	
}


int main(int argc, char * argv[]) {
	
	FILE *fp;
	int d = 0;
	int k = 0;
	int graphId = 0;

	char line[100];
	char input[20];
	
	fp = fopen("b.txt", "r");
	fscanf(fp, "%d,%d", &d, &k);
		
	while((fscanf(fp,"%s\n",&line)) != EOF) {

		if(strcmp(line,"AggiungiGrafo") == 0) {
			Graph graph = createGraph(d);
			
			for(int i = 0; i < d ; i++) {
				fscanf(fp,"%s\n",&line);
				for(int j = 0; j < d; j++) {

					int weight = line[2*j]-'0';

					if(weight != 0) {
				    	addEdge(graph, i, j, weight);
					}
				}
			}
			printGraph(graph,graphId);
			graphId++;
			// DONE riuscire a stampare i grafi con le matrici di adiacenza
			
			
			// Funzione che compara due grafi
			// Struttura dati che compara in modo efficiente i grafi
			// Guardi in base a k e ai grafi esistenti in memoria se il nuovo grafo sostituisce come classifica un precedente
			
			
		}
		
		else if(strcmp(line,"TopK") == 0) {
			topK();
		}
	}
	
	fclose(fp);
	
	return 0;
}




void printGraph(Graph graph, int graphId) {
    int v;	
    for (v = 0; v < graph->V; ++v) {
        AdjListNode pCrawl = graph->array[v].head;
        printf("\nGrafo %d \n Vertex %d\n head ",graphId,v);
        while (pCrawl) {
            printf("-> %d", pCrawl->dest);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}
 


