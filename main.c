
#include <stdio.h>
#include <stdlib.h>
 
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
    int minHeap;
    struct AdjList* array;
};


struct MinHeapNode
{
    int  v;
    int dist;
};
 
 typedef struct MinHeapNode* MinHeapNode;

struct MinHeap
{
     
    // Number of heap nodes present currently
    int size;    
   
    // Capacity of min heap
    int capacity; 
   
    // This is needed for decreaseKey()
    int *pos;   
    MinHeapNode *array;
};
 
 
 
typedef struct Graph* Graph;
typedef struct AdjListNode* AdjListNode;
typedef struct AdjList* AdjList;
typedef struct MinHeap* MinHeap;

 
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
    graph->minHeap = 0;
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


MinHeapNode newMinHeapNode(int v, int dist) {
    MinHeapNode minHeapNode = (MinHeapNode) malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}
 
MinHeap createMinHeap(int capacity) {
	MinHeap minHeap = (MinHeap) malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (MinHeapNode*) malloc(capacity *sizeof(MinHeapNode));
    return minHeap;
}

 
// Stampa il max Heap separato da spazi
void topK(int k) {
}


// TODO strutture dati : maxHeap 
//  dove nel root hai il peggior cammino quindi 
//confrontandolo con quello da aggiungere
// valuti se inserire o meno
void updateRanking(int el) {
	
//confronti e nel caso + insert e heapify
	
}

void printGraph(Graph graph, int graphId);

void decreaseKey(MinHeap minHeap, int v, int dist) {
    // Get the index of v in  heap array
    int i = minHeap->pos[v];
 
    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;
 
    // Travel up while the complete
    // tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i-1)/2;
        minHeap->pos[minHeap->array[(i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
 
        // move to parent index
        i = (i - 1) / 2;
    }
}


// A standard function to
// heapify at given idx
// This function also updates
// position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(MinHeap minHeap,int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;
 
    if (left < minHeap->size &&
        minHeap->array[left]->dist <
         minHeap->array[smallest]->dist )
      smallest = left;
 
    if (right < minHeap->size &&
        minHeap->array[right]->dist <
         minHeap->array[smallest]->dist )
      smallest = right;
 
    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        MinHeapNode smallestNode =
             minHeap->array[smallest];
        MinHeapNode idxNode =
                 minHeap->array[idx];
 
        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;
 
        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest],
                         &minHeap->array[idx]);
 
        minHeapify(minHeap, smallest);
    }
}


int isInMinHeap(MinHeap minHeap, int v) {
   if (minHeap->pos[v] < minHeap->size)
     return 1;
   return 0;
}
 
void printArr(int dist[], int n) {
	printf("Vertex   Distance from Source\n");
	for (int i = 0; i < n; ++i)
		printf("%d \t\t %d\n", i, dist[i]);
}


MinHeapNode extractMin(MinHeap minHeap)
{
    if (isEmpty(minHeap))
        return NULL;
 
    // Store the root node
    MinHeapNode root = minHeap->array[0];
 
    // Replace root node with last node
    MinHeapNode lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;
 
    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;
 
    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);
 
    return root;
}

int isEmpty(MinHeap minHeap)
{
    return minHeap->size == 0;
}

// A utility function to swap two
// nodes of min heap.
// Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a,
                     struct MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// TODO Questa funzione deve restituire un intero che indichi il cammino minimo per quel grafo
int calculateMinHeap(Graph graph, int src) {
     
    // Get the number of vertices in graph
    int V = graph->V;
   
    // dist values used to pick
    // minimum weight edge in cut
    int dist[V];    
 
    // minHeap represents set E
    MinHeap minHeap = createMinHeap(V);
 
    // Initialize min heap with all
    // vertices. dist value of all vertices
    for (int v = 0; v < V; ++v) {
        dist[v] = INT_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }
 
    // Make dist value of src vertex
    // as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]);
    minHeap->pos[src]   = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);
 
    // Initially size of min heap is equal to V
    minHeap->size = V;
 
    // In the followin loop,
    // min heap contains all nodes
    // whose shortest distance
    // is not yet finalized.
    while (!isEmpty(minHeap)) {
        // Extract the vertex with
        // minimum distance value
        MinHeapNode minHeapNode = extractMin(minHeap);
       
        // Store the extracted vertex number
        int u = minHeapNode->v;
 
        // Traverse through all adjacent
        // vertices of u (the extracted
        // vertex) and update
        // their distance values
        AdjListNode pCrawl = graph->array[u].head;
        while (pCrawl != NULL) {
            int v = pCrawl->dest;
 
            // If shortest distance to v is
            // not finalized yet, and distance to v
            // through u is less than its
            // previously calculated distance
            if (isInMinHeap(minHeap, v) == 1 && dist[u] != INT_MAX && pCrawl->weight + dist[u] < dist[v]) {
                dist[v] = dist[u] + pCrawl->weight;
 
                // update distance
                // value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }
    }
 
    // print the calculated shortest distances
    printArr(dist, V);
    return dist;
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
			
			updateRanking(calculateMinHeap(graph,0));
			
			printGraph(graph,graphId);
			graphId++;
			
		}
		else if(strcmp(line,"TopK") == 0) {
			topK(k);
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
 


