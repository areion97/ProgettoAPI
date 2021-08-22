#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

struct HeapNode
        {
    int  v;
    int dist;
        };

struct Heap
        {
    int size;
    int capacity;

    // This is needed for decreaseKey()
    int *pos;
    struct HeapNode **array;
        };

typedef struct Graph* Graph;
typedef struct AdjListNode* AdjListNode;
typedef struct AdjList* AdjList;
typedef struct HeapNode* HeapNode;
typedef struct Heap* Heap;


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

HeapNode createHeapNode(int v, int dist) {
    HeapNode minHeapNode = (HeapNode) malloc(sizeof(struct HeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

Heap createHeap(int capacity) {
    Heap minHeap = (Heap) malloc(sizeof(struct Heap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (HeapNode*) malloc(capacity *sizeof(HeapNode));
    return minHeap;
}
void swapHeapNode(HeapNode *a, HeapNode *b) {
    HeapNode t = *a;
    *a = *b;
    *b = t;
}
void decreaseKey(Heap minHeap, int v, int dist) {
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
        swapHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

void minHeapify(Heap minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist )
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->dist < minHeap->array[smallest]->dist )
        smallest = right;

    if (smallest != idx) {
        HeapNode smallestNode = minHeap->array[smallest];
        HeapNode idxNode = minHeap->array[idx];

        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        swapHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

void addChild(Heap maxHeap, int child) {
    int parent = (int)child/2;

    if (maxHeap->array[parent]->dist >= maxHeap->array[child-1]->dist ) {
        return;
    }

    if (child != parent) {

        swapHeapNode(&maxHeap->array[child-1], &maxHeap->array[parent]);
        addChild(maxHeap, parent);
    }
}

void heapify(HeapNode *arr, int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && arr[l]->dist < arr[largest]->dist)
        largest = l;

    if (r < n && arr[r]->dist < arr[largest]->dist)
        largest = r;

    if (largest != i) {
        swapHeapNode(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(HeapNode *arr, int n) {

    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for (int i = n - 1; i > 0; i--) {
        swapHeapNode(&arr[0], &arr[i]);
        heapify(arr, i, 0);
    }
}

int isEmpty(Heap heap) {
    return heap->size == 0;
}
HeapNode extractMin(Heap minHeap) {
    if(isEmpty(minHeap))
        return NULL;

    // Store the root node
    HeapNode root = minHeap->array[0];

    // Replace root node with last node
    HeapNode lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}



int isInHeap(Heap minHeap, int v) {
    if (minHeap->pos[v] < minHeap->size)
        return 1;
    return 0;
}

// TODO Questa funzione deve restituire un intero che indichi il cammino minimo per quel grafo
int dijkstraAlgorithm(Graph graph, int src) {

    int V = graph->V;

    // dist values used to pick
    // minimum weight edge in cut
    int dist[V];

    Heap minHeap = createHeap(V);

    // Initialize min heap with all
    // vertices. dist value of all vertices
    for (int v = 0; v < V; ++v) {
        dist[v] = INT_MAX;
        minHeap->array[v] = createHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex
    // as 0 so that it is extracted first
    minHeap->array[src] = createHeapNode(src, dist[src]);
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
        HeapNode minHeapNode = extractMin(minHeap);

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
            if (isInHeap(minHeap, v) == 1 && dist[u] != INT_MAX && pCrawl->weight + dist[u] < dist[v]) {
                dist[v] = dist[u] + pCrawl->weight;

                // update distance
                // value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }
    }

    // print the calculated shortest distances
    //printArr(dist, V);

    int minDistance = 0;

    return dist;
}

void topK(Heap maxHeap);
void updateRanking(Heap maxHeap, int k, int el);
void printGraph(Graph graph, int graphId);
void printArr(int dist[], int n);
void printHeap(Heap heap);

int main(int argc, char * argv[]) {

    FILE *fp;
    int d = 0;
    int k = 0;
    int graphId = 0;

    char line[100];
    int minPathValue = 0;
    int flag = 0;
    fp = fopen("c.txt", "r");
    fscanf(fp, "%d %d", &d, &k);

    Heap maxHeap = createHeap(k);

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

            //int minPathValue = dijkstraAlgorithm(graph,0);

            minPathValue = rand() % 100;

            if(flag == 0 && k == maxHeap->size) {
                heapSort(maxHeap->array, k);
                flag = 1;
            }

            updateRanking(maxHeap,k,minPathValue);

            graphId++;

        }
        else if(strcmp(line,"TopK") == 0) {
            topK(maxHeap);
        }
    }

    fclose(fp);
    return 0;
}

void topK(Heap maxHeap) {
    for(int i = 0; i < maxHeap->size; i++) {
        printf("%d ",maxHeap->array[i]->dist);
    }
}
// @Develop and @Test
void updateRanking(Heap maxHeap, int k, int el) {

    HeapNode newHeapNode = NULL;

    if(maxHeap->size > k) {
        return;
    }

    if(maxHeap->size == k) {

        HeapNode root = maxHeap->array[0];

        printf("El %d is < of root %d ? %d\n\n", el, root->dist, el < root->dist);

        if(el < root->dist) {
            newHeapNode = createHeapNode(maxHeap->size,el);
            maxHeap->size++;
            maxHeap->array[maxHeap->size-1] = newHeapNode;
            addChild(maxHeap, maxHeap->size);
            //Dealloc memory the node has been sostituted
            //maxHeap->array[maxHeap->size-1] = NULL;
        } else {
            return;
        }
    } else if(maxHeap->size < k) {
        newHeapNode = createHeapNode(maxHeap->size,el);
        maxHeap->size++;
        maxHeap->array[maxHeap->size-1] = newHeapNode;
        printf("Inserisco il nodo alla classifica %d\n", el);
    }

    printHeap(maxHeap);

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

void printArr(int dist[], int n) {
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < n; ++i)
        printf("%d \t\t %d\n", i, dist[i]);
}

void printHeap(Heap heap) {
    printf("\nHeap Representation\n");
    for (int i = 0; i < heap->size; i++)
        printf("%d \t\t %d\n", heap->array[i]->v, heap->array[i]->dist);
}

