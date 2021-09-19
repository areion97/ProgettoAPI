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
    // TODO usare una matrice ottimizza?
    struct AdjList* array;
};

struct HeapNode
{
    int  v;
    unsigned int dist;
    int graphId;
};

struct Heap
{
    int size;
    int *pos;
    struct HeapNode **array;
};

typedef struct Graph* Graph;
typedef struct AdjListNode* AdjListNode;
typedef struct AdjList* AdjList;
typedef struct HeapNode* HeapNode;
typedef struct Heap* Heap;


void topK(Heap maxHeap);
void updateRanking(Heap maxHeap, int k, int el, int graphId);
void printGraph(Graph graph, int graphId);
void printArr(unsigned int *dist, int n);
void printHeap(Heap heap);

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

HeapNode createHeapNode(int v, unsigned int dist, int graphId) {
    HeapNode minHeapNode = (HeapNode) malloc(sizeof(struct HeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    minHeapNode->graphId = graphId;
    return minHeapNode;
}

Heap createHeap(int capacity) {
    Heap minHeap = (Heap) malloc(sizeof(struct Heap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->array = (HeapNode*) malloc(capacity *sizeof(HeapNode));
    return minHeap;
}
void swapHeapNode(HeapNode *a, HeapNode *b) {
    HeapNode t = *a;
    *a = *b;
    *b = t;
}
void decreaseKey(Heap minHeap, int v, unsigned int dist) {
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

void maxHeapify(Heap heap, int n, int i)
{
    int parent = (i - 1) / 2;

    if (heap->array[parent] > 0) {

        if (heap->array[i] > heap->array[parent]) {
            swapHeapNode(&heap->array[i], &heap->array[parent]);
            maxHeapify(heap, n, n-1);
        }
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

    for(int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for(int i = n - 1; i > 0; i--) {
        swapHeapNode(&arr[0], &arr[i]);
        heapify(arr, i, 0);
    }
}
void addChild(Heap heap, HeapNode newNode) {
    heap->array[0] = newNode;
    heapSort(heap->array,heap->size);
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


int isEmpty(Heap heap) {
    return heap->size == 0;
}
HeapNode extractMin(Heap minHeap) {
    if(isEmpty(minHeap))
        return NULL;

    HeapNode root = minHeap->array[0];

    HeapNode lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;

    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}



int isInHeap(Heap minHeap, int v) {
    if (minHeap->pos[v] < minHeap->size)
        return 1;
    return 0;
}

int dijkstraAlgorithm(Graph graph, int src) {
    int sum = 0;

    int V = graph->V;

    unsigned int dist[V];
    Heap minHeap = createHeap(V);


    for (int v = 0; v < V; ++v) {
        dist[v] = INT_MAX;
        minHeap->array[v] = createHeapNode(v, dist[v],0);
        minHeap->pos[v] = v;
    }

    minHeap->array[src] = createHeapNode(src, dist[src],0);
    minHeap->pos[src] = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    minHeap->size = V;

    while (!isEmpty(minHeap)) {

        HeapNode minHeapNode = extractMin(minHeap);

        int u = minHeapNode->v;

        AdjListNode pCrawl = graph->array[u].head;
        while (pCrawl != NULL) {
            int v = pCrawl->dest;

            if (isInHeap(minHeap, v) == 1 && dist[u] != INT_MAX && pCrawl->weight + dist[u] < dist[v]) {
                dist[v] = dist[u] + pCrawl->weight;

                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }
        if(dist[u] != INT_MAX)
        sum +=dist[u];

    }
   // printArr(dist,minHeap->size);
    return sum;
}


void substr(char * string, char * newString, int start, int end) {

    int pos = 0;
    for(int i = start; i < end; i++) {
        newString[pos] = string[i];
        pos++;
    }
}

int main(int argc, char * argv[]) {

    FILE *fp;
    int d = 0;
    int k = 0;
    int graphId = 0;

    char line[100];
    int flag = 0;
    fp = fopen("test/input_5", "r");
    fscanf(fp, "%d %d", &d, &k);

    Heap maxHeap = createHeap(k);

    while((fscanf(fp,"%s\n",&line)) != EOF) {

        if(strcmp(line,"AggiungiGrafo") == 0) {
            Graph graph = createGraph(d);

            for(int i = 0; i < d ; i++) {
                fscanf(fp,"%s\n",&line);
                int start = 0;
                int weight = 0;
                int destGraph = 0;
                for(int j = 0; j <= strlen(line); j++) {
                    if(line[j] == ',' || line[j] == '\000') {
                        char newString[j - start];
                        substr(line, newString, start, j);
                        weight = atoi(newString);
                        start = j + 1;

                        if (weight != 0) {
                            addEdge(graph, i, destGraph, weight);
                        }
                        destGraph++;

                    }
                }
            }

            int minPathValue = dijkstraAlgorithm(graph,0);
            // printHeap(maxHeap);
             printf("I am trying to insert:\n%d:  %d\n",graphId, minPathValue);

            if(flag == 0 && k == maxHeap->size) {
                heapSort(maxHeap->array, k);
                printf("Heap sort!\n");
                flag = 1;
            }

            updateRanking(maxHeap,k,minPathValue,graphId);

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
    if(maxHeap->size == 0) {
        printf("\n");

    }
    else {
        for(int i = 0; i < maxHeap->size; i++) {
            printf("%d ",maxHeap->array[i]->graphId);
        }
    }
}



void updateRanking(Heap maxHeap, int k, int el, int graphId) {

    HeapNode newHeapNode = NULL;

    if(maxHeap->size > k) {
        return;
    }

    if(maxHeap->size == k) {

        HeapNode root = maxHeap->array[0];

        if(el < root->dist) {
            newHeapNode = createHeapNode(maxHeap->size,el,graphId);
            addChild(maxHeap,newHeapNode);
        } else {
            return;
        }
    } else if(maxHeap->size < k) {
        newHeapNode = createHeapNode(maxHeap->size,el,graphId);
        maxHeap->size++;
        maxHeap->array[maxHeap->size-1] = newHeapNode;
    }
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

void printArr(unsigned int *dist, int n) {
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < n; i++)
        printf("%d \t\t %d\n", i, dist[i]);
}

void printHeap(Heap heap) {
    printf("\nHeap Representation\n");
    for (int i = 0; i < heap->size; i++)
        printf("%d \t\t %d\n", heap->array[i]->graphId, heap->array[i]->dist);
}