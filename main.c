#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

struct HeapNode {
    int  v;
    unsigned int dist;
    int graphId;
};

struct Heap {
    int size;
    int *pos;
    struct HeapNode **array;
};

typedef struct Graph* Graph;
typedef struct HeapNode* HeapNode;
typedef struct Heap* Heap;


// Allocates memory and initializes Heap Node
HeapNode createHeapNode(int v, unsigned int dist, int graphId) {
    HeapNode minHeapNode = (HeapNode) malloc(sizeof(struct HeapNode));
    minHeapNode->v = v;
    minHeapNode-> dist= dist;
    minHeapNode->graphId = graphId;
    return minHeapNode;
}

// Allocates memory and initializes Max Heap
Heap createHeap(int capacity) {
    Heap minHeap = (Heap) malloc(sizeof(struct Heap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->array = (HeapNode*) malloc(capacity *sizeof(HeapNode));
    return minHeap;
}

// Deallocates memory
void destroyHeap(Heap h) {
    for(int i = 0; i < h->size; i++) {
        free(h->array[i]);
    }
    free(h->pos);
    free(h->array);
    free(h);
}

// Swaps the content of 2 nodes without allocating memory
void swapHeapNode(HeapNode *a, HeapNode *b) {
    HeapNode t = *a;
    *a = *b;
    *b = t;
}

// Auxiliary function for adding a child to Heap
void maxHeapify(HeapNode *arr, int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && arr[l]->dist > arr[largest]->dist)
        largest = l;

    if (r < n && arr[r]->dist > arr[largest]->dist)
        largest = r;

    if (largest != i) {
        swapHeapNode(&arr[i], &arr[largest]);
        maxHeapify(arr, n, largest);
    }
}

// Sort an array of HeapNodes for heapSort
void heapify(HeapNode *arr, int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if(l < n && arr[l]->dist < arr[largest]->dist)
        largest = l;

    if(r < n && arr[r]->dist < arr[largest]->dist)
        largest = r;

    if(largest != i) {
        swapHeapNode(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

// Sorts array of HeapNodes
void heapSort(HeapNode *arr, int n) {
    for(int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for(int i = n - 1; i > 0; i--) {
        swapHeapNode(&arr[0], &arr[i]);
        heapify(arr, i, 0);
    }
}

// Adds a child to Heap and sorts Heap
void addChild(Heap heap, HeapNode newNode) {
    HeapNode temp = heap->array[0];
    heap->array[0] = newNode;
    free(temp);
    maxHeapify(heap->array,heap->size,0);
}

// Decides whether to update or not MaxHeap everytime a new adjaceny matrix is read
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


// Auxiliary function used by dijkstraAlgorithm for calculating minimum distance
unsigned int minDistance(unsigned int *dist, int V, int *sptSet) {
    unsigned int min = INT_MAX, min_index = 0;

    for (int v = 0; v < V; v++) {
        if (sptSet[v] == 0 && dist[v] <= min)
            min = dist[v], min_index = v;
    }
    return min_index;
}

// Implementation of Dijkstra Algorithm
unsigned int dijkstraAlgorithm(int V, int matrix[V][V], int src) {
    unsigned int sum = 0;
    unsigned int dist[V];
    int sptSet[V];
    for (int i = 0; i < V; i++) {
        dist[i] = INT_MAX, sptSet[i] = 0;
        sptSet[V] = 0;
    }
    dist[src] = 0;
    for (int count = 0; count < V-1; count++) {
        unsigned int u = minDistance(dist, V, sptSet);
        sptSet[u] = 1;

        for (int v = 0; v < V; v++) {
            if (!sptSet[v] && matrix[u][v] && dist[u] != INT_MAX && dist[u] + matrix[u][v] < dist[v])
                dist[v] = dist[u] + matrix[u][v];
        }
    }

    for(int i = 0; i<V; i++) {
        if (dist[i] != INT_MAX)
            sum += dist[i];
    }

    return sum;
}

// Moves FILE pointer to the end of the line
void moveBufferToEnd(FILE * in) {
    char c = '\0';
    while(c != '\n' && feof(in) == 0) c = fgetc(in);
}

void topK(FILE *out, Heap maxHeap) {
    if(maxHeap->size == 0) {
        fprintf(out,"\n");
    } else {
        for(int i = 0; i < maxHeap->size; i++) {
            if(i == maxHeap->size-1) fprintf(out,"%d",maxHeap->array[i]->graphId);
            else fprintf(out,"%d ",maxHeap->array[i]->graphId);
        }
        fprintf(out,"\n");
    }
}

void printHeap(Heap heap) {
    printf("\nHeap:\n");
    for (int i = 0; i < heap->size; i++)
        printf("%d \t\t %d\n", heap->array[i]->graphId, heap->array[i]->dist);
}

int main(int argc, char * argv[]) {

    // Initialize file pointers
    FILE *in, *out;
    in = stdin;
    out = stdout;

    // Reads d and k values
    int d = 0;
    int k = 0;
    int dummyResult = fscanf(in,"%d %d\n",&d,&k);
    if(dummyResult == 0) {return 0;}

    int matrix[d][d];
    int NSTRING = d*10 + 3;
    char line[NSTRING];
    int flag = 0;
    int graphId = 0;

    Heap maxHeap = createHeap(k);

    // Loop instruction until file ends
    while(feof(in) == 0) {
        line[0] = getc(in);

        // AggiungiGrafo - Reads adjacency matrix, calculates dijkstra and updates MaxHeap
        if(line[0] == 'A') {
            char * newString;
            moveBufferToEnd(in);
            // Loop rows of adjacency matrix
            for(int i = 0; i < d; i++) {
                char *tmp = fgets(line,NSTRING,in);
                if(tmp == NULL) return 0;
                int col = 0;
                newString = ",";
                // Loops cols of adjacency matrix
                for (const char *ptr = line; *newString == ','; ptr = newString + 1) {
                    matrix[i][col] = (int)strtol(ptr, &newString, 10);
                    col++;
                    if (newString == ptr)
                        break;
                }
            }

            // Calculate dijkstra
            unsigned int minPathValue = dijkstraAlgorithm(d,matrix,0);
            if(flag == 0 && k == maxHeap->size) {
                heapSort(maxHeap->array, k);
                flag = 1;
            }

            // Updates Max Heap
            updateRanking(maxHeap,k,minPathValue,graphId);
            graphId++;
        }
        // Prints topK Graphs contained in Max Heap
        else if(line[0] == 'T') {
            moveBufferToEnd(in);
            topK(out,maxHeap);
        }

    }

    // Cleans memory and close files
    destroyHeap(maxHeap);
    fclose(in);
    fclose(out);
    return 0;
}
