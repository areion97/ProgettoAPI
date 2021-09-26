#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
typedef struct AdjListNode* AdjListNode;
typedef struct AdjList* AdjList;
typedef struct HeapNode* HeapNode;
typedef struct Heap* Heap;


void topK(Heap maxHeap);
void updateRanking(Heap maxHeap, int k, int el, int graphId);
void printGraph(Graph graph, int graphId);
void printArr(unsigned int *dist, int n);
void printHeap(Heap heap);


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

int minDistance(int *dist, int V, int *sptSet) {
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++) {
        if (sptSet[v] == 0 && dist[v] <= min)
            min = dist[v], min_index = v;
    }
    return min_index;
}

int dijkstraAlgorithm(int V, int matrix[V][V], int src) {
    int sum = 0;
    int dist[V];
    int sptSet[V];
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = 0;

    dist[src] = 0;

    for (int count = 0; count < V-1; count++) {

        int u = minDistance(dist, V, sptSet);
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




void substr(char * string, char * newString, int start, int end) {

    int pos = 0;
    for(int i = start; i < end; i++) {
        newString[pos] = string[i];
        pos++;
    }
}
void printMatrix(int V, int matrix[V][V]) {
    for(int i = 0; i < V; i++) {
        for(int j = 0; j < V; j++) {
            printf("%d, ",matrix[i][j]);
        }
        printf("\n");
    }
}

#define NSTRING 2400
int main(int argc, char * argv[]) {

    FILE *fp;
    int graphId = 0;
    char line[NSTRING];
    int d = 0;
    int k = 0;
    int c;
    int flag = 0;
    fp = fopen("../test/input.5", "r");
    fscanf(fp,"%d %d\n",&d,&k);
    Heap maxHeap = createHeap(k);

    while((fgets(line,NSTRING,fp)) != NULL) {

        if(strcmp(line,"AggiungiGrafo\n") == 0) {
            int matrix[d][d];

            for(int i = 0; i < d ; i++) {
                fgets(line,NSTRING,fp);

                int start = 0;
                int weight = 0;
                int destGraph = 0;

                for(int j = 0; j <= strlen(line); j++) {
                    if(line[j] == ',' || line[j] == '\000') {

                        char newString[j - start];
                        substr(line, newString, start, j);
                        newString[j-start] = '\0';
                        weight = atoi(newString);
                        start = j + 1;

                        if (weight != 0) {
                            matrix[i][destGraph] = weight;

                        }else {
                            matrix[i][destGraph] = 0;
                        }

                        destGraph++;
                    }
                }
            }

             int minPathValue = dijkstraAlgorithm(d,matrix,0);
            //printf("%d\n",minPathValue);
            //printf("\n%d:  %d\n",graphId, minPathValue);
            if(flag == 0 && k == maxHeap->size) {
               heapSort(maxHeap->array, k);
                //printf("Heap sort!\n");
               flag = 1;
          }

            updateRanking(maxHeap,k,minPathValue,graphId);
            graphId++;

        }
        else if(strcmp(line,"TopK\n") == 0) {
            topK(maxHeap);
        }
    }

    fclose(fp);
    return 0;
}

void topK(Heap maxHeap) {
    if(maxHeap->size == 0) {
        printf("\n");
    } else {
        for(int i = 0; i < maxHeap->size; i++) {
            printf("%d ",maxHeap->array[i]->graphId);
        }
        printf("\n");
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

            // TODO
            //Dealloc memory the node has been sostituted
            //maxHeap->array[maxHeap->size-1] = NULL;
        } else {
            return;
        }
    } else if(maxHeap->size < k) {
        newHeapNode = createHeapNode(maxHeap->size,el,graphId);
        maxHeap->size++;
        maxHeap->array[maxHeap->size-1] = newHeapNode;
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
