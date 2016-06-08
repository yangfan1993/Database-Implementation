#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <malloc/malloc.h>
#define MAXINTEGER 2147483647
#include "random.h"

int main(int argc, char* argv[]) {
    if (argc <= 3) {
        printf("No enough arguments!");
        return 1;
    }
    // set the inputs
    int K = atoi(*(argv+1));
    int P = atoi(*(argv+2));
    int L = argc - 3;
    printf("K is %d, P is %d, L is %d \n",K,P,L);
    int *fanout = malloc(L * 4);
    for (int i=0; i<L; i++) {
        *(fanout + i) = atoi(*(argv+3+i));
    }
    // the random generate of inputs (tree keys) and probes
    rand32_t *gen = rand32_init(time(NULL));
    int32_t *input = generate_sorted_unique(K, gen);
    for (int i=0; i<K; i++) {
        printf("%d ", *(input+i));
    }
    printf("\n");
    int32_t *probe = generate_sorted_unique(P, gen);
    free(gen);
    
    clock_t t1 = clock();
    int *length = calloc(L, 4);
    int *offset = calloc(K, 4);
    int *layer = calloc(K, 4);
	int x = 0;
	int y = L-1;
    // to identify the location of each entity
    if (L == 1 && K > (*fanout)-1) {
        printf("The input capacity exceed the capacity of this structure.\n");
        return 1;
    }
	for (int i=0; i<K; i++) {
		*(offset+i) = x;
		x ++;
		*(layer+i) = y;
		(*(length+y))++;
		if (y == L-1) {
			if (*(length+y) % (*(fanout+y)-1) !=0) {
                continue;
			}
			int temp = *(length+y) / (*(fanout+y)-1);
			y --;
			while (temp % (*(fanout+y)) == 0) {
                temp = temp / *(fanout+y);
                y --;
                // the case that the input number is too large
                if (y < 0 && i != K-1) {
				    printf("The input capacity exceed the capacity of this structure.\n");
                    return 1;
                }
                if (i == K-1) {
                    break;
                }
			}
			x = *(length+y);
		} else {
			y = L-1;
			x = *(length+y);
		}
	}
	// the case that the input number is too small
	if (*length == 0) {
		printf("The input capacity could not build such a tree.\n");
		return 1;
	}
    int *elength = calloc(L, 4);
    for (int i=0; i<L; i++) {
        *(elength+i) = *(length+i);
    }
	// to consider the padding of several number of MAXINTEGER
	for (int i=L-1; i>=0; i--) {
        if (*(length+i) % (*(fanout+i)-1) != 0) {
        	*(length+i) = (*(length+i)/(*(fanout+i)-1)+1) * (*(fanout+i)-1);
        }
        if (i != L-1) {
            // consider the situation that there is node in downer layer, but lack of parent node
        	int temp1 = *(length+i)/(fanout[i]-1);
        	int temp2 = *(length+i+1)/(fanout[i+1]-1);
        	if (temp2 > temp1 * fanout[i]) {
        		*(length+i) = ((*(length+i+1)/(*(fanout+i+1)-1))/ *(fanout+i)+1) * (*(fanout+i)-1);
        	}
        }
	}
    // double check the padding to avoid edge case
    for (int i=0; i<L-1; i++) {
        int div1 = *(elength+i) / (*(fanout+i)-1);
        int mod1 = *(elength+i) - div1 * (*(fanout+i)-1);
        if (mod1 == 0) {
            div1 = div1-1;
            mod1 = *(fanout+i)-1;
        }
        if (*(length+i+1) <= (div1 * (*(fanout+i)) + mod1) * (*(fanout+i+1)-1)) {
            *(length+i+1) = (div1 * (*(fanout+i)) + mod1 + 1) * (*(fanout+i+1)-1);
        }
    }
	// prepare storage space for the arrays
	void *(*p) = calloc(L, 8);
	int32_t *(*ptr) = calloc(L, 8); 
	for (int i=0; i<L; i++) {
		posix_memalign(p+i, 16, (*(length+i))*4);
        *(ptr+i) = (int32_t *) *(p+i);
	}
    printf("The storage distribution is ok\n");
	// initialize the arrays and fill in them
	for (int i=0; i<L; i++) {
        for (int j=0; j<*(length+i); j++) {
        	*(*(ptr+i) + j) = MAXINTEGER;
        }
    }
	for (int i=0; i<K; i++) {
        *(*(ptr+layer[i]) + *(offset+i)) = *(input+i);
	}
    clock_t t2 = clock();
    free(layer);
    free(offset);
	// output the array entities
    printf("Output tree structure:\n");
    for (int i=0; i<L; i++) {
        printf("layer%d: ", i);
        for (int j=0; j<*(length+i); j++) {
        	printf("%d ", *(*(ptr+i) + j));
        }
        printf("\n");
    }
    
    clock_t t3 = clock();
    int *output = malloc(P * 4);
    // use binary search to find the corresponding position
    for (int j=0; j<P; j++) {
        int32_t key = *(probe+j);
        int result = 0;
        int start = 0;
        int end = (*fanout)-1;
        for (int i=0; i<L; i++) {
            int32_t *pt = *(ptr+i);
            int pos;
            if (*(pt+start) >= key) {
                pos = start;
            } else if (*(pt+end-1) < key) {
                pos = end;
            } else {
                while (start <= end) {
                	pos = (start + end)/2;
                	if (*(pt+pos) >= key && key > *(pt+(pos-1))) {
                		break;
                	}
                	else if (key > *(pt+pos)){
                        start = pos + 1;
                	} else {
                		end = pos - 1;
                	}
                }
            }
            result += pos;
            if (i != L-1) {
                start = result * (*(fanout+i+1)-1);
                end = (result+1) * (*(fanout+i+1)-1);
            }
        }
        // store the result;
        *(output+j) = result;
	}
    clock_t t4 = clock();
    free(input);
    free(length);
    free(fanout);
    
    // output the result
    printf("Probe and results:\n");
    for (int i=0; i<P; i++) {
        printf("%d %d\n", *(probe+i), *(output+i));
    }
    double clock1 = (double)(t2-t1);
    double clock2 = (double)(t4-t3);
    double time1 = (clock1/CLOCKS_PER_SEC)*1000;
    double time2 = (clock2/CLOCKS_PER_SEC)*1000;
    printf("it takes %f clocks to generate the tree.\n", clock1);
    printf("it takes %f clocks to probe the tree.\n", clock2);
    printf("it takes %f ms to generate the tree.\n", time1);
    printf("it takes %f ms to probe the tree.\n", time2);
    free(output);
    free(probe);
	return 0;
}
