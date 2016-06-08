// author: Fan Yang(fy2207)
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc/malloc.h>
#include <x86intrin.h>
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
        int fout = atoi(*(argv+3+i));
        if ((fout!=5) && (fout!=9) && (fout!=17)) {
            printf("The current tree only accepts fanout of 5, 9, 17\n");
            return 1;
        }
        *(fanout + i) = fout;
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
                if (y < 0 && i == K-1) {
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
    int *temp = probe;
    // use binary search to find the corresponding position
    for (int j=0; j<P; j++) {
        __m128i key = _mm_cvtsi32_si128(*(temp++));
        key = _mm_shuffle_epi32(key, 0);
        int result = 0;
        for (int i=0; i<L; i++) {
            if (*(fanout+i) == 5) {
                __m128i *pt = (__m128i *) (*(ptr+i)+(result << 2));
                __m128i tempA = _mm_load_si128(pt);
                __m128i cmp = _mm_cmpgt_epi32(key, tempA);
                int mask = _mm_movemask_ps((__m128)cmp);
                int pos = _bit_scan_forward(mask^0x1F);
                result += (result << 2) + pos;
            } else if (*(fanout+i) == 9) {
                __m128i *pt = (__m128i *) (*(ptr+i)+(result << 3));
                __m128i tempA = _mm_load_si128(pt);
                __m128i tempB = _mm_load_si128(pt+1);
                __m128i cmp_1 = _mm_cmpgt_epi32(key, tempA);
                __m128i cmp_2 = _mm_cmpgt_epi32(key, tempB);
                __m128i cmp = _mm_packs_epi32(cmp_1, cmp_2);
                cmp = _mm_packs_epi16(cmp, _mm_setzero_si128());
                int mask = _mm_movemask_epi8(cmp);
                int pos = _bit_scan_forward(mask^0x1FF);
                result += (result << 3) + pos;
            } else if (*(fanout+i) == 17) {
                __m128i *pt = (__m128i *) (*(ptr+i)+(result << 4));
                __m128i tempA = _mm_load_si128(pt);
                __m128i tempB = _mm_load_si128(pt+1);
                __m128i tempC = _mm_load_si128(pt+2);
                __m128i tempD = _mm_load_si128(pt+3);
                __m128i cmp_1 = _mm_cmpgt_epi32(key, tempA);
                __m128i cmp_2 = _mm_cmpgt_epi32(key, tempB);
                __m128i cmp_3 = _mm_cmpgt_epi32(key, tempC);
                __m128i cmp_4 = _mm_cmpgt_epi32(key, tempD);
                __m128i cmp_1_2 = _mm_packs_epi32(cmp_1, cmp_2);
                __m128i cmp_3_4 = _mm_packs_epi32(cmp_3, cmp_4);
                __m128i cmp = _mm_packs_epi16(cmp_1_2, cmp_3_4);
                int mask = _mm_movemask_epi8(cmp);
                int pos = _bit_scan_forward(mask^0x1FFFF);
                result += (result << 4) + pos;
            }
        }
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
