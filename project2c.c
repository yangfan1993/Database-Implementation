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
    if (L != 3) {
        printf("It is not the hardcoded tree\n");
        return 1;
    }
    int *fanout = malloc(L * 4);
    for (int i=0; i<L; i++) {
        int fout = atoi(*(argv+3+i));
        if ((i == 0 && fout!=9) || (i == 1 && fout!=5) || (i == 2 && fout!=9)) {
            printf("It is not the hardcoded tree\n");
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
        //if (mod1 == 0) {
        //    div1 = div1-1;
        //    mod1 = *(fanout+i)-1;
        //}
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
    __m128i *pt0 = (__m128i *) (*ptr);
    __m128i temp1 = _mm_load_si128(pt0);
    __m128i temp2 = _mm_load_si128(pt0+1);
    int *output = malloc(P * 4);
    __m128i *temp = (__m128i *) probe;
    int *out = output;
    // use binary search to find the corresponding position
    for (int j=0; j<P/4; j++) {
        __m128i key = _mm_load_si128(temp++);
        __m128i k1 = _mm_shuffle_epi32(key, _MM_SHUFFLE(0,0,0,0));
        __m128i k2 = _mm_shuffle_epi32(key, _MM_SHUFFLE(1,1,1,1));
        __m128i k3 = _mm_shuffle_epi32(key, _MM_SHUFFLE(2,2,2,2));
        __m128i k4 = _mm_shuffle_epi32(key, _MM_SHUFFLE(3,3,3,3));
        // access layer 0
        __m128i cmp0_11 = _mm_cmpgt_epi32(k1, temp1);
        __m128i cmp0_21 = _mm_cmpgt_epi32(k1, temp2);
        __m128i cmp0_12 = _mm_cmpgt_epi32(k2, temp1);
        __m128i cmp0_22 = _mm_cmpgt_epi32(k2, temp2);
        __m128i cmp0_13 = _mm_cmpgt_epi32(k3, temp1);
        __m128i cmp0_23 = _mm_cmpgt_epi32(k3, temp2);
        __m128i cmp0_14 = _mm_cmpgt_epi32(k4, temp1);
        __m128i cmp0_24 = _mm_cmpgt_epi32(k4, temp2);
        __m128i cmp0_1 = _mm_packs_epi32(cmp0_11, cmp0_21);
        __m128i cmp0_2 = _mm_packs_epi32(cmp0_12, cmp0_22);
        __m128i cmp0_3 = _mm_packs_epi32(cmp0_13, cmp0_23);
        __m128i cmp0_4 = _mm_packs_epi32(cmp0_14, cmp0_24);
        cmp0_1 = _mm_packs_epi16(cmp0_1, _mm_setzero_si128());
        cmp0_2 = _mm_packs_epi16(cmp0_2, _mm_setzero_si128());
        cmp0_3 = _mm_packs_epi16(cmp0_3, _mm_setzero_si128());
        cmp0_4 = _mm_packs_epi16(cmp0_4, _mm_setzero_si128());
        int mask0_1 = _mm_movemask_epi8(cmp0_1);
        int mask0_2 = _mm_movemask_epi8(cmp0_2);
        int mask0_3 = _mm_movemask_epi8(cmp0_3);
        int mask0_4 = _mm_movemask_epi8(cmp0_4);
        int res0_1 = _bit_scan_forward(mask0_1^0x1FF);
        int res0_2 = _bit_scan_forward(mask0_2^0x1FF);
        int res0_3 = _bit_scan_forward(mask0_3^0x1FF);
        int res0_4 = _bit_scan_forward(mask0_4^0x1FF);
        // access layer 1
        __m128i *pt1 = (__m128i *) (*(ptr+1));
        __m128i temp3 = _mm_load_si128(pt1+res0_1);
        __m128i temp4 = _mm_load_si128(pt1+res0_2);
        __m128i temp5 = _mm_load_si128(pt1+res0_3);
        __m128i temp6 = _mm_load_si128(pt1+res0_4);
        __m128i cmp1_1 = _mm_cmpgt_epi32(k1, temp3);
        __m128i cmp1_2 = _mm_cmpgt_epi32(k2, temp4);
        __m128i cmp1_3 = _mm_cmpgt_epi32(k3, temp5);
        __m128i cmp1_4 = _mm_cmpgt_epi32(k4, temp6);
        int mask1_1 = _mm_movemask_ps((__m128)cmp1_1);
        int mask1_2 = _mm_movemask_ps((__m128)cmp1_2);
        int mask1_3 = _mm_movemask_ps((__m128)cmp1_3);
        int mask1_4 = _mm_movemask_ps((__m128)cmp1_4);
        int res1_1 = _bit_scan_forward(mask1_1^0x1F);
        int res1_2 = _bit_scan_forward(mask1_2^0x1F);
        int res1_3 = _bit_scan_forward(mask1_3^0x1F);
        int res1_4 = _bit_scan_forward(mask1_4^0x1F);
        res1_1 += (res0_1 << 2) + res0_1;
        res1_2 += (res0_2 << 2) + res0_2;
        res1_3 += (res0_3 << 2) + res0_3;
        res1_4 += (res0_4 << 2) + res0_4;
        // access layer 2
        __m128i *pt2 = (__m128i *) (*(ptr+2));
        __m128i temp7 = _mm_load_si128(pt2+(res1_1<<1));
        __m128i temp8 = _mm_load_si128(pt2+(res1_1<<1)+1);
        __m128i temp9 = _mm_load_si128(pt2+(res1_2<<1));
        __m128i temp10 = _mm_load_si128(pt2+(res1_2<<1)+1);
        __m128i temp11 = _mm_load_si128(pt2+(res1_3<<1));
        __m128i temp12 = _mm_load_si128(pt2+(res1_3<<1)+1);
        __m128i temp13 = _mm_load_si128(pt2+(res1_4<<1));
        __m128i temp14 = _mm_load_si128(pt2+(res1_4<<1)+1);
        __m128i cmp2_11 = _mm_cmpgt_epi32(k1, temp7);
        __m128i cmp2_21 = _mm_cmpgt_epi32(k1, temp8);
        __m128i cmp2_12 = _mm_cmpgt_epi32(k2, temp9);
        __m128i cmp2_22 = _mm_cmpgt_epi32(k2, temp10);
        __m128i cmp2_13 = _mm_cmpgt_epi32(k3, temp11);
        __m128i cmp2_23 = _mm_cmpgt_epi32(k3, temp12);
        __m128i cmp2_14 = _mm_cmpgt_epi32(k4, temp13);
        __m128i cmp2_24 = _mm_cmpgt_epi32(k4, temp14);
        __m128i cmp2_1 = _mm_packs_epi32(cmp2_11, cmp2_21);
        __m128i cmp2_2 = _mm_packs_epi32(cmp2_12, cmp2_22);
        __m128i cmp2_3 = _mm_packs_epi32(cmp2_13, cmp2_23);
        __m128i cmp2_4 = _mm_packs_epi32(cmp2_14, cmp2_24);
        cmp2_1 = _mm_packs_epi16(cmp2_1, _mm_setzero_si128());
        cmp2_2 = _mm_packs_epi16(cmp2_2, _mm_setzero_si128());
        cmp2_3 = _mm_packs_epi16(cmp2_3, _mm_setzero_si128());
        cmp2_4 = _mm_packs_epi16(cmp2_4, _mm_setzero_si128());
        int mask2_1 = _mm_movemask_epi8(cmp2_1);
        int mask2_2 = _mm_movemask_epi8(cmp2_2);
        int mask2_3 = _mm_movemask_epi8(cmp2_3);
        int mask2_4 = _mm_movemask_epi8(cmp2_4);
        int res2_1 = _bit_scan_forward(mask2_1^0x1FF);
        int res2_2 = _bit_scan_forward(mask2_2^0x1FF);
        int res2_3 = _bit_scan_forward(mask2_3^0x1FF);
        int res2_4 = _bit_scan_forward(mask2_4^0x1FF);
        res2_1 += (res1_1 << 3) + res1_1;
        res2_2 += (res1_2 << 3) + res1_2;
        res2_3 += (res1_3 << 3) + res1_3;
        res2_4 += (res1_4 << 3) + res1_4;
        //save the result
        *(out++) = res2_1;
        *(out++) = res2_2;
        *(out++) = res2_3;
        *(out++) = res2_4;
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
