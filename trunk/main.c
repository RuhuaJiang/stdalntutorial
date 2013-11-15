/*
 * 	This is example showing how to use stdaln.c provided by Heng Li.
 *  Thanks Heng Li for his great work on this alignment library.
 *  Created on: Nov 14, 2013
 *  Author: ruhua.jiang@gmail.com
 */
#include <stdio.h>
#include <string.h>
#include "stdaln.h"

#define __gen_ap(par, opt) do {									\
		int i;													\
		for (i = 0; i < 25; ++i) (par).matrix[i] = -(opt)->b;	\
		for (i = 0; i < 4; ++i) (par).matrix[i*5+i] = (opt)->a; \
		(par).gap_open = (opt)->q; (par).gap_ext = (opt)->r;	\
		(par).gap_end = (opt)->r;								\
		(par).row = 5; (par).band_width = opt->bw;				\
	} while (0)

/* char -> 5 (=4+1) nucleotides */
static unsigned char aln_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


typedef struct {
	int a; //Score of a match [1]
	int b; //Mismatch penalty [3]
    int q; //Gap open penalty [5]
	int r; //Gap extension penalty. The penalty for a contiguous gap of size k is q+k*r. [2]
	int bw; //-w INT	 Band width in the banded alignment [33]
}opt_t;

int main()
{
	//STEP1: Given two sequences
	//Input variables
	char *seq1 = strdup("TTGAATCTATTTCTTAGTTCTACATTTTTGTGTGTGTGTTCGTGAAGTTTTCAGGGTTTTCTACAAATAT");
	char *seq2 = strdup("TTG");
	int len1 = strlen(seq1); //seq length of first
	int len2 = strlen(seq2); //seq length of second
	opt_t *opt; //alignment parameters



	//Just for type cast
	unsigned char *seq11, *seq22;
	AlnParam par;

	//Temp variables
	int i,j;
	int  matrix[25];

	//Output variables
	int path_len; //path length
	int score;
	uint32_t * cigar;  //cigar string
	int n_cigar;
	path_t *path;  //For keeping track of path when doing backtrack in dynamic programming table

	int ret, subo, use_global;
    ////Make some parameters
    opt = (opt_t*)malloc(sizeof(opt_t));
	opt->a = 1;
	opt->b = 3;
	opt->q = 5;
	opt->r = 2;
	opt->bw = 33;
	par.matrix = matrix;
	__gen_ap(par, opt);


	//seq1-> seq11, seq2->seq22
	seq11 = (unsigned char*)malloc(sizeof(unsigned char) * len1);
	seq22 = (unsigned char*)malloc(sizeof(unsigned char) * len2);
	for (i = 0; i < len1; ++i)
			seq11[i] = aln_nt4_table[(int)seq1[i]];
	for (j = 0; j < len2; ++j)
			seq22[j] = aln_nt4_table[(int)seq2[j]];

	path = (path_t *)malloc((strlen(seq1)+strlen(seq2))*sizeof(path_t));

	use_global =0;

	//Do banded global alignment, pretty good
	if(use_global) //global alignment, use banded strategy, affine gap
	{
		score = aln_global_core(seq1,len1,seq2,len2, &par, path, &path_len);
		cigar = aln_path2cigar32(path,path_len, &n_cigar);


		//Get some output to see what happen
		fprintf(stderr,"score: %d  path_len: %d \n", score, path_len);
		for( i = 0; i < path_len; i++)
				fprintf(stderr,"i,j (%d,%d) \n",path[i].i,path[i].j);
		for( i = 0; i < n_cigar; i++)
			fprintf(stderr,"%llu \n",cigar[i] );
	}
	else  //local alignment, use banded strategy
	{
		score = aln_local_core(seq11, len1, seq22, len2, &par, path, &path_len, 1, &subo);
		cigar = aln_path2cigar32(path,path_len, &n_cigar);

		for( i = 0; i < path_len; i++)
					fprintf(stderr,"i,j (%d,%d) \n",path[i].i,path[i].j);
		fprintf(stderr,"score: %d  path_len: %d \n", score, path_len);

		if (ret < 0 || subo == ret) { // no hit or tandem hits
			free(path); free(cigar); free(seq1); n_cigar = 0;
			return 0;
		}
	}
}
