/* Program for fast scoring of citations.

Usage:

cscore [citations] [numcites] [numfeats] [limit] < feature_scores > results

  The [citations] file consists of [numcites] records with the following format
  in C structure notation:

  struct {
    unsigned int pmid; // PubMed ID of citation
    unsigned short nfeatures; // Number of features
    unsigned short features[nfeatures]; // Feature vector
  };

  The feature scores from standard input are a list of [numfeats]
  64-bit doubles.

  The output is a list of [limit] citation scores as score_t structures.
  
                                 
*/

/*
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
*/

#include <stdio.h>
#include <stdlib.h>

// Holds PubMed ID and score of a citation
typedef struct {
    float score;
    unsigned int pmid;
} score_t;

// For qsort, to sort the scores in decreasing order
int compare_scores(const void *a, const void *b) {
    return -( ((const score_t *)a)->score - ((const score_t *)b)->score );
}

int main (int argc, char **argv) {

    // Parameters
    char *cite_filename = argv[1]; // name of citations file
    int numcites = atoi (argv[2]); // number of citations
    int numfeats = atoi (argv[3]); // number of features
    int limit = atoi (argv[4]);  // number of results to return

    // Loop variables
    int pi = 0; // PubMed ID index in loop
    int fi = 0; // Feature index in feature vector
    unsigned short featvec_size = 0; // size of current feature vector
    unsigned short featvec[1000]; // maximum of 1000 features per citation

    // Scores of all citations 
    score_t *scores = (score_t*) malloc (numcites * sizeof(score_t));

    // Scores of all features
    double *featscores = (double*) malloc (numfeats * sizeof(double));

    // Read in feature scores
    fread(featscores, sizeof(double), numfeats, stdin);

    // Calculate citation scores
    FILE *citefile = fopen(cite_filename, "rb");
    for(pi = 0; pi < numcites; pi++) {

        // Read feature vector structure
        fread(&scores[pi].pmid, sizeof(unsigned int), 1, citefile);
        fread(&featvec_size, sizeof(unsigned short), 1, citefile);
        fread(featvec, sizeof(unsigned short), featvec_size, citefile);

        // Calculate citation score
        scores[pi].score = 0.0;
        for(fi = 0; fi < featvec_size; fi++) {
            scores[pi].score += (float)featscores[featvec[fi]];
        }
        
    }
    fclose(citefile);

    // Sort the citations
    qsort((void *)scores, numcites, sizeof(score_t), compare_scores);
    
    // Output the top [limit] citations
    fwrite(scores, sizeof(score_t), limit, stdout);

    return 0;
}
