#include <stdio.h>
#include <stdlib.h>

// Simple tests for ctypes
void double_int(int a, int *b) { *b = a*2; }
void double_array(int len, int *a) { int i;
for(i = 0; i < len; i++) a[i] *= 2; }

// Holds PubMed ID and score of a citation
typedef struct {
    float score;
    unsigned int pmid;
} score_t;

// For qsort, to sort the scores in decreasing order
int compare_scores(const void *a, const void *b) {
    return -( ((const score_t *)a)->score - ((const score_t *)b)->score );
}

void cscore(
    char *cite_filename,  // File to open for citation stream
    int numcites,         // Number of citations
    int numfeats,         // Number of features
    int limit,            // Number of pmid,score pairs to return
    float offset,        // Amount to add to citation score
    double *featscores,   // Array of feature scores
    float *o_scores,      // Output array to store returned scores
    int *o_pmids          // Output array to store returned pmids
    ) 
    {
        
    FILE *citefile = NULL; // File with citation scores

    // Loop variables
    int pi = 0; // PubMed ID index in loop
    int fi = 0; // Feature index in feature vector
    unsigned short featvec_size = 0; // size of current feature vector
    unsigned short featvec[1000]; // maximum of 1000 features per citation

    // Scores of all citations 
    score_t *scores = (score_t*) malloc (numcites * sizeof(score_t));

    // Calculate citation scores
    citefile = fopen(cite_filename, "rb");
    for(pi = 0; pi < numcites; pi++) {

        // Read feature vector structure
        fread(&scores[pi].pmid, sizeof(unsigned int), 1, citefile);
        fread(&featvec_size, sizeof(unsigned short), 1, citefile);
        fread(featvec, sizeof(unsigned short), featvec_size, citefile);

        // Calculate citation score
        scores[pi].score = offset; // start with the offset score
        for(fi = 0; fi < featvec_size; fi++) {
            scores[pi].score += (float)featscores[featvec[fi]];
        }
        
    }
    fclose(citefile);

    // Sort the citations
    qsort((void *)scores, numcites, sizeof(score_t), compare_scores);
    
    // Store top citations in o_scores and o_pmids
    for(pi = 0; pi < limit; pi++) {
        o_scores[pi] = scores[pi].score;
        o_pmids[pi] = scores[pi].pmid;
    }

    // Free the only malloc'd pointer
    free(scores);
}
