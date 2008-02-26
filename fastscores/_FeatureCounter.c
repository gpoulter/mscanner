/* Calculate the number of occurrences of each feature within a 
specified date range.

Usage:

./featcounts [citations] [numdocs] [numfeats] [mindate] [maxdate] > scores

  The [citations] file consists of [numcites] records, in the format
  used by FeatureStream.py

  struct {
    unsigned int pmid;     // PubMed ID of citation
    unsigned int date;     // Date of Medline record completion
    unsigned short nbytes; // Number of bytes in encoded feature vector
    char features[nbytes]; // Variable-byte-encoded feature vector
  };

  The output is an array of [numfeats] 32-bit integers, with the number
  of occurrences of each feature in the stream, within the specified date range.
  
                                 
*/

/* This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>. */

#include <stdio.h>
#include <stdlib.h>


// Search sorted array A of length N for needle.
// Return 1 if we find the needle, 0 if we do not
// http://en.wikipedia.org/wiki/Binary_search
int binary_search(unsigned int *A, unsigned int N, unsigned int needle) {
    int low = 0;
    int high = N-1;
    int mid = 0;
    while (low <= high) {
        mid = (low + high) / 2;
        if (A[mid] > needle) {
            high = mid - 1;
        } else if (A[mid] < needle) {
            low = mid + 1;
        } else {
            return 1;
        }
    }
    return 0;
}


int main (int argc, char **argv)
{
    // Parameters
    char *cite_filename = argv[1];             // Name of citations file
    unsigned int numcites = atoi (argv[2]);    // Number of citations
    unsigned int numfeats = atoi (argv[3]);    // Number of features
    unsigned int mindate = atoi (argv[4]);     // Minimum date to consider
    unsigned int maxdate = atoi (argv[5]);     // Maximum date to consider
    unsigned int numexcluded = atoi (argv[6]); // Number of excluded citations
    
    // Loop variables
    FILE *citefile = NULL; // File with citation scores
    unsigned int pi = 0; // Loop variable: number of PubMed ID's so far
    unsigned int fi = 0; // Loop variable: index into feature vector
    unsigned int date = 0; // Date of the current citation
    unsigned int pmid = 0; // PubMed ID of the current citation
    unsigned int ndocs = 0; // Number of documents counted
    unsigned short featvec_size = 0; // Size of current feature vector
    unsigned int featvec[1000]; // Max 1000 features per citation (16 or 32 bit)

    unsigned short featvec_nbytes = 0; // Bytes in encoded feature vector
    unsigned char bytes[4000]; // Bytes of encoded feature vector
    unsigned int gap = 0; // Gap between feature IDs
    unsigned int last = 0; // Value of previous decoded feature ID

    // Allocate space for list of excluded PMIDs 
    unsigned int *excluded = (unsigned int*) malloc (numexcluded * sizeof(int));
    
    // Allocate space for vector of feature counts 
    int *featcounts = (int*) malloc (numfeats * sizeof(int));

    // Read the excluded PMIDs from input
    fread(excluded, sizeof(int), numexcluded, stdin);

    // Initialise the feature counts to zero
    for(fi = 0; fi < numfeats; fi++) featcounts[fi] = 0;

    // Loop for calculating feature counts in the databse
    citefile = fopen(cite_filename, "rb");
    for(pi = 0; pi < numcites; pi++) {
        // Read feature vector from the binary file
        fread(&pmid, sizeof(unsigned int), 1, citefile);
        fread(&date, sizeof(unsigned int), 1, citefile);

        // Decode variable byte encoded feature vector
        fread(&featvec_nbytes, sizeof(unsigned short), 1, citefile);
        fread(bytes, sizeof(unsigned char), featvec_nbytes, citefile);
        gap = 0;
        last = 0;
        featvec_size = 0;
        // Read groups of 7 bits, low to high.  
        // Set high bit set marks end-of-number.
        for(fi = 0; fi < featvec_nbytes; fi++) {
            gap = (gap << 7) | (bytes[fi] & 0x7f);
            if(bytes[fi] & 0x80) {
                last += gap;
                featvec[featvec_size++] = last;
                gap = 0;
            }
        }
        
        // Don't bother if the date is outside the range
        if ((date < mindate) || (date > maxdate)) {
            continue;
        }
        // Don't bother if it is a member of exclude
        if (binary_search(excluded, numexcluded, pmid) == 1) {
            continue;
        }
        // Increment feature counts
        for(fi = 0; fi < featvec_size; fi++)
            featcounts[featvec[fi]]++;
        // Increment document count
        ndocs++;
    }
    fclose(citefile);

    // Print number of docs, feature counts before returning from main
    fwrite(&ndocs, sizeof(int), 1, stdout);
    fwrite(featcounts, sizeof(int), numfeats, stdout);
    return 0;
}
