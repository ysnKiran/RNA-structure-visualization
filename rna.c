#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

// Check if two bases form a valid base pair
bool isValidBasePair(char base1, char base2) {
    return (base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A') || (base1 == 'C' && base2 == 'G') || (base1 == 'G' && base2 == 'C');
}

// Check if the indices form a kink (base pairs too close)
bool isKink(int idx1, int idx2) {
    return (idx1 >= idx2 - 4);
}

typedef struct {
    int i;
    int j;
} Match;

typedef struct {
    char* sequence;
    int len;
    int maxPairs;
    int** maxMatchings;
    Match* matches;
    int curMatches;
    double execTime;
} Molecule;

// Calculate the maximum number of base pairs using dynamic programming
void calculateMaxBasePairs(Molecule* mol) {
    for (int sz = 6; sz <= mol->len; sz++) {
        for (int i = 0; i + sz - 1 < mol->len; i++) {
            int j = i + sz - 1;
            mol->maxMatchings[i][j] = mol->maxMatchings[i][j - 1];
            for (int t = i; t <= j; t++) {
                if (isValidBasePair(mol->sequence[t], mol->sequence[j]) && !isKink(t, j)) {
                    if (t == i)
                        mol->maxMatchings[i][j] = fmax(mol->maxMatchings[i][j], 1 + mol->maxMatchings[t + 1][j - 1]);
                    else
                        mol->maxMatchings[i][j] = fmax(mol->maxMatchings[i][j], 1 + mol->maxMatchings[i][t - 1] + mol->maxMatchings[t + 1][j - 1]);
                }
            }
        }
    }
    mol->maxPairs = mol->maxMatchings[0][mol->len - 1];
}

// Recursively find the base pairs that form the maximum number of matchings
void findBasePairs(Molecule* mol, int l, int r) {
    if (mol->curMatches == mol->maxPairs || mol->maxMatchings[l][r] == 0)
        return;

    if (isValidBasePair(mol->sequence[l], mol->sequence[r]) && !isKink(l, r)) {
        if (mol->maxMatchings[l + 1][r - 1] + 1 == mol->maxMatchings[l][r]) {
            mol->matches[mol->curMatches].i = l;
            mol->matches[mol->curMatches].j = r;
            mol->curMatches++;
            findBasePairs(mol, l + 1, r - 1);
            return;
        }
    }

    for (int t = l + 1; t <= r; t++) {
        if (isValidBasePair(mol->sequence[t], mol->sequence[r]) && !isKink(t, r)) {
            if (mol->maxMatchings[l][t - 1] + mol->maxMatchings[t + 1][r - 1] + 1 == mol->maxMatchings[l][r]) {
                mol->matches[mol->curMatches].i = t;
                mol->matches[mol->curMatches].j = r;
                mol->curMatches++;
                findBasePairs(mol, l, t - 1);
                findBasePairs(mol, t + 1, r - 1);
                return;
            }
        }
    }

    findBasePairs(mol, l, r - 1);
}

// Set the RNA sequence for the molecule
void setSequence(Molecule* mol, const char* seq) {
    mol->sequence = strdup(seq);
}

// Create a new molecule with the given RNA sequence
Molecule createMolecule(const char* seq) {
    Molecule mol;
    setSequence(&mol, seq);
    mol.len = strlen(seq);

    mol.maxMatchings = (int**)malloc(mol.len * sizeof(int*));
    for (int i = 0; i < mol.len; i++)
        mol.maxMatchings[i] = (int*)calloc(mol.len, sizeof(int));

    clock_t start = clock();
    calculateMaxBasePairs(&mol);
    clock_t stop = clock();
    mol.execTime = (double)(stop - start) / CLOCKS_PER_SEC * 1000;

    mol.matches = (Match*)malloc(mol.maxPairs * sizeof(Match));
    mol.curMatches = 0;
    findBasePairs(&mol, 0, mol.len - 1);

    return mol;
}

// Get the length of the RNA sequence
int getSequenceLength(Molecule* mol) {
    return mol->len;
}

// Get the RNA sequence
char* getSequence(Molecule* mol) {
    return mol->sequence;
}

// Get the maximum number of base pairs
int getMaxBasePairs(Molecule* mol) {
    return mol->maxPairs;
}

// Print the base pairs in a readable format
void printBasePairs(Molecule* mol) {
    for (int i = 0; i < mol->maxPairs; i++)
        printf("%c---%c  (%d, %d)\n", mol->sequence[mol->matches[i].i], mol->sequence[mol->matches[i].j],
               mol->matches[i].i, mol->matches[i].j);
}

// Print the base pairs in a format suitable for writing to a file
void printBasePairsToFile(Molecule* mol) {
    for (int i = 0; i < mol->maxPairs; i++)
        printf("%d %d %c %c\n", mol->matches[i].i, mol->matches[i].j,
               mol->sequence[mol->matches[i].i], mol->sequence[mol->matches[i].j]);
}

// Get the execution time for calculating the maximum base pairs
double getExecutionTime(Molecule* mol) {
    return mol->execTime;
}

// Free the memory allocated for the molecule
void freeMolecule(Molecule* mol) {
    free(mol->sequence);
    for (int i = 0; i < mol->len; i++)
        free(mol->maxMatchings[i]);
    free(mol->maxMatchings);
    free(mol->matches);
}

bool isValidSequence(const char* seq) {
    int len = strlen(seq);
    for (int i = 0; i < len; i++) {
        char base = seq[i];
        if (base != 'A' && base != 'G' && base != 'C' && base != 'U')
            return false;
    }
    return true;
}

int main() {
    char seq[2000];
    printf("Enter the RNA sequence: ");
    scanf("%s", seq);

    if (!isValidSequence(seq)) {
        printf("Invalid RNA sequence. Only A, G, C, and U are allowed.\n");
        return 1;
    }

    // Open the file for writing
    FILE* dataFile = fopen("data.txt", "w");
    if (dataFile == NULL) {
        printf("Error opening file data.txt\n");
        return 1;
    }

    // Print to both console and file
    printf("------------------------------------------------------------------------------------------------------------------------\n\n");
    
    Molecule rnaSequence = createMolecule(seq);
    
    fprintf(dataFile, "%s\n", getSequence(&rnaSequence));
    printf("RNA sequence:                      %s\n", getSequence(&rnaSequence));
    
    printf("Length of the RNA sequence:        %d\n", getSequenceLength(&rnaSequence));
    
    fprintf(dataFile, "%d\n", getMaxBasePairs(&rnaSequence));
    printf("Maximum number of base pairs:      %d\n", getMaxBasePairs(&rnaSequence));
    
    printf("Time taken to calculate max pairs: %.3fms\n", getExecutionTime(&rnaSequence));
    
    printf("The base pairs are:\n");
    
    // Print base pairs to both console and file
    for (int i = 0; i < rnaSequence.maxPairs; i++) {
        printf("%c---%c  (%d, %d)\n", rnaSequence.sequence[rnaSequence.matches[i].i], rnaSequence.sequence[rnaSequence.matches[i].j],
               rnaSequence.matches[i].i, rnaSequence.matches[i].j);
    }
    
    printf("------------------------------------------------------------------------------------------------------------------------\n\n");

    // Write base pairs to the data file
    for (int i = 0; i < rnaSequence.maxPairs; i++) {
        fprintf(dataFile, "%d %d %c %c\n", rnaSequence.matches[i].i, rnaSequence.matches[i].j,
                rnaSequence.sequence[rnaSequence.matches[i].i], rnaSequence.sequence[rnaSequence.matches[i].j]);
    }
    
    FILE* fp = fopen("rna_secondary_structure.txt", "w");
    fprintf(fp, "%s\n", getSequence(&rnaSequence));
    fprintf(fp, "%d\n", getMaxBasePairs(&rnaSequence));
    printBasePairsToFile(&rnaSequence);
    fclose(fp);

    freeMolecule(&rnaSequence);

    // Close the data file
    fclose(dataFile);

    return 0;
}