#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

bool isValidPairing(char base1, char base2) {
    return (base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A') || (base1 == 'C' && base2 == 'G') || (base1 == 'G' && base2 == 'C');
}

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
    int ans;
    int** maxMatchings;
    Match* matches;
    int curMatches;
    double exec_time;
} Molecule;

void calculateMaxMatchings(Molecule* mol) {
    for (int sz = 6; sz <= mol->len; sz++) {
        for (int i = 0; i + sz - 1 < mol->len; i++) {
            int j = i + sz - 1;
            mol->maxMatchings[i][j] = mol->maxMatchings[i][j - 1];
            for (int t = i; t <= j; t++) {
                if (isValidPairing(mol->sequence[t], mol->sequence[j]) && !isKink(t, j)) {
                    if (t == i)
                        mol->maxMatchings[i][j] = fmax(mol->maxMatchings[i][j], 1 + mol->maxMatchings[t + 1][j - 1]);
                    else
                        mol->maxMatchings[i][j] = fmax(mol->maxMatchings[i][j], 1 + mol->maxMatchings[i][t - 1] + mol->maxMatchings[t + 1][j - 1]);
                }
            }
        }
    }
    mol->ans = mol->maxMatchings[0][mol->len - 1];
}

void calculateMatchings(Molecule* mol, int l, int r) {
    if (mol->curMatches == mol->ans || mol->maxMatchings[l][r] == 0)
        return;

    if (isValidPairing(mol->sequence[l], mol->sequence[r]) && !isKink(l, r)) {
        if (mol->maxMatchings[l + 1][r - 1] + 1 == mol->maxMatchings[l][r]) {
            mol->matches[mol->curMatches].i = l;
            mol->matches[mol->curMatches].j = r;
            mol->curMatches++;
            calculateMatchings(mol, l + 1, r - 1);
            return;
        }
    }

    for (int t = l + 1; t <= r; t++) {
        if (isValidPairing(mol->sequence[t], mol->sequence[r]) && !isKink(t, r)) {
            if (mol->maxMatchings[l][t - 1] + mol->maxMatchings[t + 1][r - 1] + 1 == mol->maxMatchings[l][r]) {
                mol->matches[mol->curMatches].i = t;
                mol->matches[mol->curMatches].j = r;
                mol->curMatches++;
                calculateMatchings(mol, l, t - 1);
                calculateMatchings(mol, t + 1, r - 1);
                return;
            }
        }
    }

    calculateMatchings(mol, l, r - 1);
}

void setSequence(Molecule* mol, const char* seq) {
    mol->sequence = strdup(seq);
}

Molecule createMolecule(const char* seq) {
    Molecule mol;
    setSequence(&mol, seq);
    mol.len = strlen(seq);

    mol.maxMatchings = (int**)malloc(mol.len * sizeof(int*));
    for (int i = 0; i < mol.len; i++)
        mol.maxMatchings[i] = (int*)calloc(mol.len, sizeof(int));

    clock_t start = clock();
    calculateMaxMatchings(&mol);
    clock_t stop = clock();
    mol.exec_time = (double)(stop - start) / CLOCKS_PER_SEC * 1000;

    mol.matches = (Match*)malloc(mol.ans * sizeof(Match));
    mol.curMatches = 0;
    calculateMatchings(&mol, 0, mol.len - 1);

    return mol;
}

int getSequenceLength(Molecule* mol) {
    return mol->len;
}

char* getSequence(Molecule* mol) {
    return mol->sequence;
}

int getMaxMatchings(Molecule* mol) {
    return mol->ans;
}

void printPairs(Molecule* mol) {
    for (int i = 0; i < mol->ans; i++)
        printf("%c---%c  (%d, %d)\n", mol->sequence[mol->matches[i].i], mol->sequence[mol->matches[i].j],
               mol->matches[i].i, mol->matches[i].j);
}

void printBases(Molecule* mol) {
    for (int i = 0; i < mol->ans; i++)
        printf("%d %d %c %c\n", mol->matches[i].i, mol->matches[i].j,
               mol->sequence[mol->matches[i].i], mol->sequence[mol->matches[i].j]);
}

double getExecutionTime(Molecule* mol) {
    return mol->exec_time;
}

void freeMolecule(Molecule* mol) {
    free(mol->sequence);
    for (int i = 0; i < mol->len; i++)
        free(mol->maxMatchings[i]);
    free(mol->maxMatchings);
    free(mol->matches);
}

int main() {
    char seq[100];
    printf("Enter the RNA sequence: ");
    scanf("%s", seq);

    printf("------------------------------------------------------------------------------------------------------------------------\n\n");
    Molecule rnaSequence = createMolecule(seq);
    printf("RNA sequence:                      %s\n", getSequence(&rnaSequence));
    printf("Length of the RNA sequence:        %d\n", getSequenceLength(&rnaSequence));
    printf("Maximum number of base pairs:      %d\n", getMaxMatchings(&rnaSequence));
    printf("Time taken to calculate max pairs: %.3fms\n", getExecutionTime(&rnaSequence));
    printf("The base pairs are:\n");
    printPairs(&rnaSequence);
    printf("------------------------------------------------------------------------------------------------------------------------\n\n");

    FILE* fp = fopen("rna_secondary_structure.txt", "w");
    fprintf(fp, "%s\n", getSequence(&rnaSequence));
    fprintf(fp, "%d\n", getMaxMatchings(&rnaSequence));
    printBases(&rnaSequence);
    fclose(fp);

    freeMolecule(&rnaSequence);

    return 0;
}