#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

// Structures
typedef struct {
    int i;
    int j;
} BasePair;

typedef struct {
    char* sequence;
    int length;
    int max_pairs;
    int** max_matchings;
    BasePair* pairs;
    int current_pairs;
    double execution_time;
} RNA_Structure;

// Function prototypes
void process_rna_sequence(RNA_Structure* rna);
void calculate_max_pairs(RNA_Structure* rna);
void find_pairs(RNA_Structure* rna, int l, int r);
void set_rna_sequence(RNA_Structure* rna, const char* seq);
void create_rna(RNA_Structure* rna, const char* seq); // Change in function prototype
int get_rna_sequence_length(RNA_Structure* rna);
char* get_rna_sequence(RNA_Structure* rna);
int get_max_rna_base_pairs(RNA_Structure* rna);
double get_execution_time_rna(RNA_Structure* rna);
void print_rna_base_pairs(RNA_Structure* rna);
void print_rna_base_pairs_to_file(RNA_Structure* rna);
void free_rna(RNA_Structure* rna);

// Main function
int main() {
    printf("------------------------------------------------------------------------------------------------------------------------\n\n");

    RNA_Structure rna_molecule;
    process_rna_sequence(&rna_molecule);

    printf("RNA sequence:                      %s\n", get_rna_sequence(&rna_molecule));
    printf("Length of the RNA sequence:        %d\n", get_rna_sequence_length(&rna_molecule));
    printf("Maximum number of base pairs:      %d\n", get_max_rna_base_pairs(&rna_molecule));
    printf("Time taken to calculate max pairs: %.3fms\n", get_execution_time_rna(&rna_molecule));
    printf("The base pairs are:\n");
    print_rna_base_pairs(&rna_molecule);
    printf("------------------------------------------------------------------------------------------------------------------------\n\n");

    FILE* file_ptr = fopen("rna_secondary_structure.txt", "w");
    if (file_ptr != NULL) {
        fprintf(file_ptr, "%s\n", get_rna_sequence(&rna_molecule));
        fprintf(file_ptr, "%d\n", get_max_rna_base_pairs(&rna_molecule));
        print_rna_base_pairs_to_file(&rna_molecule);
        fclose(file_ptr);
    } else {
        printf("Error: Unable to open file for writing.\n");
    }

    free_rna(&rna_molecule);

    return 0;
}

// Process RNA sequence
void process_rna_sequence(RNA_Structure* rna) {
    char rna_seq[100];
    printf("Enter the RNA sequence: ");
    scanf("%s", rna_seq);
    create_rna(rna, rna_seq); // Change in function call
}

// Create a new molecule with the given RNA sequence
void create_rna(RNA_Structure* rna, const char* seq) {
    set_rna_sequence(rna, seq);
    rna->length = strlen(seq);

    rna->max_matchings = (int**)malloc(rna->length * sizeof(int*));
    for (int i = 0; i < rna->length; i++)
        rna->max_matchings[i] = (int*)calloc(rna->length, sizeof(int));

    clock_t start = clock();
    calculate_max_pairs(rna);
    clock_t stop = clock();
    rna->execution_time = (double)(stop - start) / CLOCKS_PER_SEC * 1000;

    rna->pairs = (BasePair*)malloc(rna->max_pairs * sizeof(BasePair));
    rna->current_pairs = 0;
    find_pairs(rna, 0, rna->length - 1);
}

// Calculate the maximum number of base pairs using dynamic programming
void calculate_max_pairs(RNA_Structure* rna) {
    for (int size = 6; size <= rna->length; size++) {
        for (int i = 0; i + size - 1 < rna->length; i++) {
            int j = i + size - 1;
            rna->max_matchings[i][j] = rna->max_matchings[i][j - 1];
            for (int t = i; t <= j; t++) {
                if ((rna->sequence[t] == 'A' && rna->sequence[j] == 'U') ||
                    (rna->sequence[t] == 'U' && rna->sequence[j] == 'A') ||
                    (rna->sequence[t] == 'C' && rna->sequence[j] == 'G') ||
                    (rna->sequence[t] == 'G' && rna->sequence[j] == 'C')) {
                    if (t == i)
                        rna->max_matchings[i][j] = fmax(rna->max_matchings[i][j], 1 + rna->max_matchings[t + 1][j - 1]);
                    else
                        rna->max_matchings[i][j] = fmax(rna->max_matchings[i][j], 1 + rna->max_matchings[i][t - 1] + rna->max_matchings[t + 1][j - 1]);
                }
            }
        }
    }
    rna->max_pairs = rna->max_matchings[0][rna->length - 1];
}

// Recursively find the base pairs that form the maximum number of matchings
void find_pairs(RNA_Structure* rna, int l, int r) {
    if (rna->current_pairs == rna->max_pairs || rna->max_matchings[l][r] == 0)
        return;

    if ((rna->sequence[l] == 'A' && rna->sequence[r] == 'U') ||
        (rna->sequence[l] == 'U' && rna->sequence[r] == 'A') ||
        (rna->sequence[l] == 'C' && rna->sequence[r] == 'G') ||
        (rna->sequence[l] == 'G' && rna->sequence[r] == 'C')) {
        if (rna->max_matchings[l + 1][r - 1] + 1 == rna->max_matchings[l][r]) {
            rna->pairs[rna->current_pairs].i = l;
            rna->pairs[rna->current_pairs].j = r;
            rna->current_pairs++;
            find_pairs(rna, l + 1, r - 1);
            return;
        }
    }

    for (int t = l + 1; t <= r; t++) {
        if ((rna->sequence[t] == 'A' && rna->sequence[r] == 'U') ||
            (rna->sequence[t] == 'U' && rna->sequence[r] == 'A') ||
            (rna->sequence[t] == 'C' && rna->sequence[r] == 'G') ||
            (rna->sequence[t] == 'G' && rna->sequence[r] == 'C')) {
            if (rna->max_matchings[l][t - 1] + rna->max_matchings[t + 1][r - 1] + 1 == rna->max_matchings[l][r]) {
                rna->pairs[rna->current_pairs].i = t;
                rna->pairs[rna->current_pairs].j = r;
                rna->current_pairs++;
                find_pairs(rna, l, t - 1);
                find_pairs(rna, t + 1, r - 1);
                return;
            }
        }
    }

    find_pairs(rna, l, r - 1);
}

// Set the RNA sequence for the molecule
void set_rna_sequence(RNA_Structure* rna, const char* seq) {
    rna->sequence = strdup(seq);
}

// Get the length of the RNA sequence
int get_rna_sequence_length(RNA_Structure* rna) {
    return rna->length;
}

// Get the RNA sequence
char* get_rna_sequence(RNA_Structure* rna) {
    return rna->sequence;
}

// Get the maximum number of base pairs
int get_max_rna_base_pairs(RNA_Structure* rna) {
    return rna->max_pairs;
}

// Get the execution time for calculating the maximum base pairs
double get_execution_time_rna(RNA_Structure* rna) {
    return rna->execution_time;
}

// Print the base pairs in a readable format
void print_rna_base_pairs(RNA_Structure* rna) {
    for (int i = 0; i < rna->max_pairs; i++)
        printf("%c---%c  (%d, %d)\n", rna->sequence[rna->pairs[i].i], rna->sequence[rna->pairs[i].j],
               rna->pairs[i].i, rna->pairs[i].j);
}

// Print the base pairs in a format suitable for writing to a file
void print_rna_base_pairs_to_file(RNA_Structure* rna) {
    for (int i = 0; i < rna->max_pairs; i++)
        printf("%d %d %c %c\n", rna->pairs[i].i, rna->pairs[i].j,
               rna->sequence[rna->pairs[i].i], rna->sequence[rna->pairs[i].j]);
}

// Free the memory allocated for the molecule
void free_rna(RNA_Structure* rna) {
    free(rna->sequence);
    for (int i = 0; i < rna->length; i++)
        free(rna->max_matchings[i]);
    free(rna->max_matchings);
    free(rna->pairs);
}

