/*RA*/
//need
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <sys/time.h>
#include <math.h>

//Assumes input clustered by Genome(Not necessarily sorted)

//need to reduce memory accesses and str function calls by creating extra variables
//need to allocate exact memory according to number of genomes and number of lines in file
//need to check <num_genomes or num_genomes+1
//need to reduce stored[] to max_protein
//need to eliminate 1 iteration loops and extra ifs for when V = 1

void computeSubResults(int *results, int level, int *genome_indices, int num_genomes, int max_protein, int *proteins, double threshold, short *seen, int *old_figfams, char **functions, int *real_proteins);
int* Prob_A_Bs(int *FFV, int Vlen, int *booleans, int *genome_indices, int num_genomes, int max_protein, int *proteins, double threshold, short *seen, int *old_figfams, char **functions, int *real_proteins);

int main(int argc, char* argv[])
{
    if (argc < 6) {
        printf("Usage: ./C1235abUnique <filename> <n=len(V)> <V={v1 ... vn}> <threshold> <level>\n");
        exit(0);
    }
    
    FILE *input;
    input = fopen(argv[1], "r");
    
    char **functions = (char**)malloc(1958965 * sizeof(char*));
    
    int Vlen = atoi(argv[2]);
    int *FFV = (int*)malloc(Vlen * sizeof(int));

    int booleans[Vlen]; // To store the state of the variables
    
    int i = 0;
    for (i = 0; i < Vlen; i++) {
        char c = argv[3+i][0];
        booleans[i] = c == '^' ? 0 : 1;
    }
    
    double threshold = atof(argv[Vlen + 3]);

    int level = atoi(argv[Vlen+4]);

    
    char GenomeID[20];
    char FFamID[256];
    char current_pr[256];

    //Clean string format to keep the number after FIG
    int j = 0;
    int k = 0;
    for (i = 0; i < Vlen; i++) {
        k = booleans[i] ? 3 : 4; // skip FIG or ^FIG
        for(j = k; j < strlen(argv[3+i]); j++) {
            current_pr[j-k] = argv[3+i][j];
        }
        FFV[i] = atoi(current_pr);
    }
    
    long long start_sum = 0;
    int num_start = 0;
    
    int line_count = 0;
    
    int max_protein = 0;
    char current_protein[256];
    int current_index = 0;
    int m = 0;
    int num_genomes = 0;
    int *genome_indices = (int*)malloc(27000 * sizeof(int)); //start,end of proteins for each genome
    int *proteins = (int*)malloc(27445281 * sizeof(int));
    int *old_figfams = (int*)malloc(1958965 * sizeof(int));
    unsigned short *stored = (unsigned short*)malloc(1958965 * sizeof(unsigned short)); //to record real proteins
    int *real_proteins = (int*)malloc(4906 * sizeof(int));

    for (i = 0; i < 27445281; i++) {
        if (i < 27000) {
            genome_indices[i] = 0;
        }
        if (i < 1958965) {
            old_figfams[i] = 0;
            functions[i] = (char*) malloc(777 * sizeof(char));
            stored[i] = 0;
        }
        proteins[i] = 0;
    }
    
    char function[777];
    char previousG[20];
    i = 0;
    int first = 1;
    int startd, endd;
    char genus[30];
    char *s = NULL;
    size_t linesize = 0;
    ssize_t linelen;
    char old_fig[256];
    char current_old_fig[256];
    int old_index = 0;
    int num_proteins = 0;
    while ((linelen = getline(&s, &linesize, input)) != -1) {
        line_count++;
     //   printf("%d\n", line_count);
        m = sscanf(s, "%s\t%s\t%d\t%d\t%[^\t]%s\t%s\n", GenomeID, genus, &startd, &endd, function, FFamID, old_fig);
        if (first) {
            strcpy(previousG, GenomeID);
            first = 0;
            num_genomes++;
        }
        if (strcmp(previousG, GenomeID) != 0) { //We found a new Genome
            genome_indices[num_genomes-1] = i;
            strcpy(previousG, GenomeID);
            num_genomes++;
        }
        if (m == 7) {
            for(k = 3; k < strlen(FFamID); k++) {
                current_protein[k-3] = FFamID[k];
            }
            current_index = atoi(current_protein);
            if (FFamID[0] == 'F' && FFamID[1] == 'I' && FFamID[2] == 'G') {
                if (current_index > max_protein) {
                    max_protein = current_index;
                }
                
                strcpy(functions[current_index], function); //store current FIGFAM function
                if (current_index == 17 && strcmp(genus, "Escherichia")) {
                    if (startd < 0) {
                        printf("Start: %d\n", startd);
                    }
                    start_sum += startd;
                    num_start++;
                }
                if (stored[current_index] == 0) {
                    real_proteins[num_proteins] = current_index;
                    num_proteins++;
                }
                for(k = 3; k < strlen(old_fig); k++) {
                    current_old_fig[k-3] = old_fig[k];
                }
                old_index = atoi(current_old_fig);
                if (stored[current_index] < num_genomes) { //current_index
                    proteins[i] = current_index;//current_index
                    stored[current_index] = num_genomes;//current_index
                    i++;
                }
                
                old_figfams[current_index] = old_index; //save old figfam
            }
        }
    }
    fclose(input);
    printf("Number of lines read = %d\n",line_count);
    printf("Number of Genomes = %d\n", num_genomes);
    printf("Max index: %d\n", max_protein);
    printf("Number of Proteins: %d\n", num_proteins);
    
    printf("total start: %lld\ttotal num: %d\n", start_sum, num_start);
    printf("Start average: %f\n", ((float)start_sum)/num_start);
    
    max_protein = max_protein + 1; //To be able to store at max_protein index
    
    struct timeval tt;
    long tsec,tusec;
    float e;
    gettimeofday(&tt,0);
    tsec = tt.tv_sec;
    tusec = tt.tv_usec;
    
    int *results;
    short *seen = (short*)malloc((max_protein+1) * sizeof(short));

    memset(seen, 0, max_protein);
    results = Prob_A_Bs(FFV, Vlen, booleans, genome_indices, num_genomes, max_protein, proteins, threshold, seen, old_figfams, functions, real_proteins);
    char fig[30];
/*    FILE *figfams = fopen("hypothetical-figfams.txt", "r");
    while (m = fscanf(figfams, "%s\n", fig) > 0) {
        if (feof(figfams)) {
            break;
        }
        for(j = 3; j < strlen(fig); j++) {
            current_pr[j-3] = fig[j];
        }
        FFV[0] = atoi(current_pr);
        results = Prob_A_Bs(FFV, Vlen, booleans, genome_indices, num_genomes, max_protein, proteins, threshold, seen, old_figfams, functions, real_proteins);
    }
    fclose(figfams);
*/
    if (level > 1) {
        printf("End of Level 1.\n");
        computeSubResults(results, --level, genome_indices, num_genomes, max_protein, proteins, threshold, seen, old_figfams, functions, real_proteins);
    }
    
    int unique = 0;
 /*   for (i = 0; i < max_protein; i++) {
        if (seen[i]) {
            unique++;
        }
    }
*/ 
    printf("Unique: %u\n", unique);
    printf("Repeats: %u\n", seen[max_protein]);

    gettimeofday(&tt,0);
    e = (float)(tt.tv_sec - tsec) + (float)(tt.tv_usec - tusec)/1000000;
    printf( "Elapsed time seq. :  %le ms\n", 1000*e);

    free(genome_indices);
    free(proteins);
    free(stored);
    free(results);
    free(seen);
    
    printf ("Proper Termination \n");

    exit(0);
}

void computeSubResults(int *results, int level, int *genome_indices, int num_genomes, int max_protein, int *proteins, double threshold, short *seen, int *old_figfams, char **functions, int *real_proteins) {
    if (results[0] < 1) {
        printf("Program converged.\n");
        exit(0);
    }
    int **new_results = (int**)malloc(results[0] * sizeof(int*));
    int i;
    int booleans[1];
    for (i = 1; i < results[0]+1; i++) {
        //            printf("V -> FIG%08d\n", results[i]);
        booleans[0] = 1;

        new_results[i] = Prob_A_Bs(&results[i], 1, booleans, genome_indices, num_genomes, max_protein, proteins, threshold, seen, old_figfams, functions, real_proteins);
        if (level > 1) {
            computeSubResults(new_results[i], level-1, genome_indices, num_genomes, max_protein, proteins, threshold, seen, old_figfams, functions, real_proteins);
        }
    }
    printf("End of Level %d.\n", level+1);
}

int* Prob_A_Bs(int *FFV, int Vlen, int *booleans, int *genome_indices, int num_genomes, int max_protein, int *proteins, double threshold, short *seen, int *old_figfams, char **functions, int *real_proteins) { //A, genome_indices..., Bs, threshold
    
    unsigned short countFFAs[Vlen]; //frequency of individual FIGFAMs in V

    //frequency of co-occurences of V with B
    unsigned short *countFFVBs = (unsigned short*)malloc((max_protein) * sizeof(unsigned short));
    int *countBs = (int*)malloc((max_protein) * sizeof(int));
    unsigned short *countFFVnBs = (unsigned short*)malloc((max_protein) * sizeof(unsigned short));
    unsigned short *countnFFVBs = (unsigned short*)malloc((max_protein) * sizeof(unsigned short));
    unsigned short *countnFFVnBs = (unsigned short*)malloc((max_protein) * sizeof(unsigned short));

    int i = 0;

    memset(countFFVBs, 0, max_protein);
    memset(countnFFVBs, 0, max_protein);
    memset(countFFVnBs, 0, max_protein);
    memset(countnFFVnBs, 0, max_protein);
    memset(countBs, 0, max_protein);


    for (i = 0; i < Vlen; i++) {
        countFFAs[i] = 0;
    }
    int countFFV = 0; //frequency of full vector V
    
    int genome_score = 0; //records members of V in current genome
    
    //records if an element of V has been seen when not supposed to
    int oops = 0;
    
    int startg = 0, endg = 0;
    int current_p;
    int n = 0, j = 0, k = 0;
    
    for (i = 0; i < num_genomes; i++) {
        genome_score = 0;
        oops = 0;
        endg = genome_indices[i];
        
        for (j = startg; j < endg; j++) {
            current_p = proteins[j];
            if (current_p >=0) {
                countBs[current_p]++; //count frequency of each protein
                for (k = 0; k < Vlen; k++) {
                    if (current_p == FFV[k]) {
                        if (booleans[k]) {
                            genome_score++;
                            countFFAs[k]++;
                        }
                        else {
                            oops = 1;
                        }
                    }
                }
                if (oops || (genome_score == Vlen)) { // V is not in this Genome or V is in this genome
               //     break;
                }
            }
        }
        if (!oops) {
            for (n = 0; n < Vlen; n++) {
                if (!booleans[n]) {
                    genome_score++; // element was rightfully absent from genome
                    countFFAs[n]++;
                }
            }
            if (genome_score == Vlen) { // V is found, look for Bs
                countFFV++;
                for (n = startg; n < endg; n++) {
                    if (proteins[n] >= 0) {
                        countFFVBs[proteins[n]]++;// = real[proteins[n]] ? countFFVBs[proteins[n]]+1 : 0; //only increment real ones
                    }
                }
            }
            /*else {
                for (n = startg; n < endg; n++) {
                    if (proteins[n] >= 0) {
                        countnFFVBs[proteins[n]]++;
                    }
                }
            }*/
        }
        
        startg = endg;
    }
    
    
    
    double *probabilities_A = (double*)malloc(Vlen * sizeof(double)); //individual FIGFAMs in V
    double *probabilities_Bs = (double*)malloc(max_protein * sizeof(double));
    
    for (i = 0; i < Vlen; i++) {
        probabilities_A[i] = ((double)countFFAs[i])/num_genomes;
    }
    double probability_V = ((double)countFFV)/num_genomes;
/*
    //Fix string format for printing by adding ^ if it was removed earlier
    for (i = 0; i < Vlen; i++) {
        if (booleans[i]) {
            printf("P(FIG%08d) = %lf\n",FFV[i], probabilities_A[i]);
        }
        else {
            printf("P(^FIG%08d) = %lf\n",FFV[i], probabilities_A[i]);
            
        }
    }
    
    printf("V =");
    for (i = 0; i < Vlen; i++) {
        if (booleans[i]) {
            printf("\tFIG%08d", FFV[i]);
        }
        else {
            printf("\t^FIG%08d", FFV[i]);
        }
    }
    printf("\n");
    
    printf("Frequency of V = %d\n", countFFV);
    printf("P(V) = %lf\n", probability_V);
*/
    double maxMI = 0;
    double max7[10];
    int max7i[10];
    for (i = 0; i < 10; i++) {
        max7[i] = 0;
        max7i[i] = 0;
    }
    
    double probability_VB = 0, probability_nVB = 0, probability_VnB = 0, probability_nVnB = 0;
    double probability_B_given_V = 0;
    double MI = 0, MI00 = 0, MI01 = 0, MI10 = 0, MI11 = 0;
    double probability_B = 0;
    int count = 0;
    int countMI = 0;

    for (j = 0; j < 4906; j++) {
        i = real_proteins[j];
        if(countBs[i] > 0){
            probability_VB = ((double)countFFVBs[i])/num_genomes;
            countnFFVBs[i] = countBs[i] - countFFVBs[i];
            probability_nVB = ((double)countnFFVBs[i])/num_genomes;
            countFFVnBs[i] = countFFV - countFFVBs[i];
            probability_VnB = ((double)countFFVnBs[i])/num_genomes;
            countnFFVnBs[i] = num_genomes - (countFFVBs[i] + countnFFVBs[i] + countFFVnBs[i]);
            probability_nVnB = ((double)countnFFVnBs[i])/num_genomes;
            
            if (probability_VB + probability_nVnB + probability_VnB + probability_nVB !=1) {
                //          printf("Total P = %lf\n", probability_VB + probability_nVnB + probability_VnB + probability_nVB);
                
            }
            if (countFFVBs[i] + countnFFVBs[i] + countnFFVnBs[i] + countFFVnBs[i] != num_genomes && countFFVBs[i] + countnFFVBs[i] + countnFFVnBs[i] + countFFVnBs[i] != 0) {
                //                printf("Total C = %lf\n", countFFVBs[i] + countnFFVBs[i] + countnFFVnBs[i] + countFFVnBs[i]);
                
            }
            
            probability_B_given_V = probability_VB/probability_V;
            probability_B = ((double)countBs[i])/num_genomes;
            
            if (probability_VB > 0 && probability_V > 0 && probability_B > 0) {
                double px = probability_VB + probability_VnB;
                double py = probability_nVB + probability_VB;
                double denomenator = px * py;
                MI11 = probability_VB * (log(probability_VB/denomenator)/log(2.0));
            }
            else {
                MI11 = 0;
            }
            if (probability_nVB > 0 && probability_V < .9999999 && probability_B > 0) {
                double px = probability_nVB + probability_nVnB;
                double py = probability_nVB + probability_VB;
                double denomenator = px * py;
                MI01 = probability_nVB * (log(probability_nVB/denomenator)/log(2.0));
            }
            else {
                MI01 = 0;
            }
            if (probability_VnB > 0 && probability_V > 0 && probability_B < .9999999) {
                double px = probability_VnB + probability_VB;
                double py = probability_nVnB + probability_VnB;
                double denomenator = px * py;
                MI10 = probability_VnB * (log(probability_VnB/denomenator)/log(2.0));
            }
            else {
                MI10 = 0;
            }
            if (probability_nVnB > 0 && probability_V < .9999999 && probability_B < .9999999) {
                double px = probability_nVnB + probability_nVB;
                double py = probability_nVnB + probability_VnB;
                double denomenator = px * py;
                MI00 = probability_nVnB * (log(probability_nVnB/denomenator)/log(2.0));
            }
            else {
                MI00 = 0;
            }
            MI = MI00 + MI01 + MI10 + MI11;
            
            if (probability_B_given_V >=threshold) {
                count++;
            }
            if (MI > 0.01) {
              //  fprintf(hypo, "FIG%08d\tFIG%08d\t%s\t%lf\n",FFV[0], old_figfams[i], functions[i], MI);
            }
            if (MI > max7[0] && i != FFV[0]) {
                max7[0] = MI;
                max7i[0] = i;
                
                for (n = 1; n < 10; n++) {
                    if (max7[n] < max7[0]) {
                        double t = max7[0];
                        max7[0] = max7[n];
                        max7[n] = t;
                        int ti = max7i[0];
                        max7i[0] = max7i[n];
                        max7i[n] = ti;
                    }
                }
            }

        }
    }
    //printf("FIG%08d: %s\n", old_figfams[FFV[0]], functions[FFV[0]]);
    count = 10;
    int *results = (int*)malloc((count+1) * sizeof(int));
    results[0] = count;
    
    for (n = 0; n < 10; n++) {
        printf( "FIG%08d\tFIG%08d\t%lf\n",old_figfams[FFV[0]], old_figfams[max7i[n]], max7[n]);
        results[n+1] = max7i[n];
    }
 //   printf("Count = %d\n", count);
 //   printf("CountMI = %d\n", countMI);

/*
    k = 1;
    for (i = 0; i < max_protein; i++) {
        probability_VB = ((double)countFFVBs[i])/num_genomes;
        probability_B_given_V = probability_VB/probability_V;
        if (probability_B_given_V >= threshold && i != FFV[0]) { //don't store my own(reduces recursion)
 //                      printf("V -> FIG%08d with probability = %lf\n", i, probability_B_given_V);
            if (!seen[i]) {
                results[k] = i;
                k++;
                seen[i] = 1;
            }
            else {
                seen[max_protein]++;
            }
        }
    }
*/
    //Free allocated memory
    free(countFFVBs);
    free(countFFVnBs);
    free(countnFFVBs);
    free(countnFFVnBs);
    free(countBs);
    free(probabilities_A);
    free(probabilities_Bs);
    
    return results;
}
