#include <stdio.h>
#include <math.h>

#define NUMATOMS 3
#define BS_Size 7
#define Num_Orbitals 9
#define num_dimensions 3

//store orbital info
struct Orbital {
    double primC[num_dimensions], expC[num_dimensions], NormC[num_dimensions], center[num_dimensions];
    int angular_momentum_vector[num_dimensions], parent_atom_Z, parent_atom_idx;
} orbital_array[Num_Orbitals];

void orbital_info(struct Orbital orb, int idx);
void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer);
void get_coefs(struct Orbital orbital_array[BS_Size], FILE *coef_pointer);
void calc_norm_const(struct Orbital orbital_array[Num_Orbitals]);
double get_norm_denominator(int angular_momentum_vector[num_dimensions]);
int factorial(int n);
double overlap(double alpha, double beta, double center_a[3], double center_b[3], int ang_mom_a[3], int ang_mom_b[3]);
double little_s(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord);

int main(){
    //get info from files.
    //geometry
    FILE *geom_pointer, *coef_pointer;
    geom_pointer = fopen("water_wolfram.xyz", "r");
    coef_pointer = fopen("cont_exp_coefs.csv", "r");
    get_geom_details(orbital_array, geom_pointer);
    get_coefs(orbital_array, coef_pointer);
    calc_norm_const(orbital_array);
    // check out all orbitals in detail
    for (int i = 0; i < Num_Orbitals; i++){
        orbital_info(orbital_array[i], i);
    }

    //overlap between first primitive of orbital 0 (1s on H1) with orbital 8 (Pz on O)
    double OV = overlap(orbital_array[0].expC[0], orbital_array[8].expC[0], orbital_array[0].center, orbital_array[8].center, orbital_array[0].angular_momentum_vector, orbital_array[8].angular_momentum_vector);
    printf("Ov of the two primatives is: %lf\n", OV);

    return 0;
}

void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer){

    if (NULL == geom_pointer){
        printf("file can't be opened\n");
    }

    int orbital_idx = 0;
    int atom_idx = 0;
    int additional_orbitals;
    // Reading coords and atomic number from file
    while (!feof(geom_pointer)){
        fscanf(geom_pointer, "%d %lf %lf %lf", 
            &orbital_array[orbital_idx].parent_atom_Z, &orbital_array[orbital_idx].center[0], 
            &orbital_array[orbital_idx].center[1], &orbital_array[orbital_idx].center[2]
        );

        // printf("atom %d has atomic number %d\n",atom_idx, orbital_array[orbital_idx].parent_atom_Z);        
        if(orbital_array[orbital_idx].parent_atom_Z > 2){
            //4 more orbitals for n=2
            additional_orbitals = 4;
        } else {
            additional_orbitals = 1; //just 2S for H
        }
        orbital_array[orbital_idx].parent_atom_idx = atom_idx;
        
        //copy this information to the other orbitals on the same atom
        for(int i = 1; i <= additional_orbitals; i++){
            //copy info to orbitals on the same atom
            orbital_array[orbital_idx+i].parent_atom_idx = atom_idx;
            orbital_array[orbital_idx+i].parent_atom_Z = orbital_array[orbital_idx].parent_atom_Z;
            orbital_array[orbital_idx+i].center[0] = orbital_array[orbital_idx].center[0];
            orbital_array[orbital_idx+i].center[1] = orbital_array[orbital_idx].center[1];
            orbital_array[orbital_idx+i].center[2] = orbital_array[orbital_idx].center[2];

            if(i > 1){ //This only triggers for the p orbitals. Not worrying about d orbitals yet :)
                //correct the angular momentum vector
                orbital_array[orbital_idx+i].angular_momentum_vector[i-2] = 1;
                // Otherwise they can stay zeros for the s orbitals
            }
        }

        //now we jump to the next orbital that was not covered previously.
        orbital_idx += additional_orbitals+1;
        atom_idx++;
    }
}

void get_coefs(struct Orbital orbital_array[BS_Size], FILE *coef_pointer){
    if (NULL == coef_pointer){
        printf("file can't be opened\n");
    }

    double basis_array[BS_Size * num_dimensions][2]; //size is total num of coeffs.
    //dimension 0 is leading coeffs and dimension 1 is exponential coeffs.

    int BS_idx = 0;
    while (!feof(coef_pointer)){
        fscanf(coef_pointer, "%lf %lf %lf %lf %lf %lf", 
            &basis_array[BS_idx][0], &basis_array[BS_idx+1][0], &basis_array[BS_idx+2][0], //contraction coefs
            &basis_array[BS_idx][1], &basis_array[BS_idx+1][1], &basis_array[BS_idx+2][1] //exponential coefs
        );
        BS_idx += 3; //go to next block of empty entries.
    }
    //the basis is defined as H1s H2s O1s O2s O2p... so we can just use direct assignment.
    //The checks within the if statements deal with each atom individually since this is closer to being 
    //extendible to larger systems.
    for (int i = 0; i < Num_Orbitals; i++){
        if (orbital_array[i].parent_atom_Z == 1 && i % 2 == 0){
            //we can't yet distinguish between 1s and 2s so we should assign coefs to both now.
            //The two actually have identical coeffs, but it is still probably a good idea to assign them independently.
            for(int j = 0; j < num_dimensions; j++){
                orbital_array[i].primC[j] = basis_array[j][0]; //assign 1s orbital coefs here
                orbital_array[i].expC[j] = basis_array[j][1]; // ^ and here.
                orbital_array[i+1].primC[j] = basis_array[j + num_dimensions][0]; //now assign 2s coefs here
                orbital_array[i+1].expC[j] = basis_array[j + num_dimensions][1]; //and here too.
            }
            continue;
        }
        if (orbital_array[i].parent_atom_Z == 8 && i % 4 == 0 && Num_Orbitals-i >= 5){ //is there potentially more room for another atom?
            //be flexible so more atoms can be added later.
            for(int j = 0; j < num_dimensions; j++){
                for(int k = 0; k < 5; k++){
                    orbital_array[i+k].primC[j] = basis_array[j + num_dimensions * (2+k)][0]; //now assign O's coeffs here and next line
                    orbital_array[i+k].expC[j] = basis_array[j + num_dimensions * (2+k)][1]; //skip the H and previously read O coeffs.
                }
            }
        }
    }
}

void orbital_info(struct Orbital orb, int idx){
    printf("orbital %d has angular momentum vector: {", idx);
    for (int j = 0; j < num_dimensions; j++){
        printf(" %d", orb.angular_momentum_vector[j]);
    }
    printf(" }\n");
    char dimensions[] = "xyz";
    printf("It has the coefficients\ncontraction     exponential\n");
    for (int j = 0; j < num_dimensions; j++){
        printf("%c: %lf     %lf\n", dimensions[j], orb.primC[j],orb.expC[j]);
    }
    printf("\n");
    printf("it is located at coordinates: (");
    for (int j = 0; j < num_dimensions; j++){
        printf(" %lf", orb.center[j]);
    }
    printf(" )\n");
    printf("it has these normalization constants:");
    for (int j = 0; j < num_dimensions; j++){
        printf(" %lf", orb.NormC[j]);
    }
    printf("\n");
    printf("it is on atom #%d, which has atomic number %d\n\n----------------\n", orb.parent_atom_idx, orb.parent_atom_Z);
}

int factorial(int n){
    int fact = 1;
    for (int i = n; i >= -1; i--){
        if(i == 0 || i == -1){
            fact *= 1;
        } else {
            fact *= i;
        }
    }
    return fact;
}

double get_norm_denominator(int angular_momentum_vector[num_dimensions]){
    double den;
    int part, fact[num_dimensions];
    //each component of the ang. mom. vector is used.
    for(int i = 0; i < num_dimensions; i++){
        part = 2 * angular_momentum_vector[i] - 1;
        fact[i] = factorial(part);
        fact[i] = factorial(fact[i]);
    }
    //the denominator is the square root of the product of all double factorials
    den = sqrt(fact[0] * fact[1] * fact[2]);
    return den;
}

void calc_norm_const(struct Orbital orbital_array[Num_Orbitals]){
    double prefactor, angular_part, alpha, denominator;
    int sum_angular_coords;
    //for each orbital
    for (int i = 0; i < Num_Orbitals; i++){
        sum_angular_coords = 0;
        //calculate the sum of the angular momentum vector as this is used in each calculation of the normalization constant.
        for (int j = 0; j < num_dimensions; j++){
            sum_angular_coords += orbital_array[i].angular_momentum_vector[j];
        }

        //calculate the denominator. This is the same regardless of alpha.
        denominator = get_norm_denominator(orbital_array[i].angular_momentum_vector);

        // finally, calculate the normalization constant for each dimension.
        for (int j = 0; j < num_dimensions; j++){
            alpha = orbital_array[i].expC[j];
            angular_part = pow((4 * alpha), (sum_angular_coords * 0.5)) / denominator;
            prefactor = pow((2 * alpha / M_PI), 0.75 );
            orbital_array[i].NormC[j] = prefactor * angular_part;
        }
    }     
}

double little_s(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord){
    double P, sum_ab;

    // printf("a: %d, b: %d\n", ang_coord_a, ang_coord_b);
    
    if (ang_coord_a == 0 && ang_coord_b == 0){
        // printf("returning 1\n");
        return 1;
    }

    sum_ab = alpha + beta;
    P = (alpha*center_a_coord + beta*center_b_coord) / sum_ab;

    if (ang_coord_a == 1 && ang_coord_b == 0){
        // printf("Simple Case\n");
        return -(center_a_coord - P);
    }
    if ( ang_coord_a != 1 && ang_coord_b == 0){ //recurrence index
        double s_prev_a, s_prev_2_a;
        s_prev_a = little_s(ang_coord_a - 1, 0, alpha, beta, center_a_coord, center_b_coord);
        s_prev_2_a = little_s(ang_coord_a - 2, 0, alpha, beta, center_a_coord, center_b_coord);
        // printf("Using intermediate values that use lower values of a: %lf and %lf\n", s_prev_a, s_prev_2_a);
        return -(center_a_coord - P) * s_prev_a + ((ang_coord_a - 1) / (2*sum_ab)) * s_prev_2_a;
    }
    if ( ang_coord_b != 0 ){ //transfer
        double s_xfer, s_prev_b;
        s_xfer = little_s(ang_coord_a + 1, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord);
        s_prev_b = little_s(ang_coord_a, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord);
        // printf("Using intermediate values that use lower values of a: %lf and %lf\n", s_xfer, s_prev_b);
        return s_xfer + (center_a_coord - center_b_coord) * s_prev_b;
    }
}

double overlap(double alpha, double beta, double center_a[3], double center_b[3], int ang_mom_a[3], int ang_mom_b[3]){
    double EAB, exponent, dist_squared, Overlap;
    dist_squared = 0;   
    for (int i = 0; i < num_dimensions; i++){
        dist_squared += pow(center_a[i]-center_b[i], 2);
    }
    exponent = -(alpha * beta / (alpha + beta)) * dist_squared;
    EAB = pow(M_E, exponent);

    Overlap = EAB * pow((M_PI / (alpha + beta)), 1.5);
    // printf("Overlap Coeff is %lf\n", Overlap);

    for (int i = 0; i < num_dimensions; i++){
        Overlap *= little_s(ang_mom_a[i], ang_mom_b[i], alpha, beta, center_a[i], center_b[i]);
        // printf("After calculating s(%d), the Overlap Integral is %lf\n", i, Overlap);
    }

    return Overlap;

}

