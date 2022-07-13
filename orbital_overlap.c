#include <stdio.h>
#include <math.h>
#include <stdlib.h> /*for system()*/

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
double primitive_overlap(int dim_a, int dim_b, struct Orbital orbital_a, struct Orbital orbital_b);
double little_s(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord);
double orbital_overlap(struct Orbital orbital_a, struct Orbital orbital_b);
void Calc_BS_OV_Matrix(struct Orbital orbital_array[Num_Orbitals], double overlap_matrix[BS_Size][BS_Size], int indicies[BS_Size]);
double little_k(int dimension, struct Orbital orbital_a, struct Orbital orbital_b, int exp_idx);
double kinetic_energy_integral(int dimension, struct Orbital orbital_a, struct Orbital orbital_b);
double dist_squared(struct Orbital orbital_a, struct Orbital orbital_b);

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
        // orbital_info(orbital_array[i], i);
    }

    //overlap between first primitive of orbital 0 (1s_x on H1) with orbital 8 (Pz_x on O)
    // double OV = primitive_overlap(0, 0, orbital_array[0], orbital_array[8]]);
    // printf("Ov of the two primatives is: %lf\n", OV);

    //overlap between orbital 0 and 8 (H1_1s and O_p_z)
    // double INTEGRAL = orbital_overlap(orbital_array[0], orbital_array[8]);
    // printf("Overlap integral is %lf\n",INTEGRAL);

    //overlap between orbital 2 and 3 (H1_1s and H1_2s)
    // double INTEGRAL = orbital_overlap(orbital_array[0], orbital_array[7]);
    // printf("Overlap integral is %lf\n", INTEGRAL);

    double BS_overlap_matrix[BS_Size][BS_Size];
    int included_indicies[BS_Size] = {0,3,4,5,6,7,8};
    Calc_BS_OV_Matrix(orbital_array, BS_overlap_matrix, included_indicies);

    system("clear"); /*clear output screen*/
    
    // KE integral for first primatives of orbital 0 (1s_x on H1) with orbital 8 (Pz_x on O) 
    double KE = kinetic_energy_integral(0, orbital_array[0], orbital_array[8]);
    printf("KE integral is %lf\n", KE);  

    // printf("problematic little s:\n");
    // double problem = little_s(orbital_array[0].angular_momentum_vector[2], orbital_array[8].angular_momentum_vector[2], orbital_array[0].expC[0], orbital_array[8].expC[0], orbital_array[0].center[2], orbital_array[8].center[2]);
    // printf(" %lf\n", problem);
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
    printf("it is on atom #%d, which has atomic number %d\n----------------\n\n", orb.parent_atom_idx, orb.parent_atom_Z);
}

double dist_squared(struct Orbital orbital_a, struct Orbital orbital_b){
    double dist_squared = 0;   
    for (int i = 0; i < num_dimensions; i++){
        dist_squared += pow(orbital_a.center[i] - orbital_b.center[i], 2);
    }
    return dist_squared;
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

    // printf("%d, %d, %lf, %lf\n", ang_coord_a, ang_coord_b, center_a_coord, P);
    if (ang_coord_a == 0 && ang_coord_b == 0){ //definition part 1
        printf("    s(0,0)=1\n");
        return 1;
    }

    sum_ab = alpha + beta;
    P = (alpha*center_a_coord + beta*center_b_coord) / sum_ab;
    // printf("alpha %lf, beta %lf, center a %lf, center b %lf\n", alpha, beta, center_a_coord, center_b_coord);
    // printf("P is %lf\n", P);

    if(ang_coord_a == 1 && ang_coord_b == 0){ //definition part 2
        printf("    s(1,0)=%lf\n", -(center_a_coord - P));
        // printf("Simple Case: -(%lf-%lf) = %lf\n",center_a_coord,P,-(center_a_coord - P));
        return -(center_a_coord - P);
    }
    if(ang_coord_a != 1 && ang_coord_b == 0){ //recurrence index
        double s_prev_a, s_prev_2_a;
        printf("s(%d,%d) requires recurrance\n",ang_coord_a,ang_coord_b);
        s_prev_a = little_s(ang_coord_a - 1, 0, alpha, beta, center_a_coord, center_b_coord);
        s_prev_2_a = little_s(ang_coord_a - 2, 0, alpha, beta, center_a_coord, center_b_coord);
        // printf("Using intermediate values that use lower values of a: %lf and %lf\n", s_prev_a, s_prev_2_a);
        double value = -(center_a_coord - P) * s_prev_a + ((ang_coord_a - 1) / (2*sum_ab)) * s_prev_2_a;
        printf("s(%d,%d) is %lf\n",ang_coord_a,ang_coord_b,value);
        return value;
    }
    if(ang_coord_b != 0){ //transfer
        double s_xfer, s_prev_b;
        printf("s(%d,%d) requires transfer\n",ang_coord_a,ang_coord_b);
        s_xfer = little_s(ang_coord_a + 1, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord);
        s_prev_b = little_s(ang_coord_a, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord);
        // printf("Using intermediate values that use different values of a and b: %lf and %lf\n", s_xfer, s_prev_b);
        double value = s_xfer + (center_a_coord - center_b_coord) * s_prev_b; 
        printf("s(%d,%d) is %lf\n",ang_coord_a,ang_coord_b,value);
        return value;
    }
}

double primitive_overlap(int dim_a, int dim_b, struct Orbital orbital_a, struct Orbital orbital_b){
    double alpha = orbital_a.expC[dim_a];
    double beta = orbital_b.expC[dim_b];
    double EAB, exponent, dist_squrd, Overlap;
    
    dist_squrd = dist_squared(orbital_a, orbital_b);
    exponent = -(alpha * beta / (alpha + beta)) * dist_squrd;
    EAB = pow(M_E, exponent);

    Overlap = EAB * pow((M_PI / (alpha + beta)), 1.5);
    // printf("Overlap Coeff is %lf\n", Overlap);

    for (int i = 0; i < num_dimensions; i++){
        Overlap *= little_s(orbital_a.angular_momentum_vector[i], orbital_b.angular_momentum_vector[i], alpha, beta, orbital_a.center[i], orbital_b.center[i]);
        // printf("After calculating s(%d), the Overlap Integral is %lf\n", i, Overlap);
    }

    return Overlap;
}

double orbital_overlap(struct Orbital orbital_a, struct Orbital orbital_b){
    double integral = 0;
    for(int i = 0; i < num_dimensions; i++){
        for(int j = 0; j < num_dimensions; j++){
            double ov = primitive_overlap(i,j,orbital_a,orbital_b);
            integral += orbital_a.NormC[i] * orbital_b.NormC[j] * orbital_a.primC[i] * orbital_b.primC[j] * ov;
            // printf("ov is %lf. Integral is now %lf after using %lf %lf %lf %lf\n", ov, integral, orbital_a.NormC[i], orbital_b.NormC[j], orbital_a.primC[i], orbital_b.primC[j]);
        }
    }
    return integral;
}

void Calc_BS_OV_Matrix(struct Orbital orbital_array[Num_Orbitals], double overlap_matrix[BS_Size][BS_Size], int indicies[BS_Size]){
    //Since there are 9 orbitals in the sytem and 7 in the BS, exclude some orbitals to avoid repeats.
    struct Orbital used_orbitals[BS_Size];
    for (int i = 0; i < BS_Size; i++){
        int next_idx = indicies[i];
        used_orbitals[i] = orbital_array[next_idx];
    }
    printf("overlap matrix:\n");
    for(int i = 0; i < BS_Size; i++){
        for(int j = 0; j < BS_Size; j++){
            overlap_matrix[i][j] = orbital_overlap(used_orbitals[i], used_orbitals[j]);
            printf("    %lf", overlap_matrix[i][j]);
        }
        printf("\n");
    }
}

double little_k(int dimension, struct Orbital orbital_a, struct Orbital orbital_b, int exp_idx){
    double alpha = orbital_a.expC[exp_idx];
    double beta = orbital_b.expC[exp_idx];
    int a = orbital_a.angular_momentum_vector[dimension];
    int b = orbital_b.angular_momentum_vector[dimension];
    double coord_a = orbital_a.center[dimension];
    double coord_b = orbital_b.center[dimension];
    double al_bet = alpha * beta;

    printf("k_%d_(%d,%d)\n",dimension,a,b);
    //initial conditions for KE:
    if(a == 0 && b == 0){
        return 2 * al_bet * little_s(1, 1, alpha, beta, coord_a, coord_b);
    }
    if(a > 0 && b == 0){
        return -a * beta * little_s(a-1, 1, alpha, beta, coord_a, coord_b) + 2 * al_bet * little_s(a+1, 1, alpha, beta, coord_a, coord_b);
    }
    if(a == 0 && b > 0){
        return -alpha * b * little_s(1, b-1, alpha, beta, coord_a, coord_b) + 2 * al_bet * little_s(1, b+1, alpha, beta, coord_a, coord_b);
    }
    //kinetic energy integral
    if(a > 0 && b > 0){
        double s_lower = little_s(a-1, b-1, alpha, beta, coord_a, coord_b);
        double s_down_a = little_s(a-1, b+1, alpha, beta, coord_a, coord_b);
        double s_up_a = little_s(a+1, b-1, alpha, beta, coord_a, coord_b);
        double s_upper = little_s(a+1, b+1, alpha, beta, coord_a, coord_b);

        return (a * b * s_lower - 2 * a * beta * s_down_a - 2 * alpha * b * s_up_a + 4 * al_bet * s_upper) / 2;
    }
}

double kinetic_energy_integral(int dimension, struct Orbital orbital_a, struct Orbital orbital_b){
    //KE integral of 2 gaussian primitives i.e. Chi_ij and K_x
    double alpha = orbital_a.expC[dimension];
    double beta = orbital_b.expC[dimension];
    double coord_a = orbital_a.center[dimension];
    double coord_b = orbital_b.center[dimension];
    double sumAB = alpha + beta;

    double dist_squrd = dist_squared(orbital_a, orbital_b);

    double EAB = pow(M_E, -(alpha*beta / sumAB) * dist_squrd);

    double coeff = pow((M_PI / sumAB), 1.5);

    double little_integrals[3][2]; //3 little s and 3 little k
    for(int i = 0; i < num_dimensions; i++){
        printf("\ni: %d\n",i);
        little_integrals[i][0] = little_s(orbital_a.angular_momentum_vector[i], orbital_b.angular_momentum_vector[i], alpha, beta, orbital_a.center[i], orbital_b.center[i]);
        little_integrals[i][1] = little_k(i, orbital_a, orbital_b, dimension);

        printf("Little k: %lf. Little s: %lf\n", i, little_integrals[i][1], little_integrals[i][0]);
    }

    double KE = EAB * coeff;
    double product;

    for(int i = 0; i < num_dimensions; i++){
        product = 1;
        // printf("%d, %d, %d\n",i,(i+1)%3, (i+2)%3);
        product *= little_integrals[i][1]; // the little k with index i
        product *= little_integrals[(i+1) % 3][0] * little_integrals[(i+2) % 3][0]; // The little s with the other indices
        printf("\nlittle integrals: k %d     s %d     s %d\n",i,(i+1)%3,(i+2)%3);
        printf("            %lf      %lf     %lf\n", little_integrals[i][1], little_integrals[(i+1) % 3][0], little_integrals[(i+2) % 3][0]);
        printf("product %lf\n", product);
        KE += product;
    }
    
    return KE;
}