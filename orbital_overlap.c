#include <stdio.h>
#include <math.h>
#include <stdlib.h> /*for system()*/

#define NUMATOMS 3
#define BS_Size 7
#define Num_Orbitals 9
#define num_dimensions 3
#define EPS 10E-10 
#define MAX_POLYNOMIAL_SIZE 10

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
int fact2(int n);
double dist_squared(double coords_A[num_dimensions], double coords_B[num_dimensions]);
double dot_product(double coords_A[num_dimensions], double coords_B[num_dimensions]);
void scalar_mult(double *coords_pt, double scalar, int size);
double primitive_overlap(int dim_a, int dim_b, struct Orbital orbital_a, struct Orbital orbital_b);
double little_s(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord);
double orbital_overlap(struct Orbital orbital_a, struct Orbital orbital_b);
void Calc_BS_OV_Matrix(struct Orbital orbital_array[Num_Orbitals], double overlap_matrix[BS_Size][BS_Size], int indicies[BS_Size]);
double little_k(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord);
double orbital_kinetic_energy_integral(struct Orbital orbital_a, struct Orbital orbital_b);
double primitives_KE(int primitive_idx_a, int primitive_idx_b, struct Orbital orbital_a, struct Orbital orbital_b);
void Calc_BS_KE_Matrix(struct Orbital orbital_array[Num_Orbitals], double KE_matrix[BS_Size][BS_Size], int indicies[BS_Size]);
double abscissa(int n, int i);
double boys_func(double x, int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]);
double omega(int n, int i);
double little_n(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord, double t, double nuc_coord);
void alt_little_n(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord, double nuc_coord, double *polynomial_pointer);
double chebychev_integral_boys(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]);
double N_e_attraction(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]);
double hyp1f1_clone(double a, double b, double x);
double hyp1f1_int_boys(double polynomial_terms[2], double alpha, double beta, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]);
double boys_func_hyp(double n, double T);
void foil_polynomials(double *polynomial_ptr_1, double *polynomial_ptr_2, double *result_ptr);

int main(){
    //get info from files.
    FILE *geom_pointer, *coef_pointer;
    geom_pointer = fopen("water_wolfram.xyz", "r");
    coef_pointer = fopen("cont_exp_coefs.csv", "r");
    //geometry
    get_geom_details(orbital_array, geom_pointer);
    //coeffs
    get_coefs(orbital_array, coef_pointer);
    calc_norm_const(orbital_array);
    // check out all orbitals in detail
    // for (int i = 0; i < Num_Orbitals; i++){
    //     orbital_info(orbital_array[i], i);
    // }

    // overlap between first primitive of orbital 0 (1s_x on H1) with orbital 8 (Pz_x on O)
    // double OV = primitive_overlap(0, 0, orbital_array[0], orbital_array[8]);
    // printf("Ov of the two primatives is: %lf\n", OV);

    // overlap between orbital 0 and 8 (H1_1s and O_p_z)
    // double INTEGRAL = orbital_overlap(orbital_array[0], orbital_array[8]);
    // printf("Overlap integral is %lf\n", INTEGRAL);

    // // overlap between orbital 2 and 3 (H1_1s and H1_2s)
    // double INTEGRAL = orbital_overlap(orbital_array[0], orbital_array[7]);
    // printf("Overlap integral is %lf\n", INTEGRAL);

    double BS_overlap_matrix[BS_Size][BS_Size];
    int included_indicies[BS_Size] = {0,3,4,5,6,7,8};
    // Calc_BS_OV_Matrix(orbital_array, BS_overlap_matrix, included_indicies);

    // system("clear"); /*clear output screen*/
    struct Orbital orbital_a = orbital_array[0];
    struct Orbital orbital_b = orbital_array[8];
    
    // KE integral for first primatives of orbital 0 (1s_x on H1) with orbital 8 (Pz_x on O) 
    // KE = primitives_KE(0, orbital_array[0], orbital_array[8]);
    // printf("\nKE integral is %lf\n", KE); //slightly off. 0.001887 instead of 0.00167343. Maybe due to different e or pi values. EAB is right for the digets shown, but there are more in the paper that aren't here.

    // K_08
    // KE = orbital_kinetic_energy_integral(orbital_a, orbital_b);
    // printf("KE integral is %lf\n", KE);
    // Int of prims 0 0 has KE of // K_17
    // double KE = orbital_kinetic_energy_integral(orbital_a, orbital_b);
    // printf("KE integral is %lf\n", KE);
    // Int of prims 0 0 has KE of -0.167203 is output. Correct

    //KE overlap integral
    // Calc_BS_KE_Matrix(orbital_array, BS_overlap_matrix, included_indicies);

    //testing NE functions
    double results;
    // results = alt_little_n(0, 0, orbital_a.expC[0], orbital_b.expC[0], orbital_a.center[0], orbital_b.center[0], 0, orbital_a.center[0]);
    // results = exp(-(orbital_a.expC[0] * orbital_b.expC[0])/(orbital_a.expC[0] + orbital_b.expC[0]) * dist_squared(orbital_a.center, orbital_b.center))
    //     * (2 * M_PI / (orbital_a.expC[0] + orbital_b.expC[0]));
    // results = little_n(0, 0, orbital_a.expC[0], orbital_b.expC[0], orbital_a.center[1], orbital_b.center[1], 0, orbital_a.center[1]);
    double polynomial_1[MAX_POLYNOMIAL_SIZE], polynomial_2[MAX_POLYNOMIAL_SIZE];
    alt_little_n(0, 1, orbital_a.expC[0], orbital_b.expC[0], orbital_a.center[2], orbital_b.center[2], orbital_a.center[2], polynomial_1);
    // printf("Results: %lf\n", results);

    // printf("----------------------\n");

    alt_little_n(1, 0, orbital_a.expC[0], orbital_b.expC[0], orbital_a.center[2], orbital_b.center[2], orbital_a.center[2], polynomial_2);
    polynomial_2[0]+=(orbital_a.center[2] - orbital_b.center[2]);

    // printf("polynomial 1 and 2, columnwise.\n");
    // for (int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
    //     printf("%lf     %lf\n", polynomial_1[i], polynomial_2[i]);
    // }

    results = N_e_attraction(0,0, orbital_a, orbital_b, orbital_a.center);
    // results = boys_func(0, 0, 0, orbital_a, orbital_b, orbital_a.center);
    // results=chebychev_integral_boys(0,0,orbital_a, orbital_b,orbital_a.center);
    // double polynomial_terms[2]={-0.4867, -0.71843};
    // results = hyp1f1_int_boys(polynomial_terms, orbital_a.expC[0], orbital_b.expC[0], orbital_a, orbital_b, orbital_a.center); //This works!!

    printf("Results: %lf\n", results);

    return 0;
}

void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer){

    if (NULL == geom_pointer){
        printf("file can't be opened\n");
        exit(1);
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
        exit(1);
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
    //print orbital info
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

double dist_squared(double coords_A[num_dimensions], double coords_B[num_dimensions]){
    double dist_squared = 0;   
    for (int i = 0; i < num_dimensions; i++){
        dist_squared += pow(coords_A[i] - coords_B[i], 2);
    }
    return dist_squared;
}

double dot_product(double coords_A[num_dimensions], double coords_B[num_dimensions]){
    double DP = 0;
    for (int i = 0; i < num_dimensions; i++){
        DP += coords_A[i]*coords_B[i];
    }
    return DP;
}

void scalar_mult(double *coords_pt, double scalar, int size){ //uses pointers to multiply vector by scalar.
    for(int i = 0; i < size; i++){
        *(coords_pt + i) = *(coords_pt + i) * scalar;
    }
}

int fact2(int n){ //like normal factorial but decriment by 2 instead of 1.
    int fact = 1;
    for (int i = n; i >= -1; i -= 2){
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
        fact[i] = fact2(part);
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

        // finally, calculate the normalization constant for each primitive.
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

    // printf("%d, %d, %lf, %lf\n", ang_coord_a, ang_coord_b, center_a_coord, center_b_coord);
    if (ang_coord_a == 0 && ang_coord_b == 0){ //definition part 1
        // printf("    s(0,0)=1\n");
        return 1;
    }

    sum_ab = alpha + beta;
    P = (alpha*center_a_coord + beta*center_b_coord) / sum_ab;

    if(ang_coord_a == 1 && ang_coord_b == 0){ //definition part 2
        // printf("    s(1,0)=%lf\n", -(center_a_coord - P));
        return -(center_a_coord - P);
    }
    if(ang_coord_a != 1 && ang_coord_b == 0){ //recurrence relation
        double s_prev_a, s_prev_2_a;
        // printf("s(%d,%d) requires recurrance\n", ang_coord_a, ang_coord_b);
        s_prev_a = little_s(ang_coord_a - 1, 0, alpha, beta, center_a_coord, center_b_coord);
        s_prev_2_a = little_s(ang_coord_a - 2, 0, alpha, beta, center_a_coord, center_b_coord);
        double value = -(center_a_coord - P) * s_prev_a + ((ang_coord_a - 1) / (2*sum_ab)) * s_prev_2_a;
        // printf("s(%d,%d) is %lf\n", ang_coord_a, ang_coord_b, value);
        return value;
    }
    if(ang_coord_b != 0){ //transfer
        double s_xfer, s_prev_b;
        // printf("s(%d,%d) requires transfer\n", ang_coord_a,ang_coord_b);
        s_xfer = little_s(ang_coord_a + 1, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord);
        s_prev_b = little_s(ang_coord_a, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord);
        double value = s_xfer + (center_a_coord - center_b_coord) * s_prev_b; 
        // printf("s(%d,%d) is %lf\n", ang_coord_a, ang_coord_b, value);
        return value;
    }
}

double primitive_overlap(int dim_a, int dim_b, struct Orbital orbital_a, struct Orbital orbital_b){
    double alpha = orbital_a.expC[dim_a];
    double beta = orbital_b.expC[dim_b];
    double al_bet = alpha + beta;
    double EAB, exponent, dist_squrd, Overlap;
    
    dist_squrd = dist_squared(orbital_a.center, orbital_b.center);
    exponent = -(alpha * beta / al_bet) * dist_squrd;
    EAB = exp(exponent);

    Overlap = EAB * pow((M_PI / al_bet), 1.5);
    // printf("Overlap Coeff is %lf\n", Overlap);

    for (int i = 0; i < num_dimensions; i++){
        Overlap *= little_s(orbital_a.angular_momentum_vector[i], orbital_b.angular_momentum_vector[i], alpha, beta, orbital_a.center[i], orbital_b.center[i]);
        // printf("After calculating s(%d), the Overlap Integral is %lf\n", i, Overlap);
    }

    return Overlap;
}

double orbital_overlap(struct Orbital orbital_a, struct Orbital orbital_b){
    //every combination of primitive overlaps for 2 given orbitals.
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
    printf("Overlap matrix:\n");
    for(int i = 0; i < BS_Size; i++){
        for(int j = 0; j < BS_Size; j++){
            overlap_matrix[i][j] = orbital_overlap(used_orbitals[i], used_orbitals[j]);
            printf("    %lf", overlap_matrix[i][j]);
        }
        printf("\n");
    }
}

//The x component of the kinetic energy integral can be written as a sum of one dimensional overlap integrals.
//Recall that S_x=constants*s_x(ax,bx)*s_y(ay,by)*s_z(az,bz)
//K_x=constants*k_x(ax,bx)sy(az,by)sz(az,bz)
//The value of k() can be obtained from s() with other factors since k is really the sum of a polynomial expansion of s() functions

double little_k(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord){
    //The sum of overlap integrals that have varying angular momentum components obtained from taking the 
    // derivative of the gaussian as the KE operator requires.
    //use the same arguments as little_s, except when changing the angular momentum numbers.
    
    // printf("%d %d %lf %lf %lf %lf\n", ang_coord_a, ang_coord_b, alpha, beta, center_a_coord, center_b_coord);

    // initial conditions
    if(ang_coord_a == 0 && ang_coord_b == 0){
        double ls = little_s(1, 1, alpha, beta, center_a_coord, center_b_coord);
        // printf("k(0,0) needs s(1,1) which is %lf\n",ls);
        return 2 * alpha * beta * ls;
    }
    if(ang_coord_a > 0 && ang_coord_b == 0 ){
        // printf("k(%d,%d) needs extra little s'\n", ang_coord_a, ang_coord_b);
        double a_up = little_s(ang_coord_a+1, 1, alpha, beta, center_a_coord, center_b_coord);
        double a_down = little_s(ang_coord_a-1, 1, alpha, beta, center_a_coord, center_b_coord);
        // printf("k(%d,%d) little s' %lf %lf\n",ang_coord_a, ang_coord_b, a_up, a_down);
        return -ang_coord_a*beta*a_down + 2*alpha*beta*a_up;
    }
    if (ang_coord_b > 0 && ang_coord_a == 0){
        // printf("k(%d,%d) needs extra little s'\n", ang_coord_a, ang_coord_b);
        double b_down = little_s(1, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord);
        double b_up = little_s(1, ang_coord_b+1, alpha, beta, center_a_coord, center_b_coord);
        // printf("k(%d,%d) little s' %lf %lf\n",ang_coord_a, ang_coord_b, b_down, b_up);
        return -ang_coord_b*alpha*b_down + 2*alpha*beta*b_up;  
    }

    if(ang_coord_a < 0 || ang_coord_b < 0){
        printf("Bad angular momentum vector. Components need to be positive.");
        exit(1);
    }

    //general case. Protected from 2 negatives triggering by above exit command.
    if(ang_coord_a * ang_coord_b > 0){
        double s_lower = little_s(ang_coord_a-1, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord);
        double s_down_a = little_s(ang_coord_a-1, ang_coord_b+1, alpha, beta, center_a_coord, center_b_coord);
        double s_up_a = little_s(ang_coord_a+1, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord);
        double s_upper = little_s(ang_coord_a+1, ang_coord_b+1, alpha, beta, center_a_coord, center_b_coord);
        return (ang_coord_b*ang_coord_a*s_lower - 2*ang_coord_a*beta*s_down_a - 2*alpha*ang_coord_b*s_up_a + 4*alpha*beta*s_upper)/2;
    }
}

double primitives_KE(int primitive_idx_a, int primitive_idx_b, struct Orbital orbital_a, struct Orbital orbital_b){
    //combine all cartesian components of the KE integral and include the prefactors
    //IE K_x=constants*k_x(ax,bx)sy(az,by)sz(az,bz)
    //This is still for one set of primitives so alpha and beta are shared all the way down.
    double alpha = orbital_a.expC[primitive_idx_a];
    double beta = orbital_b.expC[primitive_idx_b];
    double sum_ab = alpha + beta;
    //The coordinates, however, are different.
    double EAB, pi_coeff, KE;
    double little_integrals[3][2]; //3 little s and 3 little k

    for(int i = 0; i < num_dimensions; i++){
        little_integrals[i][0] = little_s(orbital_a.angular_momentum_vector[i], orbital_b.angular_momentum_vector[i], alpha, beta, orbital_a.center[i], orbital_b.center[i]);
        little_integrals[i][1] = little_k(orbital_a.angular_momentum_vector[i], orbital_b.angular_momentum_vector[i], alpha, beta, orbital_a.center[i], orbital_b.center[i]);
        // printf("Little k: %lf\n", little_integrals[i][1]);//. Little s: %lf\n", i, little_integrals[i][1], little_integrals[i][0]);
    }

    double sum = 0; //just declaring it and not setting to 0 first was bad news in some circumstances.
    for(int dim = 0; dim < num_dimensions; dim++){
        int dim_s_ii = (dim+1)%3;
        int dim_s_iii = (dim+2)%3;
        double k_i = little_integrals[dim][1]; //the little k
        double s_ii = little_integrals[dim_s_ii][0]; //the other dimensions' little s
        double s_iii = little_integrals[dim_s_iii][0];
        sum += k_i * s_ii * s_iii;
    }

    EAB = exp(-(alpha*beta/sum_ab)*dist_squared(orbital_a.center, orbital_b.center));
    pi_coeff = pow(M_PI/sum_ab, 1.5);

    KE = EAB * pi_coeff * sum;

    return KE;
}

double orbital_kinetic_energy_integral(struct Orbital orbital_a, struct Orbital orbital_b){
    //every combo of primitives with appropriate normalization constants.
    double KE = 0; //Not explicitly assigning to 0 breaks this.
    for(int i = 0; i < num_dimensions; i++){
        for(int j = 0; j < num_dimensions; j++){
            double prim_KE = primitives_KE(i, j, orbital_a, orbital_b);
            double prefactor = orbital_a.NormC[i]*orbital_b.NormC[j]*orbital_a.primC[i]*orbital_b.primC[j];
            KE += prefactor * prim_KE;
        }
    }
    return KE;
}

void Calc_BS_KE_Matrix(struct Orbital orbital_array[Num_Orbitals], double KE_matrix[BS_Size][BS_Size], int indicies[BS_Size]){
    //Since there are 9 orbitals in the sytem and 7 in the BS, exclude some orbitals to avoid repeats.
    struct Orbital used_orbitals[BS_Size];
    for (int i = 0; i < BS_Size; i++){
        int next_idx = indicies[i];
        used_orbitals[i] = orbital_array[next_idx];
    }
    printf("KE matrix:\n");
    for(int i = 0; i < BS_Size; i++){
        for(int j = 0; j < BS_Size; j++){
            KE_matrix[i][j] = orbital_kinetic_energy_integral(used_orbitals[i], used_orbitals[j]);
            printf("    %lf", KE_matrix[i][j]);
        }
        printf("\n");
    }
}

//nuclear attraction integrals can be rewritten as overlap integrals with ugly terms. We can still exploit this.
//We will need to solve the Boys integral. Guidence by Minhuey Ho's tutorial.
//NE = Eab * 2pi/p * INT_0,1((Px-t^2(Px-Rx))* exp(-Rt^2))dt
void foil_polynomials(double *polynomial_ptr_1, double *polynomial_ptr_2, double *result_ptr){
    double results[MAX_POLYNOMIAL_SIZE][MAX_POLYNOMIAL_SIZE]; //as we multiply the polynomials, we get as many polynomials as there are terms in the polynomial with the highest number of terms. Each goes up to MAX_POLYNOMIAL_SIZE terms though most terms are 0.
    for (int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
        // printf("i=%d\n",i);
        for (int j = 0; j < MAX_POLYNOMIAL_SIZE; j++){
            // printf("j=%d\n",j);
            if(i+j>=MAX_POLYNOMIAL_SIZE){
                // printf("avoiding illegal write\n");
                continue; //don't try to fill in illegal index. This may truncate the polynomial but not with our small numbers.
            }
            // printf("%lf * %lf\n", *(polynomial_ptr_1 + i) , *(polynomial_ptr_2 + j));
            results[i+j][i] = *(polynomial_ptr_1 + i) * *(polynomial_ptr_2 + j);
        }
    }
    //now we combine like terms. for all the polynomials

    for(int j = 0; j < MAX_POLYNOMIAL_SIZE; j++){ //cols
        for(int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){ //rows
            // printf("results[%d][%d] %lf\n", j,i,results[j][i]);
            *(result_ptr+j) += results[j][i];
        }
    }
    // for(int j=0; j<MAX_POLYNOMIAL_SIZE; j++){
    //     printf("%lf\n", *(result_ptr+j));
    // }
}

void alt_little_n(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord, double nuc_coord, double *polynomial_pointer){
    //nuc-elec interaction of two gaussian primitives
    //modifies an array representing a polynomial that will be integrated by term later.
    //index i in the array represents i in t^(2i). the values are the coefficients of t.
    //initial conditions
    // printf("n(%d,%d)\n", ang_coord_a, ang_coord_b);
    // printf("params n(%d,%d), a=%lf, b=%lf, A=%lf, B=%lf, t=%lf, RR=%lf\n", ang_coord_a, ang_coord_b, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
    if(ang_coord_a == 0 && ang_coord_b == 0){
        // printf("n(0,0) = 1\n");
        *(polynomial_pointer) += 1;
        return;
    }
    // double tsqrd = t * t;
    double sum_ab = alpha + beta;
    double aA_bB = alpha * center_a_coord + beta * center_b_coord;
    double PC = aA_bB/sum_ab - nuc_coord;
    double coord_salad = PC - center_a_coord + nuc_coord;
    // printf("basic integrals %lf     %lf\n", basic_int_1, basic_int_2);
    if(ang_coord_a == 1 && ang_coord_b == 0){
        // printf("Basic solution with a=1 b=0\n");
        // printf("%lf     %lf\n", basic_int_1, basic_int_2);
        *(polynomial_pointer) += coord_salad; //add like terms to existing polynomial
        *(polynomial_pointer + 1) -= PC; //add like terms to existing polynomial
        return; //-center_a_coord + aA_bB/sum_ab  - PC * tsqrd;

    }
    double dummy_pol_1[MAX_POLYNOMIAL_SIZE], dummy_pol_2[MAX_POLYNOMIAL_SIZE]; //arrays to collect terms from recursive calls
    //recurrence index
    if(ang_coord_a > 1 && ang_coord_b == 0 ){
        double sum_1[MAX_POLYNOMIAL_SIZE];
        // printf("recurrence\n");
        // printf("n(%d,%d) needs extra little n'\n", ang_coord_a, ang_coord_b);
        // double a_down = alt_little_n((ang_coord_a - 1), 0, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
        // double a_down2 = alt_little_n((ang_coord_a - 1)-1, 0, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
        alt_little_n(ang_coord_a - 1, 0, alpha, beta, center_a_coord, center_b_coord, nuc_coord, dummy_pol_1);
        alt_little_n(ang_coord_a - 2, 0, alpha, beta, center_a_coord, center_b_coord, nuc_coord, dummy_pol_2);
        double adown_q2 = (ang_coord_a - 1) / (2 * sum_ab);
        // double a_downPC = a_down * PC;
        // return (adown_q2 * a_down2 + a_down * PC) - tsqrd * (adown_q2 * a_down2 + a_down * PC);

        scalar_mult(dummy_pol_1, coord_salad, MAX_POLYNOMIAL_SIZE);
        scalar_mult(dummy_pol_2, adown_q2, MAX_POLYNOMIAL_SIZE);

        //tsquared term needs to multiply the adown polynomial by PC not coord_salad. So after this sum we find the new product
        for(int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
            // dummy_pol_1[i]=dummy_pol_1[i]+ dummy_pol_2[i];
            sum_1[i] = dummy_pol_1[i] + dummy_pol_2[i];
        }

        scalar_mult(dummy_pol_1, PC/coord_salad, MAX_POLYNOMIAL_SIZE);

        //get the second sum now
        for(int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
            dummy_pol_1[i] += dummy_pol_2[i];
        }

        for (int i = MAX_POLYNOMIAL_SIZE-2; i >= 0; i--){ //start at one slot before the last.            
            // use the first copy as is. need to "multiply each term by t^2" in the second copy.
            dummy_pol_1[i+1] = dummy_pol_1[i];
            // printf("%lf     %lf\n", dummy_pol_1[i], sum_1[i]);
            //combine with like terms in the output array.
            *(polynomial_pointer + i) += sum_1[i] - dummy_pol_1[i];
                                    //unmodified sum  -    sum * t^2
        }
        // printf("n(%d,%d) little n's %lf %lf\n",ang_coord_a, ang_coord_b, a_down, a_down2);
        // return ((ang_coord_a - 1) * a_down2 * (1-tsqrd)) / q2 + 
        //     a_down * (PC - PC * tsqrd);
        return;
    }
    //transfer equation. Fallback if other options not hit.
    if (ang_coord_b > 0){
        for (int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
            dummy_pol_1[i]=0;
            dummy_pol_2[i]=0;
        }
        // printf("transfer\n");
        // printf("n(%d,%d) needs extra little n'\n", ang_coord_a, ang_coord_b);
        // double aup_bdown = alt_little_n(ang_coord_a+1, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord, nuc_coord, dummy_pol_1);
        // double b_down = alt_little_n(ang_coord_a, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord, nuc_coord, dummy_pol_2);
        alt_little_n(ang_coord_a+1, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord, nuc_coord, dummy_pol_1);
        alt_little_n(ang_coord_a, ang_coord_b-1, alpha, beta, center_a_coord, center_b_coord, nuc_coord, dummy_pol_2);
        double center_diff = center_a_coord - center_b_coord;
        // printf("n(%d,%d) little n's %lf %lf\n", ang_coord_a, ang_coord_b, aup_bdown, b_down);
        scalar_mult(dummy_pol_2, center_diff, MAX_POLYNOMIAL_SIZE);
        // printf("")
        for (int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
            // printf("%lf     %lf\n", dummy_pol_1[i], dummy_pol_2[i]);
            *(polynomial_pointer + i) += dummy_pol_1[i] + dummy_pol_2[i];
        }
        return; // aup_bdown + center_diff * b_down;  
    }
    // in case of bad inputs.
    if(ang_coord_a < 0 || ang_coord_b < 0){
        printf("Bad angular momentum vector. Components need to be positive.");
        return;
        // exit(1);
    }
}

double hyp1f1_clone(double a, double b, double x){
    double term = 1.0;
    double result = 1.0;
    int k = 0;
    while (fabs(term) > EPS * fabs(result)){
        //since each term is almost what is needed for the next, why recalculate all from scratch?
        term *= (a+k) * x / (b+k) / (k+1); //extending the pochhammer and factorial parts of the kth term. 
        result += term; //continue the sum
        k++;
    }
    return result;
}

double boys_func_hyp(double n, double T){
    // printf("n %lf, T %lf\n", n, T);
    return hyp1f1_clone(n+0.5, n+1.5, -T) / (2*n + 1);
}

double hyp1f1_int_boys(double polynomial_terms[MAX_POLYNOMIAL_SIZE], double alpha, double beta, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]){
    double integral = 0.00;
    double P[num_dimensions];

    for (int i = 0; i < num_dimensions; i++){
        P[i] = orbital_a.center[i]*alpha + orbital_b.center[i]*beta;
        // printf("P[%d]=%lf\n",i,P[i]);
    }

    scalar_mult(P, 1/(alpha+beta), num_dimensions);
    double T = (alpha + beta) * dist_squared(P, nuc_coords);

    for(int i=0; i < MAX_POLYNOMIAL_SIZE; i++){
        integral += polynomial_terms[i] * boys_func_hyp(i, T);
    }
    // printf("%lf\n", integral);
    return integral;
}

double N_e_attraction(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]){
    double alpha = orbital_a.expC[exp_a];
    double beta = orbital_b.expC[exp_b];
    double sum_ab = alpha + beta;
    double EAB = exp(-(alpha * beta)/sum_ab * dist_squared(orbital_a.center, orbital_b.center));
    double polynomial[MAX_POLYNOMIAL_SIZE]={0,0,0,0,0,0,0,0,0,0};

    alt_little_n(0, 1, alpha, beta, orbital_a.center[2], orbital_b.center[2], orbital_a.center[2], polynomial);
    
    // for (int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
    //     printf("%lf\n",polynomial[i]);
    // }

    return EAB * (2 * M_PI / sum_ab) * hyp1f1_int_boys(polynomial, alpha, beta, orbital_a, orbital_b, nuc_coords);
}