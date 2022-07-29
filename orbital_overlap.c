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
int fact2(int n);
double dist_squared(double coords_A[num_dimensions], double coords_B[num_dimensions]);
void scalar_mult(double *coords_pt, double scalar);
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
double chebychev_integral_boys(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]);
double N_e_attraction(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]);

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
    Calc_BS_OV_Matrix(orbital_array, BS_overlap_matrix, included_indicies);

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
    Calc_BS_KE_Matrix(orbital_array, BS_overlap_matrix, included_indicies);

    //testing NE functions
    double results;
    // results = little_n(0, 0, orbital_a.expC[0], orbital_b.expC[0], orbital_a.center[0], orbital_b.center[0], 0, orbital_a.center[0]);
    // results = exp(-(orbital_a.expC[0] * orbital_b.expC[0])/(orbital_a.expC[0] + orbital_b.expC[0]) * dist_squared(orbital_a.center, orbital_b.center))
    //     * (2 * M_PI / (orbital_a.expC[0] + orbital_b.expC[0]));
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

void scalar_mult(double *coords_pt, double scalar){ //uses pointers to multiply vector by scalar.
    for(int i = 0; i < 3; i++){
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
double abscissa(int n, int i){
    double i_pi = M_PI * i;
    double n_up = n + 1;
    return (n_up -2 * i)/n_up + (2/M_PI) * (1 + (2/3) * pow(sin(i_pi/n_up), 2)) * cos(i_pi / n_up) * sin(i_pi / n_up);
}

double omega(int n, int i){
    double n_up = n + 1;
    return 16/3/n_up * pow(sin(M_PI * i / n_up), 4);
}

double little_n(int ang_coord_a, int ang_coord_b, double alpha, double beta, double center_a_coord, double center_b_coord, double t, double nuc_coord){
    //nuc-elec interaction of two gaussian primitives
    //initial conditions
    if(ang_coord_a == 0 && ang_coord_b == 0){
        // printf("n(0,0) = 1\n");
        return 1;
    }
    double sum_ab = alpha + beta;
    double aA_bB = alpha * center_a_coord + beta * center_b_coord;
    double basic_int_1 = -(center_a_coord - (aA_bB / sum_ab));
    double basic_int_2 = pow(t, 2) * ((aA_bB / sum_ab) - nuc_coord);

    if(ang_coord_a == 1 && ang_coord_b == 0){
        return basic_int_1 + basic_int_2;
    }
    //recurrence index
    if(ang_coord_a > 0 && ang_coord_b == 0 ){
        // printf("n(%d,%d) needs extra little n'\n", ang_coord_a, ang_coord_b);
        double a_down = little_n(ang_coord_a-1, 0, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
        double a_down2 = little_n(ang_coord_a-2, 0, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
        // printf("n(%d,%d) little n's %lf %lf\n",ang_coord_a, ang_coord_b, a_down, a_down2);
        return basic_int_1 + basic_int_2 * a_down + ((ang_coord_a - 1)/(2 * sum_ab)) * (1 - pow(t, 2)) * a_down2;
    }
    //transfer equation
    if (ang_coord_a > 0 && ang_coord_b > 0){
        // printf("n(%d,%d) needs extra little n'\n", ang_coord_a, ang_coord_b);
        double aup_bdown = little_n(ang_coord_a+1, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
        double b_down = little_n(ang_coord_a, ang_coord_b - 1, alpha, beta, center_a_coord, center_b_coord, t, nuc_coord);
        // printf("n(%d,%d) little n's %lf %lf\n",ang_coord_a, ang_coord_b, aup_bdown, b_down);
        return aup_bdown + (center_a_coord - center_b_coord) * b_down;  
    }

    if(ang_coord_a < 0 || ang_coord_b < 0){
        printf("Bad angular momentum vector. Components need to be positive.");
        exit(1);
    }
}


double boys_func(double x, int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]){
    double alpha = orbital_a.expC[exp_a];
    double beta = orbital_b.expC[exp_b];
    double sum_ab = alpha + beta;
    double t = (x+1)/2;

    double dummy_a[num_dimensions], dummy_b[num_dimensions], aA_bB[num_dimensions];
    for (int i = 0; i < num_dimensions; i++){
        dummy_a[i] = orbital_a.center[i];
        dummy_b[i] = orbital_b.center[i];
    }

    scalar_mult(dummy_a, alpha); // aA
    scalar_mult(dummy_b, beta); // bB

    for(int i = 0; i < num_dimensions; i++){
        aA_bB[i] = dummy_a[i] + dummy_b[i]; //vector addition
    }
    scalar_mult(aA_bB, 1/sum_ab); //scale the vector
    
    double product = 1;
    for(int i = 0; i < num_dimensions; i++){
        aA_bB[i] = aA_bB[i] - nuc_coords[i]; //vector subtraction.
        product *= little_n(orbital_a.angular_momentum_vector[i], orbital_b.angular_momentum_vector[i], alpha, beta, orbital_a.center[i], orbital_b.center[i], t, nuc_coords[i]);
    }

    double big_const = dist_squared(aA_bB, aA_bB); //((alpha * A_coords + beta * B_coords) / (sum_ab)) - nuc_coord dotted into itself (i.e. dist squared.)
    double simple_part = 0.5 * exp(-(sum_ab * pow(t, 2) * big_const));

    return simple_part * product;
}

double chebychev_integral_boys(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]){ //ideally pass boys function as arg but we only need to solve this integral so it is just baked into the function.
    //Minhuey's paper uses https://doi.org/10.1016/0010-4655(93)90035-B as a guide.
    //set some parameters
    double tolerance = pow(10, -10);
    int num_points = 50000;
    double c0 = cos(M_PI/6);
    double s0 = cos(M_PI/6);
    double c1, s1, q, p, chp, j, c, s, xp;
    float err = 10;
    int n = 3;

    double bf_a_plus = boys_func(abscissa(2,1), exp_a, exp_b, orbital_a, orbital_b, nuc_coords);
    double bf_a_minus = boys_func(-abscissa(2,1), exp_a, exp_b, orbital_a, orbital_b, nuc_coords);
    
    q = bf_a_plus * bf_a_minus * omega(2,1);
    p = boys_func(0.0, exp_a, exp_b, orbital_a, orbital_b, nuc_coords);

    chp = p + q; //value of integral
    j = 0; //oscilates between 0 and 1.
    c1 = s0;
    s1 = c0;
    while(tolerance < err && (2* n * (1 - j) + j * 4 * n / 3 - 1) <= num_points){
        j = 1 - j;
        c1 = j * c1 + (1 - j) * c0;
        s1 = j * s1 + (1 - j) * s0;
        c0 = j * c0 + (1 - j) * sqrt((1 + c0) * 0.5);
        s0 = j * s0 + (1 - j) * s0 / (c0 + c0);
        c = c0;
        s = s0;
        for(int i = 1; i < n; i += 2){
            xp = s;
            s = s * c1 + c * s1;
            c = c * c1 - xp * s1;        
        }
        n = (1 + j) * n;
        p = p + (1 - j) * (chp - q);
        err = 16 * fabs((1 - j) * (q - 3*p/2) + j * (chp -2*q)) / (3*n); //error estimates
        q = (1 - j)*q + j*chp;
        printf("tolerance < err && expression <= numpoints\n");
        printf("%lf   %lf           %lf         %d\n", tolerance, err, (2* n * (1 - j) + j * 4 * n / 3 - 1), num_points);
    }
    return 16 * q / (3 * n);
}

double N_e_attraction(int exp_a, int exp_b, struct Orbital orbital_a, struct Orbital orbital_b, double nuc_coords[num_dimensions]){
    double alpha = orbital_a.expC[exp_a];
    double beta = orbital_b.expC[exp_b];
    double sum_ab = alpha + beta;

    double EAB = exp(-(alpha * beta)/sum_ab * dist_squared(orbital_a.center, orbital_b.center));
    // printf("EAB %lf\n", EAB);
    
    return EAB * (2 * M_PI / sum_ab) * chebychev_integral_boys(exp_a, exp_b, orbital_a, orbital_b, nuc_coords);
}