#include <stdio.h>

#define NUMATOMS 3
#define BS_Size 7
#define Num_Orbitals 9
#define num_dimensions 3

//store orbital info
struct Orbital {
    double primC[num_dimensions], expC[num_dimensions];
    int angular_momentum_vector[num_dimensions];
    double center[num_dimensions];
    int parent_atom_Z, parent_atom_idx;
} orbital_array[Num_Orbitals];

void orbital_info(struct Orbital orb, int idx);
void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer);
void get_coefs(struct Orbital orbital_array[BS_Size], FILE *coef_pointer);

int main(){
    //get info from files.
    //geometry
    FILE *geom_pointer, *coef_pointer;
    geom_pointer = fopen("water_wolfram.xyz", "r");
    coef_pointer = fopen("cont_exp_coefs.csv", "r");
    get_geom_details(orbital_array, geom_pointer);
    get_coefs(orbital_array, coef_pointer);
    // check out all orbitals in detail
    for (int i = 0; i < Num_Orbitals; i++){
        orbital_info(orbital_array[i], i);
    }

    // check to see angular momentum vectors are correct.
    // for (int i = 0; i < BS_Size; i++){
    //     printf("orbital %d has angular momentum vector: {", i+1);
    //     for (int j = 0; j < num_dimensions; j++){
    //         printf(" %d",orbital_array[i].angular_momentum_vector[j]);
    //     }
    //     printf(" }\n");
    // }

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
        if (orbital_array[i].parent_atom_Z == 8 && i % 4 == 0 && Num_Orbitals-i > 4){
            //be flexible so more atoms can be added later.
            for(int j = 0; j < num_dimensions; j++){
                orbital_array[i].primC[j] = basis_array[j + num_dimensions * 2][0]; //assign 1s orbital coefs here
                orbital_array[i].expC[j] = basis_array[j + num_dimensions * 2][1]; // ^ and here.
                orbital_array[i+1].primC[j] = basis_array[j + num_dimensions * 3][0]; //now assign 2s coefs here
                orbital_array[i+1].expC[j] = basis_array[j + num_dimensions * 3][1]; //and here too.
                orbital_array[i+2].primC[j] = basis_array[j + num_dimensions * 4][0]; //now assign 2px coefs here
                orbital_array[i+2].expC[j] = basis_array[j + num_dimensions * 4][1]; //and here too.
                orbital_array[i+3].primC[j] = basis_array[j + num_dimensions * 5][0]; //now assign 2py coefs here
                orbital_array[i+3].expC[j] = basis_array[j + num_dimensions * 5][1]; //and here too.
                orbital_array[i+4].primC[j] = basis_array[j + num_dimensions * 6][0]; //now assign 2pz coefs here
                orbital_array[i+4].expC[j] = basis_array[j + num_dimensions * 6][1]; //and here too.
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
        printf(" %lf",orb.center[j]);

    }
    printf(" )\n");
    printf("it is on atom #%d, which has atomic number %d\n\n----------------\n", orb.parent_atom_idx, orb.parent_atom_Z);
}
