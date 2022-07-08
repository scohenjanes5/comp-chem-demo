#include <stdio.h>

#define NUMATOMS 3
#define BS_Size 7
#define num_dimensions 3

struct Orbital {
    double primC, expC;
    int angular_momentum_vector[num_dimensions];
    double center[num_dimensions];
    int parent_atom_Z;
} orbital_array[BS_Size];

void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer);

int main(){
    //get info from files.
    //geometry
    FILE *geom_pointer;
    geom_pointer = fopen("water_wolfram.xyz", "r");
    get_geom_details(orbital_array, geom_pointer);

    // check to see angular momentum vectors are correct.
    for (int i = 0; i < BS_Size; i++){
        printf("orbital %d has angular momentum vector: {", i+1);
        for (int j = 0; j < num_dimensions; j++){
            printf(" %d",orbital_array[i].angular_momentum_vector[j]);
        }
        printf(" }\n");
    }

    return 0;
}

void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer){
    int orbital_idx = 0;
    if (NULL == geom_pointer){
        printf("file can't be opened\n");
    }

    int additional_orbitals;
    // Reading coords and atomic number from file
    while (!feof(geom_pointer)){
        fscanf(geom_pointer, "%d %lf %lf %lf", 
            &orbital_array[orbital_idx].parent_atom_Z, &orbital_array[orbital_idx].center[0], 
            &orbital_array[orbital_idx].center[1], &orbital_array[orbital_idx].center[2]
        );
        
        if(orbital_array[orbital_idx].parent_atom_Z > 2){
            //4 more orbitals for n=2
            additional_orbitals = 4;
        } else {
            additional_orbitals = 1; //just 2S for H
        }
        
        for(int i = 1; i <= additional_orbitals; i++){
            //copy info to orbitals on the same atom
            orbital_array[orbital_idx+i].parent_atom_Z = orbital_array[orbital_idx].parent_atom_Z;
            orbital_array[orbital_idx+i].center[0] = orbital_array[orbital_idx].center[0];
            orbital_array[orbital_idx+i].center[1] = orbital_array[orbital_idx].center[1];
            orbital_array[orbital_idx+i].center[2] = orbital_array[orbital_idx].center[2];

            if(i > 1){ //This only triggers for the p orbitals. Not worrying about d orbitals yet :)
                //correct the angular momentum vector
                orbital_array[orbital_idx+i].angular_momentum_vector[i-2] = 1;
            }
        }
        orbital_idx += additional_orbitals;
    }
}
