#include <stdio.h>

#define NUMATOMS 3
#define BS_Size 7
#define Num_Orbitals 9
#define num_dimensions 3

//store orbital info
struct Orbital {
    double primC, expC;
    int angular_momentum_vector[num_dimensions];
    double center[num_dimensions];
    int parent_atom_Z, parent_atom_idx;
} orbital_array[Num_Orbitals];

void orbital_info(struct Orbital orb, int idx);
void get_geom_details(struct Orbital orbital_array[BS_Size], FILE *geom_pointer);

int main(){
    //get info from files.
    //geometry
    FILE *geom_pointer;
    geom_pointer = fopen("water_wolfram.xyz", "r");
    get_geom_details(orbital_array, geom_pointer);

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
        printf("atom number %d\n", atom_idx);
        fscanf(geom_pointer, "%d %lf %lf %lf", 
            &orbital_array[orbital_idx].parent_atom_Z, &orbital_array[orbital_idx].center[0], 
            &orbital_array[orbital_idx].center[1], &orbital_array[orbital_idx].center[2]
        );

        printf("atom %d has atomic number %d\n",atom_idx, orbital_array[orbital_idx].parent_atom_Z);
        
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
    printf("\n");

    // check out all orbitals in detail
    for (int i = 0; i < Num_Orbitals; i++){
        orbital_info(orbital_array[i], i);
    }

}

void orbital_info(struct Orbital orb, int idx){
    printf("orbital %d has angular momentum vector: {", idx+1);
    for (int j = 0; j < num_dimensions; j++){
        printf(" %d", orb.angular_momentum_vector[j]);
    }
    printf(" }\n");
    printf("it is located at coordinates: (");
    for (int j = 0; j < num_dimensions; j++){
        printf(" %lf",orb.center[j]);

    }
    printf(" )\n");
    printf("it is on atom #%d, which has atomic number %d\n", orb.parent_atom_idx, orb.parent_atom_Z);
}
