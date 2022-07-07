#include <stdio.h>
#include <math.h>

#define MAXATOMS 20

struct Atom {
    int Z;
    double atom_x, atom_y, atom_z;
};

// void add_atom_details(struct Atom *atom, int nuc_charge, float x, float y, float z);

int main() {   
    struct Atom atom_array[MAXATOMS];
    int atom_idx = 0;

    FILE *file_pointer;
    file_pointer = fopen("geom.xyz", "r");

    if (NULL == file_pointer) {
        printf("file can't be opened \n");
    }

    while (!feof(file_pointer)) {
        fscanf(file_pointer, "%d %lf %lf %lf", 
            &atom_array[atom_idx].Z, &atom_array[atom_idx].atom_x, 
            &atom_array[atom_idx].atom_y, &atom_array[atom_idx].atom_z);
        atom_idx++;
    }

    fclose(file_pointer);
    
    double nuclear_repulsion = 0;
    for(int i = 0; i < atom_idx - 1; i++){
        //calculate distance to atoms with greater index.
        for(int j = i+1; j < atom_idx; j++){
            double distance = sqrt(
               pow(atom_array[i].atom_x - atom_array[j].atom_x, 2)
             + pow(atom_array[i].atom_y - atom_array[j].atom_y, 2)
             + pow(atom_array[i].atom_z - atom_array[j].atom_z, 2)
            );
        
            nuclear_repulsion += atom_array[i].Z * atom_array[j].Z / distance;
            printf("The nuclear distance between atom %d and %d is %lf, increasing nuclear repulsion to %lf\n",
                i, j, distance, nuclear_repulsion);
        }
        // printf("Atom %d has charge %d coords %lf %lf %lf\n", i, 
        // atom_array[i].Z, atom_array[i].atom_z, atom_array[i].atom_y, atom_array[i].atom_z);
    } 
    
}

// void add_atom_details(struct Atom *atom, int nuc_charge, float x, float y, float z) {
//     atom->Z = nuc_charge;
//     atom->atom_x = x;
//     atom->atom_y = y;
//     atom->atom_z = z;
// }
