#include <stdio.h>
#define MAX_POLYNOMIAL_SIZE 4

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

int main(){
    
    double pol_a[4]={1.00, 1.00,0 ,0};
    double pol_b[4]={1.00, 1.00,0 ,0};
    double results[MAX_POLYNOMIAL_SIZE];
    foil_polynomials(pol_a,pol_b,results);
    for (int i = 0; i < MAX_POLYNOMIAL_SIZE; i++){
        printf("%lf ", results[i]);
    }
    printf("\n");
    return 0;
}