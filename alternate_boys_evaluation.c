#include <stdio.h>
#include <math.h>
#define EPS 10E-10
double pochhammer(double a, double n);
double hyp1f1_clone(double a, double b, double x);

int main(){
    double result = hyp1f1_clone(1,2,3);
    printf("%lf\n",result);
    return 0;
}

double hyp1f1_clone(double a, double b, double x){
    double term = 1.0;
    double result = 1.0;
    for (int k = 0; k < 500; k++){
        //since each term is almost what is needed for the next, why recalculate all from scratch?
        term *= (a+k) * x / (b+k) / (k+1); //extending the pochhammer and factorial parts of the kth term. 
        result += term; //continue the sum
        if(fabs(term) <= EPS * fabs(result)){ //i.e. is this change extremely small?
            break;
        }
    }
    return result;
}

double pochhammer(double a, double n){
    if(n == 1){
        return 1.0;
    } else {
        return (a + n - 1) * pochhammer(a, n - 1);
    }
}

