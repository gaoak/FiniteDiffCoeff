#include<iostream>
#include<vector>
#include "reconstruction.h"
using namespace std;
int testcoefficient() {
    vector<double> coeff(10);
    int k;
    printf("=====================\n");
    k=1;
    for(int r=-1;r<k;++r) {
        recon_coef_u_(&r, &k, coeff.data());
        for(int j=0;j<k;++j) {
            printf("%f,",coeff[j]);
        }
        printf("\n");
    }
    printf("=====================\n");
    k=2;
    for(int r=-1;r<k;++r) {
        recon_coef_u_(&r, &k, coeff.data());
        for(int j=0;j<k;++j) {
            printf("%f,",coeff[j]);
        }
        printf("\n");
    }
    printf("=====================\n");
    k=3;
    for(int r=-1;r<k;++r) {
        recon_coef_u_(&r, &k, coeff.data());
        for(int j=0;j<k;++j) {
            printf("%f,",coeff[j]);
        }
        printf("\n");
    }
    printf("=====================\n");
    k=4;
    for(int r=-1;r<k;++r) {
        recon_coef_u_(&r, &k, coeff.data());
        for(int j=0;j<k;++j) {
            printf("%f,",coeff[j]);
        }
        printf("\n");
    }
    printf("=====================\n");
    return 0;
}

int main() {
    ;
}