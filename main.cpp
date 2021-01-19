#include<iostream>
#include<vector>
#include "reconstruction.h"
using namespace std;
int main() {
    int k=3;
    vector<double> coeff(k);
    for(int r=0;r<3;++r) {
        recon_coef_u_(&r, &k, coeff.data());
        for(int j=0;j<k;++j) {
            printf("%f,",coeff[j]);
        }
        printf("\n");
    }
}