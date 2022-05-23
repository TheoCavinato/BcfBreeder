#ifndef _RECOMBINATION_H
#define _RECOMBINATION_H

#include <utils/otools.h>
#include <utils/random_number.h>

class recombination {
public:
    //DATA
    vector <int> bcfPosBp;
    vector <double> bcfRecRateBtwPos;
    vector <vector <bool>> haploToSelect;

    //CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
    recombination();
    ~recombination();

    //ROUTINES
    void readBcfPos(string fvcf, int lower_limit, int upper_limit);
    void bpToRecRate(vector <int> &gmap_bp_pos, vector <double> &gmap_cm_pos);
    void simulateRecombination(int n_samples, string recvalid);
};

#endif