/*******************************************************************************
 * Copyright (C) 2020 Théo Cavinato, University of Lausanne
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

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
    void simulateRecombination(int n_samples, string recvalid, string phasevalid);
};

#endif