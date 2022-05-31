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

#ifndef _PED_READER_H
#define _PED_READER_H

#include <utils/otools.h>

class ped_reader{
public:

    //DATA
    vector <string> offspring;
    set <string> founders;
    vector <unsigned int> parents_idxs;

    //CONSTRUCTOR/DESTRUCTOR
    ped_reader();
    ~ped_reader();

    //ROUTINES
    void readPed(string ped_path);
    void calculate_idxs();

private:
    set <string> parents;
    vector <string> raw_parents;
    vector <string> raw_offspring;
};

#endif
