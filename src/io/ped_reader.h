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
