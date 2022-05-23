#ifndef _GENOTYPE_READER_WRITER_H
#define _GENOTYPE_READER_WRITER_H

#include <utils/otools.h>

class genotype_reader_writer{
public:
    //CONSTRUCTOR/DESCTRUCTOR
    genotype_reader_writer();
    ~genotype_reader_writer();

    //ROUTINES
    //void readAndWriteGenotypes(string fvcfin, string fvcfout, string region);
    void readAndWriteGenotypes(string fvcfin, string fvcfout, string region, set <string> &founders, vector <string> &offspring, vector <unsigned int> &parents_idxs, vector < vector <bool> > &parents_haplotypes, int lower_limit, int upper_limi);

};

#endif