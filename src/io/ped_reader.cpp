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

#include <io/ped_reader.h>

ped_reader::ped_reader(){

}

ped_reader::~ped_reader(){

}

void ped_reader::readPed(string ped_path){
    tac.clock();

	string buffer;
	int line = 0;
    vector <string> tokens;
	input_file fd_ped(ped_path);
	if (fd_ped.fail()) vrb.error("Cannot open ped file");
	while (getline(fd_ped, buffer, '\n')) {
		stb.split(buffer, tokens);
        if (tokens[2]!="NA" && tokens[3]!="NA" ) {
            //offspring.insert(tokens[1]);
            parents.insert(tokens[2]);
            parents.insert(tokens[3]);
            //raw_offspring.push_back(tokens[1]);
            offspring.push_back(tokens[1]);
            raw_parents.push_back(tokens[2]);
            raw_parents.push_back(tokens[3]);
        }
        line++;
    }
    fd_ped.close();
	vrb.bullet("PED parsing [n=" + stb.str(line) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void ped_reader::calculate_idxs(){

    //step 1: order as in the output file (founders at the beginning, offpsring a the end)
    for(string i : parents) {
        if (!count(offspring.begin(), offspring.end(), i)){
            founders.insert(i);
        }
    }

    vector <string> output_line;
    for(string i : founders) output_line.push_back(i);
    for(string i : offspring) output_line.push_back(i);

    //step 2: find idx of the parent of each offspring in output_line
    string child_name, father_name, mother_name;
    std::vector<string>::iterator it_child, it_father, it_mother;
    unsigned int idx_child, idx_father, idx_mother;
    //populate parents_idxs with 0s
    for (int i = 0; i<offspring.size()*2; i++) parents_idxs.push_back(0);

    for(int i=0; i<offspring.size(); i++){
        // get raw data positions
        child_name = offspring[i];
        father_name = raw_parents[i*2];
        mother_name = raw_parents[i*2+1];

        it_child = find(output_line.begin(), output_line.end(), child_name);
        it_father = find(output_line.begin(), output_line.end(), father_name);
        it_mother = find(output_line.begin(), output_line.end(), mother_name);

        idx_child = it_child-output_line.begin();
        idx_father = it_father-output_line.begin();
        idx_mother = it_mother-output_line.begin();

        parents_idxs[(idx_child-founders.size())*2] = idx_father;
        parents_idxs[(idx_child-founders.size())*2+1] = idx_mother;
    }

	vrb.bullet("PED calculation of idxs (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}