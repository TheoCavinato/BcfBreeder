#include <tool/tool_header.h>


void tool::read_files_and_initialise() {
	vrb.title("Initialization:");

	//step0: Initialize seed and multi-threading
	rng.setSeed(options["seed"].as < int > ());
	vrb.bullet("Seed	: " + to_string(rng.getSeed()));

	//step2: Read input files
	//read pedigree
	PED.readPed(options["ped"].as < string > ());
	PED.calculate_idxs();
	int n_offspring = PED.offspring.size();

	//read the genetic map for the correct chromosome -> OK
	GMAP.readGeneticMapFile(options["map"].as < string > ());

}
