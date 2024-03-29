#include <tool/tool_header.h>

#define MODE_MAF	0
#define MODE_MAC	1
#define MODE_AF		2
#define MODE_AC		3

void tool::runMainTask() {
	vrb.title("Compute main TASK");

	int n_offspring = PED.offspring.size();

	//simulate recombination sites
	REC_SITES.readBcfPos(options["vcf"].as < string > (), GMAP.pos_bp[0], GMAP.pos_bp.back());
	REC_SITES.bpToRecRate(GMAP.pos_bp, GMAP.pos_cm);

	string valid_path="None";
	string phase_path="None";
	if (options.count("recvalid")) valid_path = options["recvalid"].as < string > ();
	if (options.count("out-phasing")) phase_path = options["out-phasing"].as < string > ();
	REC_SITES.simulateRecombination(n_offspring, valid_path, phase_path);

	//read bcf and write simulated individuals
	GEN.readAndWriteGenotypes(options["vcf"].as < string > (), options["output"].as <string> (), options["region"].as < string > (), PED.founders, PED.offspring, PED.parents_idxs, REC_SITES.haploToSelect, GMAP.pos_bp[0], GMAP.pos_bp.back());



}
