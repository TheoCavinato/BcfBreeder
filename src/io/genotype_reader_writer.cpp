#include <io/genotype_reader_writer.h>

genotype_reader_writer::genotype_reader_writer(){}

genotype_reader_writer::~genotype_reader_writer(){}

void genotype_reader_writer::readAndWriteGenotypes(string fvcfin, string fvcfout, string region, set <string> &founders, vector <string> &offspring, vector <unsigned int> &parents_idxs, vector < vector <bool> > &parents_haplotypes, int lower_limit, int upper_limit){
	int n_samples_out = founders.size() + offspring.size();
	unsigned int n_variants=0;

	//-----------INITIALISE VCF TO READ---------------//
	//Initialize HTSlib VCF/BCF reader
	htsFile * fin= hts_open(fvcfin.c_str(),"r");
	bcf_hdr_t * hdr_in = bcf_hdr_read(fin);

    //Retrieve founders idxs
    int n_samples_in = bcf_hdr_nsamples(hdr_in);
	vector <unsigned int> founders_idxs;
	for (string founder : founders){
        for (int j =0; j<n_samples_in; j++){
            if (founder == string(hdr_in->samples[j])) {
                founders_idxs.push_back(j*2);
                founders_idxs.push_back(j*2+1);
            }
        }
    }

	//for (unsigned int i : founders_idxs) cout << i << " ";
	//cout << endl;
	//for (unsigned int i : parents_idxs) cout << i << " "; 
	//cout << endl;

	//To read VCF records
	unsigned int nset = 0;
	int ngt, *gt_arr = NULL, ngt_arr = 0;
	bcf1_t * rec_in = bcf_init();

	//Data to verbose every 10 % written
	double ten_percent_comparator = (float) parents_haplotypes.size() * 10 / 100;
	double ten_percent = ten_percent_comparator;
	int percentage = 10;

	//-----------INITIALISE VCF TO WRITE---------------//
	string fname = fvcfout;
	string file_format="w";
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; }
	bcf_hdr_t * hdr_out = bcf_hdr_init("w");
	htsFile * fout = hts_open(fname.c_str(),file_format.c_str());

	//Create VCF header
	bcf_hdr_append(hdr_out, string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr_out, "##source=shapeit4.1.3");
	bcf_hdr_append(hdr_out, string("##contig=<ID="+region+">").c_str());
	bcf_hdr_append(hdr_out, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr_out, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr_out, "##INFO=<ID=CM,Number=A,Type=Float,Description=\"Interpolated cM position\">");
	bcf_hdr_append(hdr_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	for(string founder : founders ){
		bcf_hdr_add_sample(hdr_out, founder.c_str());
	}
	for(string offspr : offspring){
		bcf_hdr_add_sample(hdr_out, offspr.c_str());
	}
	bcf_hdr_add_sample(hdr_out, NULL);      // to update internal structures
	if (bcf_hdr_write(fout, hdr_out) < 0) vrb.error("Failing to write VCF/header");

	//To write VCF record
	bcf1_t *rec_out = bcf_init();
	int * genotypes_out = (int*)malloc(bcf_hdr_nsamples(hdr_out)*2*sizeof(int));

	//-----------READ fvcfin AND WRITE fvcfout---------------//
	tac.clock();
	while((nset=bcf_read(fin, hdr_in, rec_in))==0) {
		bcf_clear(rec_out);
		bcf_unpack(rec_in, BCF_UN_STR);

        if (rec_in->pos > lower_limit && rec_in->pos < upper_limit ){ //to only read positions that are between in the gmap limits

			//Populate genotypes_out with founders -> validated
			bool a0;
			bool a1;
			unsigned int a0_founder_idx_in;
			unsigned int a1_founder_idx_in;
			int count_alt = 0;
			ngt = bcf_get_genotypes(hdr_in, rec_in, &gt_arr, &ngt_arr);
			for(int i = 0; i < founders_idxs.size(); i+=2){
				a0_founder_idx_in = founders_idxs[i];
				a1_founder_idx_in = founders_idxs[i+1];
				a0 = (bcf_gt_allele(gt_arr[a0_founder_idx_in])==1);
				a1 = (bcf_gt_allele(gt_arr[a1_founder_idx_in])==1);
				genotypes_out[i] = bcf_gt_phased(a0);
				genotypes_out[i+1] = bcf_gt_phased(a1);
				count_alt += a0+a1;
			}

			//Populate genotypes_out with individuals to simulate 
			unsigned int first_parent_idx_out;
			unsigned int second_parent_idx_out;
			bool first_haplo;
			bool second_haplo;
			int vector_itr = 0;
			for(int i = founders_idxs.size(); i<founders_idxs.size()+parents_idxs.size(); i+=2){
				first_parent_idx_out = parents_idxs[vector_itr]*2;
				second_parent_idx_out = parents_idxs[vector_itr+1]*2;

				first_haplo = parents_haplotypes[n_variants][vector_itr];
				second_haplo = parents_haplotypes[n_variants][vector_itr+1];

				//cout << "Output position " << i << " and " << i+1 <<
				//" | " << first_haplo << " " << second_haplo <<
				//" | " << first_parent_idx_out+first_haplo << " " << second_parent_idx_out+second_haplo << endl;

				a0 = (bcf_gt_allele(genotypes_out[first_parent_idx_out+first_haplo])==1);
				a1 = (bcf_gt_allele(genotypes_out[second_parent_idx_out+second_haplo])==1);

				genotypes_out[i] = bcf_gt_phased(a0);
				genotypes_out[i+1] = bcf_gt_phased(a1);
				count_alt += a0+a1;
				vector_itr+=2;
			}
			//cout << "--------------------------" << endl;

			//Add all the 9 first column to genotypes_out
			rec_out->pos = rec_in->pos;
			rec_in->rid = rec_in->rid;

			string ref = string(rec_in->d.allele[0]);
			string alt = string(rec_in->d.allele[1]);
			string alleles = ref + "," + alt;
			bcf_update_alleles_str(hdr_out, rec_out, alleles.c_str());

			bcf_update_info_int32(hdr_out, rec_out, "AC", &count_alt, 1);
			float freq_alt = count_alt * 1.0 / (2 * n_samples_out);
			bcf_update_info_float(hdr_out, rec_out, "AF", &freq_alt, 1);

			bcf_update_id(hdr_out, rec_out, rec_in->d.id);

			//Write genotypes_out to the fout
			bcf_update_genotypes(hdr_out, rec_out, genotypes_out, bcf_hdr_nsamples(hdr_out)*2);
			if (bcf_write(fout, hdr_out, rec_out) < 0) vrb.error("Failing to write VCF/record");

			//verbose each 10% of the file written
			n_variants++;

			if(n_variants > ten_percent_comparator){
				vrb.bullet("[ " + to_string(percentage) + " % of the file written (" +  stb.str(tac.rel_time()*1.0/1000, 2) + "s) ]");
				ten_percent_comparator+=ten_percent;
				percentage+=10;
				tac.clock();
			}
		}
	}

	//Free allocated memory
	free(genotypes_out);
	bcf_hdr_destroy(hdr_out);
	bcf_destroy(rec_out);

	free(gt_arr);

	if (hts_close(fin)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	if (hts_close(fout)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	vrb.bullet("[ " + to_string(percentage) + " % of the file written (" +  stb.str(tac.rel_time()*1.0/1000, 2) + "s) ]");
}
