Scripts for the analysis of metagenomic and genomic data searched against the eggNOG database of orthologous groups.

Directories:
/data/
contains output from HMMer. In order to run scripts in the project, the following must be downloaded from eggNOG:
bactNOG
bctoNOG
defNOG
firmNOG

The script get_phylum_nog_dicts.py generates the output in /data/phylum_nogs, which consists of python dictionaries,
saved in json format, detailing the NOGs present in each species within eggNOG.

Additionally, asf_hmmer_output and metagenome_hmmer_output contain raw HMMer output in table format.

/bin/

Contains the following scripts:
asf_coverage_per_sample.py:
	Requires: /data/bactNOG/bactNOG.annotations.tsv
		  /results/metagenome_samples_per_NOG.tsv
		  /results/metagenome_coverage/asf_mg_coverage.tsv
	Output:
		/results/mg_cover_by_asf.svg (boxplot with NOG coverage of each metagenome by the ASF)

asf_NOG_analysis.py:
	Requires:  /results/asf_unique_nogs.tsv
		   /results/asf_nogs_in_core_metagenome.tsv
	Output:
		/results/figure_6.svg (stacked bar plot of unique contribution of each ASF species in each NOG category)
		/results/figure_7.svg (stacked bar plot of unique contribution of each ASF species to core metagenome)
		/results/fig6_data.tsv
		/results/fig7_data.tsv

parse_hmmer_output.py: parses ASF/metagenomic HMMer output to write asf_nog_lists and metagenome_nog_lists in /results/
	Requires:  /data/asf_hammer_output/*_scan.tbl  (hmmer output table for each ASF species)
		   /data/metagenome_hmmer_output/ <-- MUST CONTAIN ONE SUBDIRECTORY FOR EACH SAMPLE. EACH SUBDIRECTORY MUST CONTAIN .tbl HMMER OUTPUT FILES
	Output:
		/results/asf_NOG_lists/all_nogs/*_nogs.txt  (for each ASF species)
		/results/metagenome_NOG_lists/*_nogs.txt  (for each metagenome)

get_phylum_nog_dicts.py: parses bctoNOG, firmNOG, and defNOG to save species:[NOG] dicts and save in json format.
	Requires:  /data/bactNOG/bactNOG.annotations.tsv
		   /results/metagenome_samples_per_NOG.tsv
		   /data/bctoNOG/bctoNOG.members.tsv
		   /data/firmNOG/firmNOG.members.tsv
		   /data/bactNOG/bactNOG.members.tsv
	Output:
		/data/phylum_nogs/bcto_nogs_dict.json
		/data/phylum_nogs/firm_nogs_dict.json

filter_for_metagenome_coverage.py: process species:[NOG] dicts to remove NOGs not in any metagenome, then generates
     coverage matrices for ASF, bctoNOG, and firmNOG.
	Requires:  /data/bactNOG/bactNOG.annotations.tsv
		   /results/metagenome_samples_per_NOG.tsv
		   /data/phylum_nogs/bcto_nogs_dict.json
		   /data/phylum_nogs/firm_nogs_dict.json
		   /results/asf_NOG_lists/asf_nog_dicts.json
	Output:
		/results/metagenome_coverage/asf_mg_coverage.tsv
		/results/metagenome_coverage/bcto_mg_coverage.tsv
		/results/metagenome_coverage/firm_mg_coverage.tsv
		/results/metagenome_coverage/asf_bcto_firm_mg_coverage.tsv

metagenome_nog_analysis.py: generates several figures to help interpret the metagenomic NOG data. Only figure S1 is included in publication.
	Requires:  /data/bactNOG/bactNOG.annotations.tsv
		   /results/metagenome_samples_per_NOG.tsv
	Output:
		/results/figure_S1.eps
		/results/figure_S2.eps
		/results/figure_S3.eps
		/results/figure_S4.eps
		/results/figure_S5.eps

get_unique_asf_nogs.py: finds and saves NOGs that are unique to each ASF species and creates:
			1. asf_nog_dicts.json: an ASF nog dictionary for filter_for_metagenome_coverage.py
			2. asf_nogs_in_core_metagenome.tsv: breakdown of ASF NOGs present in all metagenomes, for asf_nog_analysis
			3. asf_nog_dicts.json: dictionary of species:list_of_nogs
	Requires:  /data/bactNOG/bactNOG.annotations.tsv
		   /results/metagenome_samples_per_NOG.tsv
		   A file for each ASF species in /results/asf_NOG_lists/all_nogs/
	Output:
		/results/asf_unique_nogs.tsv'
		/results/asf_nogs_in_core_metagenome.tsv
		/results/asf_NOG_lists/asf_nog_dicts.json

perform_draws_and_coverage.py: uses coverage matrices to create random draws and compare metagenomic coverage of ASF vs. draws. This may take a while to run if the number of draws is high (~100 draws should complete in ~minutes)
	Requires:  /data/bactNOG/bactNOG.annotations.tsv
		   /results/metagenome_samples_per_NOG.tsv
		   /results/metagenome_coverage/asf_mg_coverage.tsv
		   /results/metagenome_coverage/bcto_mg_coverage.tsv
		   /results/metagenome_coverage/firm_mg_coverage.tsv
	Output:
		/results/asf_v_draws.svg
		   
