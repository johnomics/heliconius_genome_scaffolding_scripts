These scripts were written to build a linkage map of Heliconius melpomene and use the linkage map to produce a chromosomal scaffold of the Heliconius melpomene genome.

The scripts are provided FOR INFORMATION ONLY and are not intended for use. They contain many Heliconius-specific details and have not been adapted for general use. Input and output files for each script are given in the comments for each file.

Please see paper for full details of results and methods, particularly Supplementary Information Section S4.

Publication: Heliconius Genome Consortium, "Butterfly genome reveals promiscuous exchange of mimicry adaptations among species", Nature, DOI: 10.1038/nature11041.

Author: John Davey, john.davey@ed.ac.uk, July 2011 - February 2012.



Pipelines:

Generate chromosomal AGP file and revise scaffold AGP file based on linkage:
1. Generate VCF file for RAD-sequenced cross
2. convert_vcf_to_joinmap_heliconius.pl
   Generate JoinMap-format markers and positions of those markers on scaffolds
3. validate_markers_for_mapping.pl
   Error-correct and convert intercross markers for mapping
4. Generate linkage maps using JoinMap or similar software
5. scaffold_heliconius_genome.pl
   Synthesise linkage maps and scaffolds to produce chromosome AGP file

Revise chromosomal AGP file based on mate pair information, using Bambus2:
1. generate_bambus2_agps.pl
   Run Bambus2 on matepair data to find additional scaffold bridges
2. revise_chromosome_agp.pl
   Use Bambus2 output to revise chromosome AGP

Generate scaffolding statistics and figure S4.6.1:
1. generate_genome_stats.pl
2. data_for_rad_scaffold_figure.pl
3. draw_rad_scaffold_figure.R