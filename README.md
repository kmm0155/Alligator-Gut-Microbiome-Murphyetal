Title: Xenobiotic estradiol-17ß alters the microbial gut communities of hatchling American Alligators (Alligator mississippiensis)

Kaitlyn M. Murphy1, Madison M. Watkins1, John W. Finger1, Meghan D. Kelley1, Ruth M. Elsey2, Daniel A. Warner1, Mary T. Mendonça1

1Department of Biological Sciences, Auburn University, Auburn, AL 36849
2Louisiana Department of Wildlife and Fisheries, Grand Chenier, LA 

#This script (along with the attached text version) will give you the output for future analyses involved in this manuscript. For additional questions, feel free to email me at kmm0155@auburn.edu.

#Use this if running on your base computer
#!/bin/sh
#cd Desktop
#source activate qiime2-2020.2

#To demultiplex your samples, use the following code. You will need to download
#Dr. Demuxy from GitHub at the following address: https://github.com/lefeverde/Mr_Demuxy
#You need to make sure that Mr_Demuxy is in your path. If you installed Miniconda, for #QIIME, you #should already have pip. This is what you need to install Mr_Demux:
#pe_demuxer.py  -r1 R1.fastq -r2 R2.fastq -r2_bc reverse_bcs.txt -r1_bc forward_bcs.txt 

#I have paired end sequences and need to list the file as: #SampleData[PairedEndSequencesWithQuality]. R1 is forward and R2 is reverse.
#To see other types of 'type' files for QIIME tools import, see https://github.com/qiime2/q2cli/#issues/124
#The demultiplexer did not work, I asked for it to be done. In the future, always have your samples #demultiplexed.

#RUN THE CODE! Also, I usually run this in short step-by-step parts, not as an entire code. I used the Alabama Super Computer and ran "conda activate qiime2-dev"
		
#Create qiime aritfact
#https://docs.qiime2.org/2019.7/tutorials/importing/ 
#You should already have your manifest made before reaching this next step.

qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Manifest.csv \
	--output-path paired-end-seqs.qza \
	--input-format PairedEndFastqManifestPhred33
	 
#Use this to peek at the data
qiime tools peek paired-end-seqs.qza

#Denoises sequences, dereplicates them, and filters chimeras
#https://docs.qiime2.org/2018.6/plugins/available/dada2/denoise-pyro/
#You should have plugged in a few files into usegalaxy.org and figured out where to truncate your data before running the next line of code. 

qiime dada2 denoise-paired \
	--i-demultiplexed-seqs paired-end-seqs.qza \
	--p-trim-left-f 20 \
	--p-trim-left-r 20 \
	--p-trunc-len-f 235 \
	--p-trunc-len-r 235 \
	--p-trunc-q 15 \
	--p-chimera-method consensus \
	--o-representative-sequences rep-seqs-denoise.qza \
	--o-table rep_seq_feature_table.qza \
	--o-denoising-stats denoising-stats.gza \
	--verbose
	
#Summary stats of denoise and quality filtering

qiime feature-table summarize \
	--i-table rep_seq_feature_table.qza \
	--o-visualization rep_seq_feature_table-view.qzv \
	--m-sample-metadata-file /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata.txt

#Feature table!	

qiime feature-table tabulate-seqs \
	--i-data rep-seqs-denoise.qza \
	--o-visualization rep-seqs-view.qzv

#Filter features from feature table
#features must be a minimum sum of 8 across all samples and must be present in at least 4 samples: https://docs.qiime2.org/2019.7/tutorials/filtering/
	
qiime feature-table filter-features \
	--i-table rep_seq_feature_table.qza \
	--p-min-frequency 8 \
	--p-min-samples 4 \
	--o-filtered-table rep_seq_feature_table2.qza
	
#Now filter sequences to match table 
#https://docs.qiime2.org/2018.8/plugins/available/feature-table/filter-seqs/

qiime feature-table filter-seqs \
	--i-data rep-seqs-denoise.qza \
	--i-table rep_seq_feature_table2.qza \
	--o-filtered-data rep-seqs-filtered.qza
	
#Summary stats of filtering

qiime feature-table summarize \
	--i-table rep_seq_feature_table2.qza \
	--o-visualization rep_seq_feature_table_filter-view.qzv \
	--m-sample-metadata-file /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata.txt

#Rarefaction curve

qiime diversity alpha-rarefaction \
	--i-table rep_seq_feature_table2.qza \
	--p-max-depth 197 \
	--m-metadata-file /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata.txt \
	--o-visualization alpha-rarefaction.qzv

#Taxonomy Classification and taxonomic analyses
#https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
#https://docs.qiime2.org/2018.6/tutorials/feature-classifier/
#https://forum.qiime2.org/t/silva-132-classifiers/3698
#https://www.dropbox.com/s/5tv5uk95pk3ukwf/7_level_taxonomy.qza?dl=0

#Consensus must match 100% where majority is 90% to the cluster; 7 levels refers to taxonomic levels

qiime feature-classifier extract-reads \
	--i-sequences /mnt/beegfs/home/aubkmm/99_otus.qza \
	--p-f-primer GTGCCAGCMGCCGCGGTAA \
	--p-r-primer GGACTACHVGGGTWTCTAAT \
	--o-reads gg-ref-seqs.qza
	
#training sklearn

qiime feature-classifier fit-classifier-naive-bayes \
	--i-reference-reads gg-ref-seqs.qza \
	--i-reference-taxonomy /mnt/beegfs/home/aubkmm/99_otu_taxonomy.qza \
	--o-classifier gg-99-classifier.qza
  
qiime feature-classifier classify-sklearn \
	--i-classifier gg-99-classifier.qza \
	--i-reads rep-seqs-filtered.qza \
	--o-classification classified_taxonomy_table.qza
	
qiime metadata tabulate \
	--m-input-file classified_taxonomy_table.qza \
	--o-visualization classified_taxonomy.qzv

#Plug the below barplot into view.qiime2 to see what it looks like!

qiime taxa barplot \
	--i-table rep_seq_feature_table2.qza \
	--i-taxonomy classified_taxonomy_table.qza \
	--m-metadata-file /Users/kaitlynmurphy/Desktop/Alligator_Data/Alligator_Metadata.txt \
	--o-visualization taxa-barplots.qzv
	
#Qiime taxa collapse 

qiime taxa collapse \
	--i-table rep_seq_feature_table2.qza \
	--i-taxonomy classified_taxonomy_table.qza \
	--p-level 7 \
	--output-dir rep_seq_feature_table3_collapsed
	
#Export data

#exports the feature-table.biom
qiime tools export \
	--input-path /mnt/beegfs/home/aubkmm/rep_seq_feature_table3_collapsed/collapsed_table.qza \
	--output-path exported

#exports the taxonomy.tsv
qiime tools export \
	--input-path classified_taxonomy_table.qza \
	--output-path exported

#Go back into the MetaData file and add a '#' to the title line. This will allow biom to detect the header
	
biom add-metadata \
	-i /mnt/beegfs/home/aubkmm/exported/feature-table.biom \
	-o Alligator_w_taxonomy.biom \
	--observation-metadata-fp /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata1.txt \
	--sc-separated taxonomy

biom convert -i Alligator_w_taxonomy.biom -o Alligator_w_taxonomy.tsv --to-tsv

#Filter/separate treatment groups based on treatment in order to determine differential abundance
#Remove '#' from metadata file

qiime feature-table filter-samples \
	--i-table /mnt/beegfs/home/aubkmm/rep_seq_feature_table3_collapsed/collapsed_table.qza \
	--m-metadata-file /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata.txt \
	--p-where "subject='Control'" \
	--o-filtered-table Control-abundance-table.qza

#Save as a biom file
qiime tools export \
  --input-path Control-abundance-table.qza \
  --output-path Control

#Convert to tsv file
biom convert -i /mnt/beegfs/home/aubkmm/Control/feature-table.biom -o Alligator_Control.tsv --to-tsv

#Do the sam for Low treatment groups
qiime feature-table filter-samples \
	--i-table /mnt/beegfs/home/aubkmm/rep_seq_feature_table3_collapsed/collapsed_table.qza \
	--m-metadata-file /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata.txt \
	--p-where "subject='Low'" \
	--o-filtered-table Low-abundance-table.qza

#Save as a biom file
qiime tools export \
  --input-path Low-abundance-table.qza \
  --output-path Low

#Convert to a tsv file
biom convert -i /mnt/beegfs/home/aubkmm/Low/feature-table.biom -o Alligator_Low.tsv --to-tsv

#Lastly, do the same for the high treatment group
qiime feature-table filter-samples \
	--i-table /mnt/beegfs/home/aubkmm/rep_seq_feature_table3_collapsed/collapsed_table.qza \
	--m-metadata-file /mnt/beegfs/home/aubkmm/Alligator_Data_2/Alligator_Metadata.txt \
	--p-where "subject='High'" \
	--o-filtered-table High-abundance-table.qza

#Save as a tsv file
qiime tools export \
  --input-path High-abundance-table.qza \
  --output-path High

#Take out of 'high' file and put onto Desktop
biom convert -i /mnt/beegfs/home/aubkmm/High/feature-table.biom -o Alligator_High.tsv --to-tsv

##DIVERSITY##
#https://docs.qiime2.org/2018.11/plugins/available/diversity/

qiime diversity alpha \
  --i-table rep_seq_feature_table2.qza \
  --p-metric shannon \
  --o-alpha-diversity shannon_vector.qza

#Exports file in tsv format
qiime tools export \
  --input-path shannon_vector.qza \
  --output-path shannons
    	
