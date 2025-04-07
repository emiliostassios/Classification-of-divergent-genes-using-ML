HOW TO RUN THE SCRIPTS AND A SHORT DESCRIPTION OF WHAT EVERYTHING DOES.  

You should run these scripts with the conda environment, that is provided, enabled.

The diamond_feature_extractor.py is the script used to extract the features from the DIAMOND output.  
    The way to run it is: python diamond_feature_extractor.py $DIAMOND_data_frame $simulated_sequences.fasta $simulated_orphans.fasta $directory_to_save_the_feature_data_frame
    The output of this script is saved to a directory called 'diamond_features' and called 'diamond_features.tsv'. This file contains rows with NAs that should be cleaned. Two more scripts should be run.
    The first one is normalization.py that normalizes the alignment_count feature with the length of the query sequence. The way to run it is:
          python normalization.py $DIAMOND_FEATURE_DATA_FRAME $path/to/simulated_fasta $output_file
    The second script, class.py, is for creating the CLASS column. The way to run it is:
          python class.py $DIAMOND_FEATURE_DATA_FRAME $class (neg or pos)
  
