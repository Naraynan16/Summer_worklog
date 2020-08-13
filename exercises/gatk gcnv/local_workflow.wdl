workflow gcnv{
	call preprocess
}

task preprocess {
command {  
    /mnt/c/Users/Narayanan/Downloads/Thesis/gatk-4.1.8.1/gatk PreprocessIntervals \
    -R /mnt/c/Users/Narayanan/Downloads/Thesis/ref/Homo_sapiens_assembly19.fasta \
        -L /mnt/c/Users/Narayanan/Downloads/Thesis/SSV4.bed \
		--padding 0 \
        -imr OVERLAPPING_ONLY \
		-O Preprocessed_SSV4.interval_list
  }
  output {
    File interval_list = "Preprocessed_SSV4.interval_list"
  }
}

