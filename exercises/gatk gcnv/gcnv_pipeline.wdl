workflow gcnv{
	call preprocess
}

task preprocess {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  File inputInterval
command { ${GATK} PreprocessIntervals \
        -R ${RefFasta} \
        -L ${inputInterval} \
		    --padding 0 \
        -imr OVERLAPPING_ONLY \
		    -O ${inputInterval}.interval_list
  }
  output {
    File interval_list = "${inputInterval}.interval_list"
  }
}