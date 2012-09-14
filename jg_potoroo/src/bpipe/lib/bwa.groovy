@Transform("bwa")
bwa_aln = {
  doc title: "Map reads using bwa."
  exec "bwa aln $optargs $genome $input > ${input}.bwa"
  exec """bwa samse $genome ${input}.bwa $input > $output"""
  forward input
}

@Transform("sam")
bwa_samse = {
  doc title: "Convert bwa to a sam file"
  exec "bwa samse $genome $input.bwa $input > $output"
}

bwa_index = {
  exec "bwa index -p $vectordir -a is $vector"
}  
