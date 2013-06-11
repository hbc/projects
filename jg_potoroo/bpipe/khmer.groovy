@Transform("keep")
khmer_digital_normalization = {
  exec """$khmer/normalized_by_median $optargs $input"""
  optargs=""
}

@Transform("kh")
khmer_count = {
  doc title: "Builds a counting hash using Khmer."
desc: """ This is a necessary requirement for many of the other steps. """
  
  exec """$KHMER/load-into-counting.py $optargs $output $input"""
  optargs=""
  forward input
}

@Transform("hist")
khmer_histogram = {
  doc title: "Creates an abundance histogram of kmers."
  exec """$KHMER/abundance-dist.py -z `basename $input .fasta`.kh $input $output"""
  forward input
}
