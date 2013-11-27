PROJDIR="/n/scratch00/hsph/projects/js_trinity"
KHMER="~/opt/share/khmer/scripts"
KHMER_ARGS="-N 4 -x 4e9"

@Transform("normalized")
khmer_digital_normalization = {
  exec "$khmer/normalize-by-median.py $optargs $input"
}

@Transform("kh")
khmer_count = {
  doc title: "Builds a counting hash using Khmer."
desc: """ This is a necessary requirement for many of the other steps. """
author: "rory.kirchner@gmail.com"
  
  exec """$KHMER/load-into-counting.py $KHMER_ARGS $output $input"""
  forward input
}

@Transform("hist")
khmer_histogram = {
  doc title: "Creates an abundance histogram of kmers."
  exec """$KHMER/abundance-dist.py -z 
`basename $input .fasta`.kh $input $output
"""
  forward input
}

khmer_load_graph = {
  exec """$KHMER/load-graph.py $KHMER_ARGS -k 32 `basename $input .fasta` $input"""
  forward input
}

khmer_partition_graph = {
  exec """$KHMER/partition-graph.py --threads 4 -s `basename $input .fasta`"""
  exec """$KHMER/merge-partitions.py `basename $input .fasta`"""
  forward input
}

Bpipe.run {
  //khmer_count + khmer_histogram
  //khmer_load_graph + khmer_partition_graph
//khmer_digital_normalization.using(khmer:"~/opt/share/khmer/scripts", optargs:"-N 4 -x 4e9")
khmer_digital_normalization.using(khmer:"~/opt/share/khmer/scripts", optargs:"-N 4 -x 4e9 -C 20")
}
