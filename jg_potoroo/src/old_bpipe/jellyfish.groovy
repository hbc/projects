@Transform("jellycount")
jellyfish_count = {
  exec "jellyfish count $userargs -o $output $input"
  exec "jellyfish merge -o $output ${output}'/_*'"
}

@Transform("jellydump")
jellyfish_dump = {
  exec """jellyfish dump -o $output -L 1000000 $input"""
}

@Transform("jellyhist")
jellyfish_hist = {
  exec """jellyfish histo -o $output $input"""
}

@Transform("jellystats")
jellyfish_stats = {
  exec """jellyfish stats -o $output $input"""
}

Bpipe.run{ jellyfish_count.using(userargs: "-m 31 -t 8 -s 10G") }
  
//Bpipe.run { jellyfish_count }
//Bpipe.run { jellyfish_dump }
//Bpipe.run { jellyfish_hist }
//Bpipe.run { jellyfish_stats }
