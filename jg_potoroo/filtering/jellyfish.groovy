@Transform("count")
jellyfish_count = {
  exec """jellyfish count -m 31 -o $output -t 8 $input"""
}

Bpipe.run { jellyfish_count }
