@Transform("count")
jellyfish_count = {
  exec """jellyfish count -m 31 -o $projectdir/data/jellyfish/$output -t 8 -s 10G $input"""
}

jellyfish_merge = {
  exec """jellyfish  """
}

Bpipe.run { jellyfish_count }
