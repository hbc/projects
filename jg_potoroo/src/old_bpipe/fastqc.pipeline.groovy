PROJDIR="/n/scratch00/hsph/projects/js_trinity"

@Transform("fastqc")
fastqc = {
  exec "mkdir $output"
  exec "fastqc -o $output -t 6 $input"
}

Bpipe.run {
  "530_%_s" * [fastqc]
}
