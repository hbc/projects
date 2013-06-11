@Transform("fastqc")
fastqc = {
  exec "mkdir $output"
  exec "fastqc -o $output $optargs $input"
}
