@Filter("sickled")
sickle = {
  doc title: "Trim ends of reads based on quality using sickle."
  exec "sickle $pair -f $input -o $output -t $qtype $optargs"
  pair=""
  qtype=""
  optargs=""
}
