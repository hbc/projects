Filter("dusted")
tagdust = {
  doc title: "Filter out possible contaminants with Tagdust.",
desc: """ Filters out possible contaminants listed in the fasta
file pointed to by "$contam"""",
constraints: """ example usage: tagdust.using(contam="contamination_file") """,
author: "rory.kirchner@gmail.com"
  
  exec "tagdust -s $optargs -o ${output} $contam ${input}"
}
