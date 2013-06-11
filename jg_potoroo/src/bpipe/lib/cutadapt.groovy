@Filter("nopolya") 
trim_polya = {
  doc title: "Trims off polyA from 3' ends of reads using cutadapt"
  exec "cutadapt $optargs -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAA $input > $output"
  optargs=""
}

@Filter("cut")
cutadapt = {
  doc title: "Generic cutadapt to remove arbitrary contaminants,"
  exec """cutadapt $optargs $adaptor $input > $output"""
  optargs=""
  adaptor=""
}

@Filter("noillum")
trim_illum = {
  doc title: "Trims off Illumina adaptors from ends of reads (NOT IMPLEMENTED)"
  exec "echo NOT IMPLEMENTED YET"
}
