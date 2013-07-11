import sys
import pysam

#chr10	102016028	102016066	VXHF8:1806:1874	 gene_id "ENSG00000095485"

def main(tmp_file, bam_file):
    print "#gffTags"
    with open(tmp_file) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            lsplit = line.split("\t")
            read_name = lsplit[3]
            read = find_read(bam_file, read_name)
            if read:
                lsplit[3] = "read_id=" + ";".join([read_name, "=".join(lsplit[4].split())])
                lsplit[3] += ';sequence="%s"' % (read.seq)
                del lsplit[4]
                print "\t".join(lsplit)

def find_read(bam_file, read_name):
    with pysam.Samfile(bam_file, "rb") as in_handle:
        for read in in_handle:
            if read.qname == read_name:
                return read
        return None


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
