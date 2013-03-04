import argparse
import sh


def main(coordinate_file, vcf_file):
    tabix = _get_tabix_cmd(vcf_file)

    with open(coordinate_file) as coordinate_handle:
        for line in coordinate_handle:
            tabix_line = _lookup_line(tabix, line)
            if tabix_line:
                print tabix_line[0]

def _lookup_line(tabix, line):
    coordinate, rsid = _get_coordinate_and_rsid(line)
    tabix_output = tabix(coordinate)
    return _parse_tabix_output(tabix_output)
    #return coordinate, rsid, _parse_tabix_output(tabix_output)

def _parse_tabix_output(tabix_output):
    variants = tabix_output.strip().split("\n")
    valid_variants = filter(_is_valid_variant, variants)
    kept_variants = filter(_is_not_structural_variant, valid_variants)
    return kept_variants

def _is_valid_variant(variant):
    return len(variant.split("\t")) == 8

def _is_not_structural_variant(variant):
    return not "esv" in variant.split("\t")[2]

def _get_tabix_cmd(vcf_file):
    return sh.tabix.bake(vcf_file)

def _get_coordinate_and_rsid(line):
    coordinate = line.split("\t")[0]
    rsid = line.split("\t")[1]
    return _convert_coordinate(coordinate), rsid

def _convert_coordinate(coordinate):
    chromosome = coordinate.split(":")[0]
    base_string = coordinate.split(":")[1]
    base = int(base_string.split("-")[0]) + 1
    return "{chromosome}:{base}-{base}".format(chromosome=chromosome,
                                               base=base)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Looks up variants matching '
                                     'the coordinates in a file in a VCF file.')
    parser.add_argument('coordinate_file',
                        help="Name of file with coordinates to look up.")
    parser.add_argument('vcf_file',
                        help="Name of indexed file to look variants up in.")
    args = parser.parse_args()
    main(args.coordinate_file, args.vcf_file)
