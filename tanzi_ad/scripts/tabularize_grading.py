import argparse
import yaml

TYPES = ['snp', 'indel', 'mixed']

def print_header():
    print 'sample', 'type', 'concordant', 'extra', 'missing', 'vardiff', 'hethom'

def print_grades(grading, var_type):
    concordant = grading['concordant']['concordant']
    discordant = grading['discordant'].get(var_type, {})
    shared = discordant.get('shared', {})
    extra = discordant.get('extra', {})
    missing = discordant.get('missing', {})
    print grading['sample'], var_type, concordant.get(var_type, 0), \
        extra.get('low-coverage', 0), missing.get('low-coverage', 0), \
        shared.get('vardiff', 0), shared.get('hethom', 0)



if __name__ == "__main__":
    description = "Flatten the grading files to something loadable by R."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('files', nargs='+', type=argparse.FileType('r'))
    args = parser.parse_args()


    print_header()
    for f in args.files:
        grading = yaml.load(f)[0]
        for t in TYPES:
            print_grades(grading, t)
