from argparse import ArgumentParser
import glob
import os
import subprocess
import time

BCBIO_VARIATION = ("/n/hsphS10/hsphfs1/chb/local/share/java/bcbio_variation/"
                   "bcbio.variation-0.1.6-SNAPSHOT-standalone.jar")

def get_family_lookup(famfile):
    with open(famfile) as in_handle:
        family_lookup = {line.split("\t")[1]: line.split("\t")[0]
                         for line in in_handle}
    return family_lookup

def get_vcfs(vcf_dir):
    vcfs = list(glob.glob(os.path.join(vcf_dir, "*.vcf.gz")))
    return [vcf for vcf in vcfs if "passonly" not in vcf]

def sample_from_vcf(vcf):
    return os.path.basename(vcf).split(".vcf.gz")[0]

def family_from_vcf(vcf):
    return os.path.basename(vcf).split("-")[0]

def get_casava_samples(casava):
    casava_vcf_files = get_vcfs(casava)
    return [sample_from_vcf(x) for x in casava_vcf_files]

def get_recalled_families(recalled):
    recalled_vcf_files = get_vcfs(recalled)
    return [os.path.basename(x).split("-")[0] for x in recalled_vcf_files]

def get_samples_to_check(casava, recalled, family_lookup):
    samples = get_casava_samples(casava)
    families = get_recalled_families(recalled)
    to_check = [x for x in samples if x in family_lookup
                and family_lookup[x] in families]
    done = [x.split("-")[0] for x in glob.glob("*-discordant.vcf")]
    to_check = [x for x in to_check if x not in done]
    return to_check


def get_family_vcf_lookup(vcfs):
    return {family_from_vcf(vcf): vcf for vcf in vcfs}

def run(cmd):
    import subprocess
    return subprocess.check_call(cmd, shell=True)


if __name__ == "__main__":
    description = ("Perform concordance checks between two sets of variants.")
    parser = ArgumentParser(description=description)
    parser.add_argument("--bed", help="A list of regions to check.")
    parser.add_argument("--ref", help="Path to reference sequence.")
    parser.add_argument("casava", help="A directory of CASAVA calls.")
    parser.add_argument("recalled", help="A directory of recalled family-wise calls.")
    parser.add_argument("fam", help="FAM file.")
    args = parser.parse_args()
    family_lookup = get_family_lookup(args.fam)
    to_check = get_samples_to_check(args.casava, args.recalled, family_lookup)

    casava_vcf_files = get_vcfs(args.casava)
    family_vcf_lookup = get_family_vcf_lookup(get_vcfs(args.recalled))
    cmds = []

    for casava_vcf in casava_vcf_files:
        sample = sample_from_vcf(casava_vcf)
        if sample not in family_lookup:
            print "Skipping %s because it is not in the family table." % sample
            continue
        family = family_lookup[sample]
        if sample not in to_check:
            continue
        else:
            recalled_vcf = family_vcf_lookup[family]
            print "Processing %s and %s." % (casava_vcf, recalled_vcf)
            cmd = ("sbatch -p general --mem=10000 --wrap='java -jar -Xms750m -Xmx8000m {BCBIO_VARIATION} variant-utils comparetwo "
                   "{casava_vcf} {recalled_vcf} {args.ref} {args.bed}'").format(**locals())
            print "Running %s." % cmd
            cmds.append(cmd)
    for cmd in cmds:
        subprocess.check_call(cmd, shell=True)
        time.sleep(30)
