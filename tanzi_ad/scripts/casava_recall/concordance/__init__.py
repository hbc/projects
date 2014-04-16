import vcf
import sys
from itertools import izip_longest, islice

class Variant(object):
    def __init__(self, chrom, start, end, ref, alt, samples):
        self.chrom = str(chrom)
        self.start = int(start)
        self.end = int(end)
        self.ref = str(ref)
        self.alt = str(alt)
        self.samples = set(samples)
    def __eq__(self, other):
        return (self.chrom == other.chrom and self.start == other.start and
                self.end == other.end and self.ref == other.ref and
                self.alt == other.alt)
#                self.alt.intersection(other.alt))
    def __repr__(self):
        return "\t".join(map(str, [self.__class__, self.chrom, self.start, self.end,
                          self.ref, self.alt]))
    def __add__(self, other):
        if other:
            self.samples = self.samples.union(other.samples)


def _get_vcf_handles(vcf_files):
    return [vcf.Reader(filename=v) for v in vcf_files]

def _lookup_variant(handle, chrom, start, end):
    chrom = chrom.split("chr")[1]
    try:
        variants = handle.fetch(chrom, int(start), int(end))
    except ValueError:
        sys.stderr.write("\t".join(map(str[chrom, start, end])))
        variants = None
    return variants

def _concordance(gemini, recalled):
    if not recalled:
        concordant = set([])
        recalled_only = set([])
        gemini_only = gemini.samples
        mismatch = set([])

    elif gemini == recalled:
        concordant = gemini.samples.intersection(recalled.samples)
        recalled_only = recalled.samples.difference(gemini.samples)
        gemini_only = gemini.samples.difference(recalled.samples)
        mismatch = set([])
    else:
        concordant = set([])
        mismatch = gemini.samples.intersection(recalled.samples)
        recalled_only = recalled.samples.difference(gemini.samples)
        gemini_only = gemini.samples.difference(recalled.samples)
    return {"concordant": "|".join(concordant),
            "recalled": "|".join(recalled_only),
            "gemini": "|".join(gemini_only),
            "mismatch": "|".join(mismatch)}

def match_gemini_variant(gemini_variant, variants):
    samples = []
    for variant in variants:
        alleles = map(str, variant.alleles)
        if gemini_variant.alt not in alleles:
            continue
        match = str(alleles.index(gemini_variant.alt))
#        import ipdb
#        ipdb.set_trace()
        for v in variant.samples:
            try:
                v.gt_alleles
            except:
                continue
            if match in v.gt_alleles:
                samples.append(v.sample)
    return samples


def concordance_calculator(line, vcf_handles, to_consider):
    # skip variants with more than one allele
    if len(line['alt'].split(",")) > 1:
        return ""
    gemini_samples = [x for x in line['variant_samples'].split(",") if x in to_consider]
    gemini_variant = Variant(line['chrom'].split("chr")[1],
                             line['start'], line['end'],
                             line['ref'], line['alt'],
                             gemini_samples)
    same_variant = None
    mismatch_variants = []
    matches = []
    for handle in vcf_handles:
        variants = _lookup_variant(handle, line['chrom'], line['start'], line['end'])
        matches += match_gemini_variant(gemini_variant, variants)
        # for variant in variants:
        #     if variant.FILTER:
        #         continue
        #     samples = [sample.sample for sample in variant.samples
        #                if sample.is_variant and sample.sample in to_consider]
        #     new_variant = Variant(variant.CHROM, variant.POS-1, variant.end,
        #                           variant.REF, "|".join(map(str, variant.ALT)),
        #                           samples=samples)
        #     import ipdb
        #     ipdb.set_trace()
        #     if new_variant.ref != gemini_variant.ref:
        #         continue
        #     if new_variant == gemini_variant:
        #         if not same_variant:
        #             same_variant = new_variant
        #         else:
        #             same_variant += new_variant
        #     else:
        #         mismatch_variants.append(new_variant)
    matches = set(to_consider).intersection(set(matches))
    gemini_variant.samples = gemini_variant.samples.intersection(set(to_consider))
    concordant = gemini_variant.samples.intersection(matches)
    num_concordant = len(concordant)
    recalled_only = matches.difference(gemini_variant.samples)
    num_recalled_only = len(recalled_only)
    casava_only = gemini_variant.samples.difference(matches)
    num_casava_only = len(casava_only)


#    concordance = _concordance(gemini_variant, same_variant)
#    mismatch = "|".join(["|".join(v.samples) for v in mismatch_variants])
#    mismatch_alleles = "|".join([v.ref + "/" + str(list(v.alt))  for v in mismatch_variants])
#    concordant = concordance["concordant"].split("|")
#    num_concordant = len(concordance["concordant"].split("|")) if concordance["concordant"] else 0
#    num_gemini = len(concordance["gemini"].split("|")) if concordance["gemini"] else 0
#    num_recalled = len(concordance["recalled"].split("|")) if concordance["recalled"] else 0
#    num_mismatch = len(mismatch.split("|")) if mismatch else 0
    return ",".join(map(str, [gemini_variant.chrom, gemini_variant.start,
                              gemini_variant.end, gemini_variant.ref,
                              gemini_variant.alt,
                              num_concordant, num_casava_only,
                              num_recalled_only]))


                                 # "|".join(gemini_variant.alt), concordance["concordant"],
                                 # concordance["gemini"], concordance["recalled"],
                                 # mismatch, mismatch_alleles,
                                 # str(len(line['variant_samples'].split(","))),
                                 # str(num_concordant),
                                 # str(num_gemini),
                                 # str(num_recalled),
                                 # str(num_mismatch)]))

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))
