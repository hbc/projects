import scipy.stats as stats
import sys
import pybedtools
from itertools import izip

def main(control_bam, exp_bam, windows_bed):
    control_tool = pybedtools.BedTool(control_bam)
    exp_tool = pybedtools.BedTool(exp_bam)
    windows_tool = pybedtools.BedTool(windows_bed)
    control_coverage = control_tool.coverage(windows_tool)
    exp_coverage = exp_tool.coverage(windows_tool)

    norm_ratio = float(exp_tool.count()) / float(control_tool.count())
    for cwindow, ewindow in izip(control_coverage, exp_coverage):
        control_count = float(cwindow[8]) * norm_ratio
        exp_count = float(ewindow[8])
        if control_count == 0 and (exp_count > 2):
            print "\t".join(map(str, [ewindow[0], ewindow[1], ewindow[2],
                                      ewindow[3], control_count,
                                      exp_count, "0.0"]))
        elif (1 - stats.poisson.cdf(exp_count, control_count)) < 0.05:
            pval = (1 - stats.poisson.cdf(exp_count, control_count))
            print "\t".join(map(str, [ewindow[0], ewindow[1], ewindow[2],
                                      ewindow[3], control_count,
                                      exp_count, pval]))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
