import pysam
import sys
import pandas as pd
import rpy2.robjects as robjects
from rkinf.utils import replace_suffix
import pandas.rpy.common as com


if __name__ == "__main__":
    filename = sys.argv[1]
    samfile = pysam.Samfile(filename, "rb")
    multi_counts = {}
    read_lengths = {}
    for read in samfile:
        multi_counts[read.qname] = multi_counts.get(read.qname, 0) + 1
        read_lengths[read.qname] = len(read.seq)

    df = pd.DataFrame([multi_counts, read_lengths], index=["count", "length"])
    df = df.T

    r = robjects.r
    r.assign('plot.file', replace_suffix(filename, ".pdf"))
    r.assign('r_df', com.convert_to_r_dataframe(df))

    r('''
    library(ggplot2)
    p = ggplot(r_df, aes(count, length)) + geom_point()
    ggsave(filename=plot.file)
    ''')
