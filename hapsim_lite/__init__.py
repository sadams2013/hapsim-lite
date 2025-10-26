"""
HapSim Lite
"""

__author__ = "Solomon M. Adams, PharmD, PhD"
__version__ = "0.0.1"

VCF_HEADER = """##fileformat=VCFv4.2
##source=HapSim_Lite
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
{contigs}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{samples}"""
