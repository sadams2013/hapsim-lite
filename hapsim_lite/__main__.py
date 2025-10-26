"""
Entrypoint for hapsim-lite CLI
"""

import argparse
import logging
import sys

from .population import PopulationData
from .generate import HaplotypeGenerator
from .vcf import write_vcf, generate_sample_ids, generate_call_matrix


def main():
    """CLI Entrypoint"""
    parser = argparse.ArgumentParser(
        prog="HapSim Lite",
        description="Simulate phased variant calls from pre-calculated frequencies and LD",
    )
    parser.add_argument(
        "-f", "--maf-file", type=str, help="Plink2 afreq file with allele frequences"
    )
    parser.add_argument(
        "-l",
        "--ld-file",
        type=str,
        help="Plink2 vcor file with unphased R values (not R^2)",
    )
    parser.add_argument(
        "-n", "--num-samples", type=int, help="Number of samples to output"
    )
    parser.add_argument(
        "-p", "--ploidy", type=int, help="Ploidy in output VCF", default=2
    )
    parser.add_argument(
        "--maf-only",
        action="store_true",
        help="Skip Markov-Chain and just simulate based on allele frequencies",
    )
    # output options
    parser.add_argument(
        "--unphased", action="store_true", help="output unphased results"
    )
    # hyperparameters
    parser.add_argument(
        "-w", "--ld-window", type=int, help="LD window in adjacent variants", default=2
    )
    parser.add_argument(
        "-t", type=float, help="LD smoothing hyperparameter", default=0.3
    )
    parser.add_argument(
        "-d",
        type=float,
        help="Distance decay for LD - meaured in kb. \
            Higher numbers allows for higher contributions of long-distance ld",
        default=7.5,
    )
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stderr, level=logging.INFO)

    logging.info("Parsing population data inputs. ")
    pop_data = PopulationData.from_plink2_afreq_vcor(  # pylint: disable=invalid-name
        args.maf_file, args.ld_file
    )
    logging.info("Generating Haplotypes")
    hap = HaplotypeGenerator(
        pop_data,
        args.num_samples * args.ploidy,
        window=args.ld_window,
        tau=args.t,
        lam=args.d,
    )
    if args.maf_only:
        logging.info(
            "Skipping Markov-Chain walk, simulation will not have realistic LD patterns"
        )
    else:
        logging.info("Generating initial forward pass")
        hap.forward_pass()
        logging.info("Generating final reverse pass")
        hap.reverse_pass()
    logging.info("Writing VCF")
    call_matrix = generate_call_matrix(
        hap.hap_matrix, args.ploidy, args.num_samples, args.unphased
    )
    sample_ids = generate_sample_ids(args.num_samples)
    write_vcf(call_matrix, pop_data, sample_ids)
