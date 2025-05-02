import argparse, sys, os, pysam
from typing import List


def create_parser() -> argparse.ArgumentParser:
    # :return: creates parser with all command line arguments arguments
    MSMuTect_intro = "MSMuTect\n Version 4.0\n Authors: Yossi Maruvka, Avraham Kahan, and the Maruvka Lab at Technion"
    parser = argparse.ArgumentParser(description=MSMuTect_intro)
    parser.add_argument("-T", "--tumor_file", help="Tumor BAM file")
    parser.add_argument("-N", "--normal_file", help="Non-tumor BAM file")
    parser.add_argument("-S", "--single_file", help="Analyze a single file for histogram and/or alleles")
    parser.add_argument("-l", "--loci_file", help="File of loci to be processed and included in the output", required=True)
    parser.add_argument("-O", "--output_prefix", help="prefix for all output files", required=True)
    parser.add_argument("-c", "--cores", help="Number of cores to run MSMuTect on", type=int, default=1)
    parser.add_argument("-b", "--batch_start", help="1-indexed number locus to begin analyzing at (Inclusive)", default=1, type=int)
    parser.add_argument("-e", "--batch_end", help="1-indexed number locus to stop analyzing at (Inclusive)", type=int)
    parser.add_argument("-H", "--histogram", help="Output a Histogram File", action='store_true')
    parser.add_argument("-A", "--allele", help="Output allele file", action='store_true')
    parser.add_argument("-m", "--mutation", help="Output mutation file", action='store_true')
    parser.add_argument("-F", "--flanking", help="Length of flanking on both sides of an accepted read", type=int, default=10)
    parser.add_argument("-r", "--read_level", help="Minimum number of reads to call allele", type=int, default=5)
    parser.add_argument("-f", "--force", help="Overwrite pre-existing files", action='store_true')
    parser.add_argument("--integer", help="Only use indels of integer deletions/insertions of the repeat unit when calling alleles", action='store_true')
    parser.add_argument("--vcf", help="Output VCF file in addition to TSV file. IMPORTANT NOTE: the location of the indel may be anywhere in the locus, but the VCF will always have it at the beginning of the locus. In addition, all loci with non-reference alleles are listed but only mutations are listed as 'PASS'", action='store_true')
    return parser


def exit_on(message: str, status: int = 1):
    # print message, and exit
    print("ERROR: " + message)
    sys.exit(status)


def simple_index_check(bam: str):
    bam_bai_path = bam + ".bai"
    bai_path = bam[:-4] + ".bai"
    bam_bai_file_exists = os.path.exists(bam_bai_path)
    bai_file_exists = os.path.exists(bai_path)
    index_file_older_than_bam_message = "Index file older than BAM file. Index file must be younger than BAM file. If you are sure the index file is correct, run 'touch [index_file]'"
    if not bam_bai_file_exists and not bai_file_exists:
        exit_on("Given BAM file/s are not sorted and/or indexed")


def validate_indexing(bam_files: List[str]) -> None:
    """ validates that given BAM files are indexed"""
    prefixes = ['', 'chr', 'Chr']
    for bam in bam_files:
        simple_index_check(bam) # simply checks for .bai file
        validated = False
        current_handle = pysam.AlignmentFile(open(bam, 'rb'))
        for prefix in prefixes:
            try:
                _ = current_handle.fetch(f"{prefix}{1}", start=10_000, multiple_iterators=False)
                validated = True
                break  # verified
            except ValueError:  # different prefix, or not indexed
                continue
        if not validated:
            exit_on("Given BAM file/s are not sorted and/or indexed, or contain an unsusual prefix (not 'chr', 'Chr', or nothing")


def validate_bams(arguments: argparse.Namespace):
    if (bool(arguments.tumor_file) or bool(arguments.normal_file)) == bool(arguments.single_file): #  XOR
        exit_on("Provide Single file, or both Normal and Tumor file")
    elif bool(arguments.tumor_file) != bool(arguments.normal_file): #  XOR
        exit_on("Provide Single file, or both Normal and Tumor file")
    elif arguments.single_file:
        if not os.path.exists(arguments.single_file):
            exit_on("Provided single file path does not exist")
    else:
        if not os.path.exists(arguments.tumor_file) or not os.path.exists(arguments.normal_file):
            exit_on("Provided Normal or Tumor BAM path does not exist")
    validate_indexing([bam_file for bam_file in [arguments.tumor_file, arguments.normal_file, arguments.single_file] if bool(bam_file)])


def validate_output_files(arguments: argparse.Namespace):
    overwrite_files_mssg = "Files would be overwritten by this run. To force overwrite, use -f flag"
    if os.path.sep not in arguments.output_prefix:
        arguments.output_prefix = os.path.join(os.getcwd(), arguments.output_prefix)
    if not os.path.exists(os.path.dirname(arguments.output_prefix)):
        exit_on("Output directory does not exist")
    if arguments.force:
        return
    elif arguments.single_file:
        if arguments.histogram and not arguments.allele:
            if os.path.exists(arguments.output_prefix + ".hist.tsv"):
                exit_on(overwrite_files_mssg)
        else:
            if os.path.exists(arguments.output_prefix + ".all.tsv"):
                exit_on(overwrite_files_mssg)
    else:  # pair file
        if arguments.mutation and arguments.vcf and os.path.exists(arguments.output_prefix+".vcf"):
            exit_on(overwrite_files_mssg)
        if (arguments.histogram or arguments.allele) and arguments.mutation and os.path.exists(arguments.output_prefix + ".full.mut.tsv"):
                    exit_on(overwrite_files_mssg)
        elif not arguments.histogram and not arguments.allele and os.path.exists(arguments.output_prefix + ".partial.mut.tsv"):
                exit_on(overwrite_files_mssg)
        elif arguments.allele and (os.path.exists(arguments.output_prefix + ".tumor.all.tsv") or
        os.path.exists(arguments.output_prefix + ".normal.all.tsv")):
            exit_on(overwrite_files_mssg)
        elif arguments.histogram and (os.path.exists(arguments.output_prefix + ".tumor.hist.tsv") or
                                      os.path.exists(arguments.output_prefix + ".normal.hist.tsv")):
            exit_on(overwrite_files_mssg)


def validate_input(arguments: argparse.Namespace):
    validate_bams(arguments)
    validate_output_files(arguments)
    if not os.path.exists(arguments.loci_file):
        exit_on("Provided loci file does not exist")
    if arguments.mutation and arguments.single_file:
        exit_on("Pair of files must be provided to call mutations")
    elif arguments.batch_start <= 0:
        exit_on("Batch Start must be equal to or greater than 1")
    elif arguments.cores <= 0:
        exit_on("Cores must be equal to or greater than 1")
    elif arguments.flanking < 0:
        exit_on("Flanking must be equal to or greater than 0")
    elif arguments.read_level < 1:
        exit_on("Minimum Read Level for calling alleles must be equal to or greater than 1")
    elif not os.path.exists(arguments.loci_file):
        exit_on("Loci file path does not exist")
    elif arguments.vcf and not arguments.mutation:
        exit_on("VCF file can only be generated for mutation calls")

