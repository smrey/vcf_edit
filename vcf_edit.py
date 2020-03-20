import vcf
import argparse
import textwrap


def get_args():
    argument_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent(
            '''
            summary:
            Describe software!!!
            '''))

    # Add arguments
    # REQUIRED: Path to config file containing authorisation information
    argument_parser.add_argument(
        '-v', '--path_to_vcf', action='store',
        help=textwrap.dedent(
            '''
            File path to vcf file to be edited. REQUIRED.
            '''))

    # OPTIONAL:
    optional_options = argument_parser.add_mutually_exclusive_group()
    return argument_parser.parse_args()


def load_vcf(vcf_file_path):
    print(f"vcf loaded")
    read_vcf = vcf.Reader(open(vcf_file_path, 'r'))
    vcf_headers = read_vcf.metadata
    print(read_vcf.infos)
    print(read_vcf.filters)
    print(read_vcf.formats)
    print(read_vcf.samples)
    for call in read_vcf:
        # To read out
        print(call)
        print(call.ALT)
        print(call.CHROM)
        print(call.POS)
        print(call.ID)
        print(call.REF)
        print(call.QUAL)
        print(call.FILTER)
        print(call.INFO)
        if len(call.ALT) != 1:
            raise Exception(f"Only SNPs supported")


def main():
    __version__ = '0.0.1'
    __updated__ = '20/03/2020'
    args = get_args()
    load_vcf(args.path_to_vcf)


if __name__ == '__main__':
    main()