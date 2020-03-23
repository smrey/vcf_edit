import vcf
import argparse
import textwrap

# Dictionary of IUPAC codes
iupac_dict = {('A', 'G'): 'R', ('C', 'T'): 'Y', ('C', 'G'): 'S', ('A', 'T'): 'W', ('G', 'T'): 'K', ('A', 'C'): 'M',
              ('C', 'G', 'T'): 'B', ('A', 'G', 'T'): 'D', ('A', 'C', 'T'): 'H', ('A', 'C', 'G'): 'V',
              ('A', 'C', 'G', 'T'): 'N'}

# This will use a lot of memory for large vcfs


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

   # OPTIONAL: TODO Support these options
    optional_options = argument_parser.add_mutually_exclusive_group()
    optional_options.add_argument(
        '-o', '--output_file', action='store_true', default=False,
        help=textwrap.dedent(
            '''
            Add output filename. Default is {input filename}_edited.vcf. OPTIONAL.
            '''))
    optional_options.add_argument(
        '-b', '--bcftools', action='store_true', default=False,
        help=textwrap.dedent(
            '''
           Specify that bcftools variant caller that has been used to generate the VCF. OPTIONAL.
            '''))

    return argument_parser.parse_args()


def write_vcf(records_to_write, read_vcf, output_vcf="test.vcf"):
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), read_vcf)
    for record in records_to_write:
        vcf_writer.write_record(record)


def single_snp_call(ref_snp, called_snp):
    # Add appropriate IUPAC code
    find_iupac = [ref_snp, str(called_snp[0])]
    find_iupac.sort()
    find_iupac = tuple(find_iupac)
    called_snp = iupac_dict.get(find_iupac)
    if not called_snp:
        # Combination of nucleotides not in dictionary- call uncertain- set to N
        called_snp = 'N'
    return called_snp


def multiple_snp_calls(ref_snp, called_snps):
    snps = [ref_snp]
    [snps.append(str(s)) for s in called_snps]
    # Handle cases where
    snps.sort()
    snps = tuple(snps)
    # Add appropriate IUPAC code
    find_iupac = (snps)
    called_snps = iupac_dict.get(find_iupac)
    if not called_snps:
        # Combination of nucleotides not in dictionary- call uncertain- set to N
        called_snps = 'N'
    return called_snps


def load_vcf(vcf_file_path):
    try:
        read_vcf = vcf.Reader(open(vcf_file_path, 'r'))
    except TypeError:
        raise FileNotFoundError("VCF file could not be loaded. Check path to VCF has been correctly entered on"
                                " the command line (flag -v) and that it is a VCF file.")
    # Update vcf metadata to add processing
    read_vcf.metadata.update({"vcf_edit": ['VCF edited to apply IUPAC uncertainty codes to ALT snp calls.']})
    return read_vcf


def parse_vcf(vcf, AF_data_location):
    vcf_record = []
    # Update ALT base for SNPs
    # Do not process indel calls
    for record in vcf:
        # Don't process indels
        if record.POS == 28988: # TEMP
            break # TEMP
        if record.is_indel:
            continue
        # Multisample vcfs not supported- TODO if needed
        if len(record.samples) != 1:
            raise NotImplementedError(f"Multisample VCFs not supported. Only one sample per record supported.")
        print(record.samples[0].data[AF_data_location])
        # Process snps
        if len(record.ALT) == 1:
            record.ALT = single_snp_call(record.REF, record.ALT)
        elif len(record.ALT) > 1:
            record.ALT = multiple_snp_calls(record.REF, record.ALT)
        else:
            raise Exception(f"No alt call at position {record.POS}. Cannot parse VCF") # TODO check this is actually doing anything
        vcf_record.append(record)
    return vcf_record



def main():
    __version__ = '0.0.1'
    __updated__ = '23/03/2020'
    args = get_args()
    vcf_data = load_vcf(args.path_to_vcf)
    # from bcftools the location of the allele depth (per allele)- TODO link to option for variant caller
    index_of_AD = 2
    updated_vcf_data = parse_vcf(vcf_data, index_of_AD)
    write_vcf(updated_vcf_data, vcf_data)


if __name__ == '__main__':
    main()
