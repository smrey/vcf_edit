import vcf
import argparse
import textwrap

# Dictionary of IUPAC codes
iupac_dict = {('A', 'G'): 'R', ('C', 'T'): 'Y', ('C', 'G'): 'S', ('A', 'T'): 'W', ('G', 'T'): 'K', ('A', 'C'): 'M',
              ('C', 'G', 'T'): 'B', ('A', 'G', 'T'): 'D', ('A', 'C', 'T'): 'H', ('A', 'C', 'G'): 'V',
              ('A', 'C', 'G', 'T'): 'N', 'uncertainty': 'N'}

# This will use a lot of memory for large vcfs

minimum_depth = 20
minimum_allele_support_proportion = 0.25
confident_allele_support_proportion = 0.75

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
        '-d', '--min_depth', action='store_true', default=False,
        help=textwrap.dedent(
            '''
            Add user configurable minimum depth. Bases below this depth will be set to N Default value 20. OPTIONAL.
            '''))
    optional_options.add_argument(
        '-b', '--bcftools', action='store_true', default=False,
        help=textwrap.dedent(
            '''
           Specify that bcftools variant caller that has been used to generate the VCF. OPTIONAL.
            '''))

    return argument_parser.parse_args()


def write_vcf(records_to_write, read_vcf, output_vcf="output.vcf"):
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
        called_snp = iupac_dict.get('uncertainty')
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
        called_snps = iupac_dict.get('uncertainty')
    return called_snps


def apply_iupac(bases):
    # Order in alphabetical order for lookup
    bases.sort()
    bases = tuple(bases)
    find_iupac = bases
    called = iupac_dict.get(find_iupac)
    if not called:
        # Combination of nucleotides not in dictionary- call uncertain- set to N
        called = iupac_dict.get('uncertainty')
    return called



def reads_to_proportions(total_depth, depths_per_base):
    proportions = []
    for base_depth in depths_per_base:
        proportions.append(base_depth/total_depth)
    return proportions


def proportion_per_base(total_reads, reads_per_base):
    print(total_reads, reads_per_base)
    return None


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
        # Inbuild is-indel function finds reference calls (no alt) to be indels- bespoke identify indels
        is_indel = record.INFO.get('INDEL', False)
        # Do not process indels further
        if is_indel:
            vcf_record.append(record)
            continue

        # Process all sites that are single nucleotides
        # Sites with minimum depth of coverage over all alleles at below threshold value are uncertain and set to N
        # TODO Double check this threshold is per site not per base
        if record.INFO['DP'] < minimum_depth:
            record.ALT = iupac_dict.get('uncertainty')
            vcf_record.append(record)
            continue

        # Do not process reference where no alt was called sites further (note, these sites have no PL annotation)
        if record.is_indel:
            vcf_record.append(record)
            continue

        # Multisample vcfs not supported- TODO if needed
        if len(record.samples) != 1:
            raise NotImplementedError(f"Multisample VCFs not supported. Only one sample per record supported.")

        # For variant sites passing minimum depth threshold, obtain proportion of reads supporting each base
        all_bases = [str(base) for base in record.ALT]
        all_bases.insert(0, record.REF)
        #print(record.samples) TODO delete
        supporting_reads = record.samples[0].data[AF_data_location]
        supporting_proportions = reads_to_proportions(record.INFO['DP'], supporting_reads)

        # Filtering on proportion of supporting bases
        bases_to_process = []
        iupac_uncertainty = False
        for i, b in enumerate(all_bases):
            # Pass through only bases representing 0.75 or more of the total reads
            if supporting_proportions[i] > confident_allele_support_proportion:
                # Only one base can meet the above criteria
                if b == record.REF:
                    record.ALT = '.'
                else:
                    record.ALT = b
            # Sites with bases with support from >0.25 but <0.75 of the total reads are uncertain
            elif minimum_allele_support_proportion < supporting_proportions[i] < confident_allele_support_proportion:
                bases_to_process.append(b)
                iupac_uncertainty = True

        # Process uncertainty sites
        if iupac_uncertainty:
            record.ALT = apply_iupac(bases_to_process)

        # Update vcf records
        vcf_record.append(record)
    return vcf_record


def main():
    __version__ = '0.0.1'
    __updated__ = '23/03/2020'
    args = get_args()
    vcf_data = load_vcf(args.path_to_vcf)
    # from bcftools the location of the allele depth (per allele)- TODO link to option for variant caller
    index_of_AD = 2 # TODO unpack these values to a dict so don't need to hardcode this in and can look up AD
    updated_vcf_data = parse_vcf(vcf_data, index_of_AD)
    write_vcf(updated_vcf_data, vcf_data)


if __name__ == '__main__':
    main()
