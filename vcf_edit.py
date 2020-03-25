#!/usr/bin/env python

import vcf
import argparse
import textwrap
import os
from collections import namedtuple

# Sub-class private class so can be edited
from vcf.model import _Call, make_calldata_tuple


# Sub-class private class so can be edited
class Caller(_Call):
    def caller(self):
        self.Parent.model_Call()


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
    try:
        for base_depth in depths_per_base:
            proportions.append(base_depth/total_depth)
    except TypeError:
        proportions.append(depths_per_base/total_depth)
    return proportions


def parse_out_sample_format(format_field, call_data):
    sample_format_data_dict = {}
    format_field_annotations = format_field.split(':')
    for i, annotation in enumerate(format_field_annotations):
        sample_format_data_dict[annotation] = call_data[i]
    return sample_format_data_dict


def get_total_bases_at_site(info_field):
    try:
        high_quality_bases = sum(info_field['DP4'])
    except TypeError:
        try:
            high_quality_bases = int(info_field['DP4'])
        except TypeError or ValueError:
            raise ValueError(f"Total number of high quality bases could not be determined from DP4 value. Proportion"
                             f" of bases supporting each read cannot be determined.")
    return high_quality_bases


def get_total_supporting_reads(reads):
    try:
        high_quality_bases = sum(reads)
    except TypeError:
        try:
            high_quality_bases = int(reads)
        except TypeError or ValueError:
            raise ValueError(f"Total number of high quality bases could not be determined from DP4 value. Proportion"
                             f" of bases supporting each read cannot be determined.")
    return high_quality_bases


def load_vcf(vcf_file_path):
    try:
        read_vcf = vcf.Reader(open(vcf_file_path, 'r'))
    except TypeError:
        raise FileNotFoundError("VCF file could not be loaded. Check path to VCF has been correctly entered on"
                                " the command line (flag -v) and that it is a VCF file.")
    # Update vcf metadata to add processing
    read_vcf.metadata.update({"vcf_edit": ['VCF edited to apply IUPAC uncertainty codes to ALT snp calls.']})
    return read_vcf


def parse_vcf(vcf):
    vcf_record = []
    # Update ALT base for SNPs
    # Do not process indel calls
    for record in vcf:
        # Don't process indels
        # Inbuild is-indel function finds reference calls (no alt) to be indels- bespoke identify indels
        is_indel = record.INFO.get('INDEL', False)
        # Do not process indels further
        if is_indel:
            # Do not propagate indels called through this workflow
            #vcf_record.append(record)
            continue

        # Process all sites that are single nucleotides
        # DP CANNOT BE USED if quality filtering is included in mpileup- calculate total depth from sun of DP4
        # Handle cases where there is only one DP4 value
        total_bases_at_site = get_total_bases_at_site(record.INFO)

        # Sites with minimum depth of coverage over all alleles at below threshold value are uncertain and set to N
        if total_bases_at_site < minimum_depth:
            record.ALT = iupac_dict.get('uncertainty')
            vcf_record.append(record)
            continue

        # Calculate high quality depth at site
        # Do not process reference where no alt was called sites further (note, these sites have no PL annotation)
        if record.is_indel:
            vcf_record.append(record)
            continue

        # Multisample vcfs not supported- TODO if needed
        if len(record.samples) != 1:
            raise NotImplementedError(f"Multisample VCFs not currently supported. "
                                      f" Only one sample per record supported.")

        # For variant sites passing minimum depth threshold, obtain proportion of reads supporting each base
        # Do not add an alt base if there is not one (and it is therefore None)
        all_bases = [str(base) for base in record.ALT if base is not None]
        all_bases.insert(0, record.REF)

        # Obtain correct FORMAT data for site
        format_field_dictionary = parse_out_sample_format(record.FORMAT, record.samples[0].data)
        try:
            supporting_reads = format_field_dictionary['AD']
        except KeyError as e:
            raise KeyError(f"Allele depth annotation not found for VCF. Unable to determine supporting reads per"
                           f" allele. Cannot filter VCFs without the AD annotation.")
        # Sanity Check
        total_supporting_reads = get_total_supporting_reads(supporting_reads)
        if total_bases_at_site != total_supporting_reads:
            raise ValueError(f"Total number of high quality bases could not be unambiguously determined. Proportion"
                             f" of bases supporting each read cannot be determined.")

        # Obtain proportion of reads supporting each base from numbers
        supporting_proportions = reads_to_proportions(total_bases_at_site, supporting_reads)

        # Filtering on proportion of supporting bases
        # TODO remove depth annotations for alleles that are no longer coming through- URGENT
        bases_to_process = []
        iupac_uncertainty = False
        for i, b in enumerate(all_bases):
            # Pass through only bases representing 0.75 or more of the total reads
            if supporting_proportions[i] > confident_allele_support_proportion:
                # Only one base can meet the above criteria
                print(record.POS)
                print(supporting_reads)
                print(supporting_reads[i])
                print(record.samples[0].data)
                print(record.samples)
                print(record.samples[0].site)
                caller_data = record.samples[0].data
                #caller_data = make_calldata_tuple(('GT', 'PL')) # requires input in the 'fields' format, which is the samp_fmt.split(':')
                print(caller_data)
                print(record.INFO)
                print(f"fields are {caller_data._fields}")
                caller_data = make_calldata_tuple(('GT', 'PL', 'AD'))
                ca = namedtuple('CallData', ['GT', 'PL', 'AD'])
                caller_data = ca(3, 1, 677)
                print(caller_data)
                record.samples[0].data = caller_data
                print(record.samples)
                c = Caller(record.samples[0].site, record.samples[0].sample, caller_data)
                print(c.data)
                print(c.data.GT)
                print(c._format_cache)
                c.data.GT = 3 #c.data.GT_replace(GT=3)
                print(c.data)
                #record.samples[0].data = {'GT':0, 'PL':3, 'AD':45}
                print(record.samples[0].data)
                # Edit INFO and FORMAT fields to remove reads where bases have been removed

                # If final call is the same as the reference base alt should be set to .
                if b == record.REF:
                    record.ALT = '.'
                else:
                    record.ALT = b
            # Sites with bases with support from >0.25 but <0.75 of the total reads are uncertain
            elif minimum_allele_support_proportion < supporting_proportions[i] < confident_allele_support_proportion:
                # Edit INFO and FORMAT fields to supply evidence for ambiguity base (summation of reads supporting)

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
    __updated__ = '24/03/2020'
    args = get_args()
    vcf_data = load_vcf(args.path_to_vcf)
    # Obtain default output filename
    outfile = f"{os.path.splitext(args.path_to_vcf)[0]}_edited.vcf"
    # Update/Filter vcf data where required
    updated_vcf_data = parse_vcf(vcf_data)
    # Write updated vcf file
    write_vcf(updated_vcf_data, vcf_data, outfile)


if __name__ == '__main__':
    main()
