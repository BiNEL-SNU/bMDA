#!/usr/bin/env python
import re
import os
import sys

def main():
    """
        import vcf file of SV and print out support vectors for each samples
    """
    #1: vcf file name #2: output_file_name

    curdir=os.getcwd()
    vcf_file_name= sys.argv[1]
    output_file_name = sys.argv[2]

    of = open(output_file_name, 'w')
    select_cols_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    select_cols_info = ['SUPP', 'SUPP_VEC', 'SVLEN', 'SVTYPE' ,'CHR2', 'END']

    header_dict = dict()
    col_names = []
    sample_names, sample_idx_start = [], 9
    with open(vcf_file_name) as fp:
        for line in fp:
            l = line.strip()
            if l.startswith('##'):
                # h_name, h_val = l[2:].split("=", 1)
                # if h_val.startswith('<') and h_val.endswith('>'):
                #     h_value_dict = {x.split("=", 1)[0] : x.split("=", 1)[1] for x in h_val[1:-1].split(",")}
                #     header_dict[h_name] = h_value_dict
                # else:
                #     header_dict[h_name] = h_val
                continue
            elif l.startswith('#'):
                col_names = l[1:].split("\t")
                print(f"found header : ", col_names[:sample_idx_start])
                sample_names = col_names[sample_idx_start:]
                print(f"found {len(sample_names)} samples in vcf : ", sample_names)
                of.write(",".join(select_cols_header + select_cols_info + sample_names) + "\n")
            else:
                vals = l.split("\t")
                qual, filter = vals[col_names.index("QUAL")], vals[col_names.index("FILTER")]
                info, format = vals[col_names.index("INFO")], vals[col_names.index("FORMAT")]
                info_dict = {x.split("=", 1)[0] : x.split("=", 1)[1] for x in info.split(";") if "=" in x} #ignore some cases
                format_list = format.split(":")
                # print(info_dict)
                print_vals = [vals[col_names.index(x)] for x in select_cols_header] + [info_dict[x] for x in select_cols_info]
                for s_name in sample_names:
                    sample_vals = vals[col_names.index(s_name)].split(":")

                    # if sample_vals[format_list.index("GT")] != "./.":
                    if sample_vals[format_list.index("PSV")] != "NaN":
                        print_vals.append('1')
                    else:
                        print_vals.append('0')
                of.write(",".join(print_vals) + "\n")



    of.close()
    return



if __name__ == "__main__":
    main()
