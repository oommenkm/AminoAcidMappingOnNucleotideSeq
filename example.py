#!/usr/bin/env python

import AminoAcidMapper

amino_acid_mapper_obj = AminoAcidMapper.AminoAcidMapper()
mapping_details = amino_acid_mapper_obj.map("ATGCATGCATGCCGTAGCTAGCTAGCTAGTCGATCGTAGTCGATCGT", "MHACRS", 6)
for mapping_detail in mapping_details:
    for my_key in mapping_detail.keys():
        print my_key, mapping_detail[my_key]
    print "\n"