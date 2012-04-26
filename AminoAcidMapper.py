"""

"""

import sys
#import OMORFFinder
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#from Bio import Translate

class AminoAcidMapper():
    
    def __init__ (self):
        pass
    
    def __find_orfs_with_trans(self, seq, trans_table, min_protein_length):
        answer = []
        seq_len = len(seq)
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                seq_obj = Seq(str(nuc[frame:]))
                prot_obj = seq_obj.translate(trans_table)
                prot_seq = str(prot_obj)
                trans_len = len(prot_seq)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = prot_seq.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end - aa_start >= min_protein_length:
                        if strand == 1:
                            start = frame + aa_start * 3
                            end = min(seq_len, frame + aa_end * 3 + 3)
                        else:
                            start = seq_len - frame - aa_end * 3 - 3
                            end = seq_len - frame - aa_start * 3
                        answer.append((start, end, strand, prot_seq[aa_start:aa_end]))
                    aa_start = aa_end + 1
        answer.sort()
        return answer

    def map(self, dna = "ATGCATGCATGCCGTAGCTAGCTAGCTAGTCGATCGTAGTCGATCGT", prt_seq = "MHACRS", nucleotide_num = 6):
        """
        This will return a list of dictionaries contains frame_number, codon_num, nucleotide_position and amino_acid
        """
        dna = re.sub("\W", "", dna)
        
        my_seq = Seq(dna, IUPAC.unambiguous_dna)
        
        my_protein = my_seq.translate()
        
        table = 11
        min_pro_len = 2
        
        n = nucleotide_num
        
        cur_fragment = my_seq
        orf_list = self.__find_orfs_with_trans(cur_fragment, table, min_pro_len)
        
        aam_result_list = []
        
        for start, end, strand, pro in orf_list:
            cur_ind = pro.find(prt_seq)
            if cur_ind == 0:
                aam_result_dict = {}
                codon_num = (nucleotide_num - 1 - start)  / 3
                aam_result_dict["frame_number"] = strand
                aam_result_dict["codon_num"] = codon_num + 1
                aam_result_dict["nucleotide_position"] = codon_num * 3
                aam_result_dict["amino_acid"] = pro[codon_num]
                aam_result_list.append(aam_result_dict)
        return aam_result_list