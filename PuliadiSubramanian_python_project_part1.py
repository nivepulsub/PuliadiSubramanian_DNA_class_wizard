#!/usr/bin/env python
# coding: utf-8

# Make a class called seq


class seq:

    # Call instance attributes
    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print(self.name)
        print(self.organism)
        print(self.type)
        print(self.sequence)

    # define the length function
    def length(self):
        print(len(self.sequence))

    # define the fasta_out function
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# Write the new protein child class here
class protein(seq):

    # Call instance attributes
    def __init__(self, name, organism, sequence, type, size):
        self.size = size
        super().__init__(name, organism, sequence, type)

    # define the fasta_out function
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_"
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()

    # define the mol_weight function
    def mol_weight(self):
        aa_mol_weights = {
            "A": 89.09,
            "C": 121.15,
            "D": 133.1,
            "E": 147.13,
            "F": 165.19,
            "G": 75.07,
            "H": 155.16,
            "I": 131.17,
            "K": 146.19,
            "L": 131.17,
            "M": 149.21,
            "N": 132.12,
            "P": 115.13,
            "Q": 146.15,
            "R": 174.2,
            "S": 105.09,
            "T": 119.12,
            "V": 117.15,
            "W": 204.23,
            "X": 0,
            "Y": 181.19,
        }
        total_mol_weight = 0
        for aa in self.sequence:
            total_mol_weight += aa_mol_weights[aa]
        print("The total molecular weight of the protein sequence:", total_mol_weight)


# Write the new nucleotide class here
class nucleotide(seq):

    # call the instance attributes
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    # define gc_content function
    def gc_content(self):
        count = 0
        for i in self.sequence:
            if i == "G" or i == "C":
                count = count + 1
        GC_content_percentage = (count / len(self.sequence)) * 100
        print("The GC_content_percentage is : ", str(GC_content_percentage))


# Write the DNA class here
class DNA(nucleotide):

    # call the instance attributes
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    # define transcribe function
    def transcribe(self):
        transcribed_sequence = self.sequence.replace("T", "U")
        print("The transcribed sequence is:", transcribed_sequence)

    # define print_set_of_3_nucleotides function to show codons in the reading frame
    def print_set_of_3_nucleotides(self, description, frame):
        out_frame = ""
        for i in range(len(frame)):
            if i % 3 == 0:
                out_frame = out_frame + " "
            out_frame = out_frame + frame[i]
        print(description + out_frame.strip())

    # define six_frames function
    def six_frames(self):
        forward_first_reading_frame = self.sequence[0:]
        forward_second_reading_frame = self.sequence[1:]
        forward_third_reading_frame = self.sequence[2:]
        self.print_set_of_3_nucleotides(
            "The first forward strand reading frame is:", forward_first_reading_frame
        )
        self.print_set_of_3_nucleotides(
            "The second forward strand reading frame is:", forward_second_reading_frame
        )
        self.print_set_of_3_nucleotides(
            "The third forward strand reading frame is:", forward_third_reading_frame
        )

        reverse_complement_seq = self.reverse_complement()
        print("The reverse complement of the sequence:", reverse_complement_seq)
        reverse_first_reading_frame = reverse_complement_seq[0:]
        reverse_second_reading_frame = reverse_complement_seq[1:]
        reverse_third_reading_frame = reverse_complement_seq[2:]
        self.print_set_of_3_nucleotides(
            "The first reverse strand reading frame is:", reverse_first_reading_frame
        )
        self.print_set_of_3_nucleotides(
            "The second reverse strand reading frame is:", reverse_second_reading_frame
        )
        self.print_set_of_3_nucleotides(
            "The third reverse strand reading frame is:", reverse_third_reading_frame
        )

    # define reverse_complement function
    def reverse_complement(self):
        complementary_base_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        reverse_complement_seq = ""
        for base in self.sequence:
            complementary_base = complementary_base_dict[base]
            # print("The complementary base is:",complementary_base)
            reverse_complement_seq = complementary_base + reverse_complement_seq
        # print("The reverse complement of the sequence:",reverse_complement_seq)
        return reverse_complement_seq

    # define fasta_out function
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# Write the RNA class here
class RNA(nucleotide):

    # call the instance attributes
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    # define start function
    def start(self):
        start_codon_position = self.sequence.find("AUG")
        print("The start codon starts at index:", start_codon_position)

    # define translate function
    def translate(self):
        standard_code = {
            "UUU": "F",
            "UUC": "F",
            "UUA": "L",
            "UUG": "L",
            "UCU": "S",
            "UCC": "S",
            "UCA": "S",
            "UCG": "S",
            "UAU": "Y",
            "UAC": "Y",
            "UAA": "*",
            "UAG": "*",
            "UGA": "*",
            "UGU": "C",
            "UGC": "C",
            "UGG": "W",
            "CUU": "L",
            "CUC": "L",
            "CUA": "L",
            "CUG": "L",
            "CCU": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "CAU": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGU": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "AUU": "I",
            "AUC": "I",
            "AUA": "I",
            "AUG": "M",
            "ACU": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "AAU": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "AGU": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GUU": "V",
            "GUC": "V",
            "GUA": "V",
            "GUG": "V",
            "GCU": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "GAU": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "GGU": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
        }

        start_codon_position = self.sequence.find("AUG")
        print("The start codon starts at index:", start_codon_position)
        prot = ""
        for i in range(start_codon_position, len(self.sequence), 3):
            codon = self.sequence[i : i + 3]
            if len(codon) < 3:
                break
            try:
                aa = standard_code[codon]
            except:
                aa = "?"
            prot += aa
            if aa == "*":
                break
        print("The protein sequence is:", prot)

    # define fasta_out function
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# test out methods for DNA
uidA_DNA = DNA(
    name="uidA_DNA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="Bacteria",
    type="DNA",
)
uidA_DNA.transcribe()
uidA_DNA.six_frames()
uidA_DNA.reverse_complement()
uidA_DNA.fasta_out()


# test out methods for RNA
uidA_RNA = RNA(
    name="uidA_RNA",
    sequence="CGCAUGUUACGUCCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
    organism="Bacteria",
    type="RNA",
)
uidA_RNA.fasta_out()
uidA_RNA.translate()


# test out methods for protein
uidA_protein = protein(
    name="uidA_protein",
    sequence="MLRPVETPTREIKK",
    organism="Bacteria",
    type="protein",
    size="50",
)
uidA_protein.fasta_out()
uidA_protein.mol_weight()
