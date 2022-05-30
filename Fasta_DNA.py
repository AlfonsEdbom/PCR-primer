from os import path

from Trie import Trie


class Fasta_DNA:
    """"
    Reads Fasta file and stores the forward and reverse strands of the genome sequence
    """

    def __init__(self, fasta_file: str, length) -> None:

        self.length = length
        self.forward_strand = self.read_fasta(fasta_file)
        self.reverse_strand = self.calculate_reverse_DNA(self.forward_strand)



    def read_fasta(self, file_path: str) -> str:
        """
        Reads a Fasta file located inside the folder sequences
        Returns a string with forward strand
        """
        with open(file_path, 'r') as f:  # Read the file
            content = f.readlines()

        # Create string of DNA
        self.forward_strand = ""
        for i in range(1, len(content)):
            self.forward_strand += content[i].strip()

        return self.forward_strand[:self.length]

    def calculate_reverse_DNA(self, DNA_string: str) -> str:
        """
        Calculates the complimentary strand of a DNA sequence
        Returns the string in the 5'-3' direction
        """
        # Loop over DNA string from behind to create complimentary strand
        reversed_DNA = []
        for i in reversed(DNA_string):  # Match with correct opposing base
            if i == "A":
                reversed_DNA.append("T")
            elif i == "T":
                reversed_DNA.append("A")
            elif i == "G":
                reversed_DNA.append("C")
            elif i == "C":
                reversed_DNA.append("G")

        self.reverse_strand = "".join(reversed_DNA)  # Convert to string

        return self.reverse_strand

    def get_forward_strand(self):
        """Returns the forward strand of the genome in the 5'-3' direction"""
        return self.forward_strand

    def get_reverse_strand(self):
        """Returns the reverse strand of the genome in the 5'-3' direction"""
        return self.reverse_strand

    def build_primer_Trie(self, primer_length: int):
        """
        Split the genome into all possible primers of a specific length
        Returns a trie containing all primers found at least once in the genome
        """
        t = Trie()
        # Creates lists containing all primers found in forward and reverse strand
        forward_list = [self.forward_strand[i:i + primer_length]
                        for i in range(0, len(self.forward_strand) - primer_length + 1)]
        reverse_list = [self.reverse_strand[i:i + primer_length]
                        for i in range(0, len(self.reverse_strand) - primer_length + 1)]

        # Loops over both lists and add to trie
        for forward_primer, reverse_primer in zip(forward_list, reverse_list):
            t.insert(forward_primer)
            t.insert(reverse_primer)

        return t
