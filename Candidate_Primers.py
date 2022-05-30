from Fasta_DNA import Fasta_DNA
from Trie import Trie
from Filters import Filters, calc_anneal_temp


class Primer:
    """
    Primer objects with information about each primer
    """

    def __init__(self, sequence: str, start: int, end: int) -> None:
        self.sequence = sequence
        self.start = start
        self.end = end
        self.anneal_temp = calc_anneal_temp(sequence)

    def __repr__(self):
        return f"{self.sequence}"

    def __str__(self):
        return f"{self.sequence}: ({self.start}-{self.end}): Annealing temperature: {self.anneal_temp}"


class Candidate_Primers:
    """
    Class with methods for removing primers that are determined to
    not be a candidate PCR-primer
    """

    def __init__(self, genome: Fasta_DNA):
        self.forward = genome.get_forward_strand()  # forward strand of DNA
        self.reverse = genome.get_reverse_strand()  # reverse strand of DNA
        self.filters = Filters()  # Filter object, with contains the logics for the filters
        self.forward_primers = []  # List of primers deemed to candidate forward PCR-primers
        self.reverse_primers = []  # List of primers deemed to be candidate reverse PCR-primers

    def filter_GC_content(self, window_size: int, GC_min: int, GC_max: int) -> None:
        """
        Uses a non-overlapping window to replace regions of the DNA with too high/low GC-content with '-'
        These regions will not be considered by the other functions
        """
        # Create non-overlapping windows for GC_content analysis
        forward_list = [self.forward[i:i + window_size]
                        for i in range(0, len(self.forward), window_size)]
        reverse_list = [self.reverse[i:i + window_size]
                        for i in range(0, len(self.reverse), window_size)]

        # Replaces bases with '-' if outside min and max
        for i, primer in enumerate(forward_list):
            primer = self.filters.remove_GCC(primer, GC_min, GC_max)
            forward_list[i] = primer

        for i, primer in enumerate(reverse_list):
            primer = self.filters.remove_GCC(primer, GC_min, GC_max)
            reverse_list[i] = primer

        # Rewrites the DNA-string with the GC-regions removed
        self.forward = "".join(forward_list)
        self.reverse = "".join(reverse_list)

    def apply_filters(self, primer_length: int, min_T: int, max_T: int, is_circular: bool = True):
        """
        Applies the filters found in Filters.py to all primers found in the genome
        Saves them in the list of candidate PCR-primers
        """

        for i in range(len(self.forward) - primer_length + 1):  # - primer_length + 1 to not get short primer in the end
            cur_primer = self.forward[i:i + primer_length]
            if ("-" in cur_primer):
                continue
            if not self.filters.GC_clamp(cur_primer):
                continue
            if not self.filters.annealing_temp(cur_primer, min_T, max_T):
                continue
            if not self.filters.GC_end(cur_primer):
                continue
            if not self.filters.self_dimerisation(cur_primer, 10):
                continue

            primer = Primer(cur_primer, i, i + primer_length)
            self.forward_primers.append(primer)


        for i in range(len(self.reverse) - primer_length + 1):
            cur_primer = self.reverse[i:i + primer_length]
            if ("-" in cur_primer):
                continue
            if not self.filters.GC_clamp(cur_primer):
                continue
            if not self.filters.annealing_temp(cur_primer, min_T, max_T):
                continue
            if not self.filters.GC_end(cur_primer):
                continue
            if not self.filters.self_dimerisation(cur_primer, 10):
                continue
            primer = Primer(cur_primer, -i, -(i + primer_length))
            self.reverse_primers.append(primer)

        if is_circular:  # If circular, add the end in front of the beginning
            # Start with forward strand
            for i in range(1, primer_length):
                end_part = self.forward[len(self.forward) - i:]  # Gets the i last bases in the sequence
                start_part = self.forward[:primer_length - i]  # Gets the rest of bases in primer from the beginning
                cur_primer = end_part + start_part  # combine end and start into string
                if ("-" in cur_primer):
                    continue
                if not self.filters.GC_clamp(cur_primer):
                    continue
                if not self.filters.annealing_temp(cur_primer, min_T, max_T):
                    continue
                if not self.filters.GC_end(cur_primer):
                    continue
                if not self.filters.self_dimerisation(cur_primer, 10):
                    continue

                primer = Primer(cur_primer, len(self.forward) - i, primer_length - i)  # indices ranging from the end over to beginning e.g. (999-18)
                self.forward_primers.append(primer)
            # Reverse strand
            for i in range(1, primer_length):
                end_part = self.reverse[len(self.reverse) - i:]
                start_part = self.reverse[:primer_length - i]
                cur_primer = end_part + start_part
                if ("-" in cur_primer):
                    continue
                if not self.filters.GC_clamp(cur_primer):
                    continue
                if not self.filters.annealing_temp(cur_primer, min_T, max_T):
                    continue
                if not self.filters.GC_end(cur_primer):
                    continue
                if not self.filters.self_dimerisation(cur_primer, 10):
                    continue
                primer = Primer(cur_primer, -len(self.reverse) - i,-primer_length - i)  # indicies ranging from the end over to the beginning, but negative e.g. ((-999)-(-18))
                self.reverse_primers.append(primer)

    def remove_non_unique(self, trie: Trie):
        all_primers = trie.query("")

        non_unique_primers = []
        for primer in all_primers:
            if primer[1] > 1:
                non_unique_primers.append(primer[0])

        tmp_flist = []
        for forward_primer in self.forward_primers:
            if forward_primer.sequence not in non_unique_primers:
                tmp_flist.append(forward_primer)

        tmp_revlist = []
        for reverse_primer in self.reverse_primers:
            if reverse_primer.sequence not in non_unique_primers:
                tmp_revlist.append(reverse_primer)

        self.forward_primers = tmp_flist
        self.reverse_primers = tmp_revlist

    def remove_low_complexity_primers(self):
        tmp_forward = []
        for primer in self.forward_primers:
            seq = primer.sequence
            if self.filters.low_complexity(seq):
                tmp_forward.append(primer)

        tmp_reverse = []
        for primer in self.reverse_primers:
            seq = primer.sequence
            if self.filters.low_complexity(seq):
                tmp_reverse.append(primer)

        self.forward_primers = tmp_forward
        self.reverse_primers = tmp_reverse

    def remove_similar(self, trie: Trie, max_mismatches: int):
        """
        Removes primers from self.candidate_primers that have too binding temp.
        """

        tmp_forward = []
        for primer in self.forward_primers:
            sequence = primer.sequence
            similar_list = trie.search_hamming_dist(sequence, max_mismatches)
            if not similar_list:
                tmp_forward.append(primer)

        tmp_reverse = []
        for primer in self.reverse_primers:
            sequence = primer.sequence
            similar_list = trie.search_hamming_dist(sequence, max_mismatches)
            if not similar_list:
                tmp_reverse.append(primer)

        self.forward_primers = tmp_forward
        self.reverse_primers = tmp_reverse


    def filter_primer_pairs(self, primer_pairs: list[tuple[Primer, Primer, int]], min_GC: int, max_GC: int):
        tmp_primers = []
        for primer_pair in primer_pairs:
            forward_primer = primer_pair[0]
            reverse_primer = primer_pair[1]
            distance = primer_pair[2]

            if self.filters.inter_dimerisation(forward_primer.sequence, reverse_primer.sequence):
                start_pos = forward_primer.start
                end_pos = forward_primer.start + distance

                contig = self.forward[start_pos:end_pos]
                num_G = contig.count("G")
                num_C = contig.count("C")
                num_A = contig.count("A")
                num_T = contig.count("T")

                num_GC = num_G + num_C
                num_AT = num_A + num_T

                GCC = num_GC / (num_GC + num_AT)

                if (min_GC / 100) < GCC < (max_GC / 100):
                    tmp_primers.append(primer_pair)

        primer_pairs = tmp_primers

        return sorted(primer_pairs, key=lambda x: x[2], reverse=True)
