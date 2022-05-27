import time
from datetime import timedelta

import json
import logging

from Fasta_DNA import Fasta_DNA
from Candidate_Primers import Candidate_Primers
from Primer_Pairs import Primer_Pairs



def get_config(config_file: str) -> dict:
    with open(config_file, "r") as c:
        config = json.load(c)

    return config

def get_logger(log_file: str):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger

def main():
    start_time = time.monotonic()  # Start timer for program

    # Load config settings from config.json and initiate logger
    config = get_config("config.json")
    logger = get_logger("log.log")
    # Set settings from config file to variables
    primer_length = config["settings"]["length"]
    GC_min = config["settings"]["GC_min"]
    GC_max = config["settings"]["GC_max"]
    T_min = config["settings"]["T_min"]
    T_max = config["settings"]["T_max"]
    deltaT = config["settings"]["deltaT"]
    GC_window = config["settings"]["GC_window"]

    # DNA-object
    genome_DNA = Fasta_DNA(config["files"]["SARS"])


    # Build trie
    t = genome_DNA.build_primer_Trie(primer_length)

    logger.debug(f"The Trie contains {len(t.query(''))} primers!")
    logger.debug(f"Primer Trie built successfully: {timedelta(seconds=time.monotonic() - start_time)}")

    # Remove bad primers
    candidate_primers = Candidate_Primers(genome_DNA)  # initiate filter object

    candidate_primers.filter_GC_content(GC_window, GC_min, GC_max)  # remove windows with too high GC-content

    candidate_primers.apply_filters(primer_length, T_min, T_max)  # Apply the rest of the filteres
    candidate_primers.remove_non_unique(t)
    candidate_primers.remove_low_complexity_primers()

    candidate_primers.remove_similar(t, deltaT)  # Remove primers with too low deltaT

    # Create primer pairs
    PPs = Primer_Pairs(genome_DNA, candidate_primers.forward_primers, candidate_primers.reverse_primers)
    PPs.find_primer_pairs(300, 1500, is_circular=True)
    PPs.filter_primer_pairs(GC_min, GC_max)

    # Cut with restriction enzymes
    PPs.restriction_enzymes_cut()

    primer_pairs = PPs.get_primer_pairs(3)


    print(f"The program took {timedelta(seconds=time.monotonic() - start_time)} to execute)")


if __name__ == '__main__':
    main()