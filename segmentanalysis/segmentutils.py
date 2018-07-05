# This module contains set of small service functions


def complement(sequence):
    """
    Make compliment of input sequence
    :param sequence: DNA sequence to make compliment
    :return:
    """
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([basecomplement[base] for base in sequence])


def revcomp(sequence):
    """
    Make reverse-compliment of input sequence
    :param sequence: DNA sequence to revcomp
    :return:
    """
    return complement(sequence)[::-1]
