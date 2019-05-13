#!/usr/bin/env python3

# Small script to test how Aho-corasick implementation works with multiple identical sequences

import ahocorasick as aho_corasick  # Fast and memory efficient library for exact or approximate multi-pattern string search

sequence = "atatatat"
fragments = ["tata", "tata"]  # Two identical fragments

A = aho_corasick.Automaton()
for idx, key in enumerate(fragments):
    A.add_word(key, (idx, key))
A.make_automaton()

for end_index, (insert_order, fragment) in A.iter(sequence):
    position = end_index - len(fragment) + 1
    print("Found fragment #{} - {}, position {}".format(insert_order, fragment, position))
    
