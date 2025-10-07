#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 18:04:49 2025

@author: Nebeolisa Ernest Onwuneme
"""

import time
from Levenshtein import median, median_improve, distance as lev_distance

def read_sigma_from_file(I_10_20_0):
    
    """
    Reads n, m, and sigma (binary strings) from the given file.
    Expects format:
        n = length of the string
        m = number of string
        sigma = <binary rows>
    """
    with open(I_10_20_0, "r") as f:
        lines = f.read().strip().splitlines()

    # Extract n and m
    n = int(lines[0].split("=")[1].strip())
    m = int(lines[1].split("=")[1].strip())
    # Extract sigma rows (from line 3 onwards)
    sigma = []
    for line in lines[3:]:
        row = "".join(line.split())  # join digits into a single string
        sigma.append(row)

    return n, m, sigma 


 
#Function to find an approximate median string using Levenshtein.median. 
def find_median_string(str_list, weights=None, improve_steps=0):
    med = median(str_list, wlist=weights) if weights else median(str_list)
    for _ in range(improve_steps):
        med = median_improve(med, str_list, wlist=weights) if weights else median_improve(med, str_list)
    return med

#finding the total edit distance of median candidate and the strings
def total_distance(candidate, strings):
    """Sum of Levenshtein distances from candidate to all strings in S."""
    return sum(lev_distance(candidate, s) for s in strings)


def one_edit_neighbors(s, alphabet):
    """
    Generate all unique strings at Levenshtein distance 1 from s:
      - deletions
      - insertions (for every symbol in alphabet)
      - substitutions (for every symbol != s[i])
    """
    seen = set()
    # deletions
    for i in range(len(s)):
        t = s[:i] + s[i+1:] #take everything before position i and everything after position i, join them → the character at index i is removed.
        if t not in seen:
            seen.add(t)
            yield t

    # insertions
    for i in range(len(s) + 1):
        for c in alphabet:
            t = s[:i] + c + s[i:]
            if t not in seen:
                seen.add(t)
                yield t

    # substitutions
    for i in range(len(s)):
        orig = s[i]
        for c in alphabet:
            if c != orig:
                t = s[:i] + c + s[i+1:]
                if t not in seen:
                    seen.add(t)
                    yield t


def hybrid_greedy_iterative(m0, strings, alphabet, max_passes=1000):
    """
    Hybrid Greedy iterative local search:
    - Start from m0
    - Repeatedly scan the 1-edit neighborhood
    - Move to the best strictly improving neighbor
    - Stop at local optimum
    Returns: (best_string, best_cost, passes)
    """
    best = m0 
    best_cost = total_distance(best, strings)
    passes = 0

    while passes < max_passes:
        passes += 1
        improved = False
        candidate_best = best
        candidate_cost = best_cost

        for nbr in one_edit_neighbors(best, alphabet):
            cst = total_distance(nbr, strings)
            if cst < candidate_cost:
                candidate_best = nbr
                candidate_cost = cst
                improved = True

        if not improved:
            break

        best, best_cost = candidate_best, candidate_cost

    return best, best_cost, passes

"""
Calculates the Levenshtein edit distance between a median string
and each string in a list, and displays the results.
"""
def calculate_and_display_distances(median_string, string_list):
    print("\n Edit Distances for Hybrid Greedy Median String")
    total_distance = 0
    for original_string in string_list:
        distance = lev_distance(median_string, original_string)
        total_distance += distance
        print(f"Distance between '{median_string}' and '{original_string}': {distance}")
    print(f"Total distance for '{median_string}': {total_distance}") 

# ---------- Main flow ----------

if __name__ == "__main__":
    filename = "I_10_20_0.txt"
    n, m, sigma = read_sigma_from_file(filename)
    alphabet = sorted(set("".join(sigma)))  # e.g. ['0','1'] for binary

    print(f"n = {n}, m = {m}")
    print("Sigma =", sigma)
    print("Alphabet Σ =", alphabet)
    print()

    # Greedy median
    start = time.process_time()
    median_str = find_median_string(sigma)
    t_greedy = time.process_time() - start
    cost_greedy = total_distance(median_str, sigma)
    print(f"[Greedy] median = {median_str} | total cost = {cost_greedy} | CPU time = {t_greedy:.6f}s")

    # Refined median
    start = time.process_time()
    median_refined = find_median_string(sigma, improve_steps=2)
    t_refined = time.process_time() - start
    cost_refined = total_distance(median_refined, sigma)
    print(f"[Refined] median = {median_refined} | total cost = {cost_refined} | CPU time = {t_refined:.6f}s")

    # Hybrid greedy starting ONLY from refined median
    start = time.process_time()
    hybrid_med, hybrid_cost, passes = hybrid_greedy_iterative(median_refined, sigma, alphabet, max_passes=1000)
    t_hybrid = time.process_time() - start
    print(f"[Hybrid from refined] median = {hybrid_med} | total cost = {hybrid_cost} | passes = {passes} | CPU time = {t_hybrid:.6f}s")
    
    calculate_and_display_distances(hybrid_med, sigma)