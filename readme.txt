
Repository for the Verification of Integral Distinguishers:

This repository contains the resources used for verifying integral distinguishers
as described in the paper "Mind the Composition of Toffoli Gates:
Structural Algebraic Distinguishers of ARADI."

Contents:

- C++ Code for ARADI Encryption and Distinguisher Verification
  This includes the implementation for ARADI encryption and the
  verification of integral distinguishers for low dimension cubes.

  Command to compile and run:
  g++ -std=c++11 aradi_encryption.cpp

- Division Property Models for ARADI Block Cipher:
  1. Compute the algebraic degree of Boolean functions for output bits after r rounds.
  2. Compute the degrees of products of two Boolean functions for output bits after r rounds.
  3. Use the provided Makefile to run the above computations.
