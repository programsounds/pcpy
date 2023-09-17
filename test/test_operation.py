# operation_test.py
# Testing script for functions defined in the operation module

from pcsets.operation import *


lst = [3, 5, 1, 10, 23, 16, 8]
print("\nPitch intervals:")
print(pitchInterval(lst))  # Unordered (directed) interval
print(pitchInterval(lst, ordered=True))  # Ordered (absolute) interval
print("\nOrdered pitch-class intervals:")
print(interval(lst))
print("\nUnordered pitch-class intervals (ics):")
print(intervalClass(lst))

s = [4, 0, 9, 11, 8]
print("\nTransposition in pc space:")
print(transpose(s, 5), opT(s, 5))
print("\nInversion:")
print(invert(s))
print(opTnI(s, 5))
print(opIxy(s, 2, 3))

s = [0, 1, 2, 6, 8]
print("\nInterval-class vector:")
print(icv(s))

s = [0, 1, 4, 8]
print("\nIndex vector")
print(indexVector(s))

s = [4, 0, 9, 11, 8]
print("\nNormal form:")
print(normalForm(s))

s = [0, 4, 6, 11, 1, 7, 10]
print("\nPrime form:")
print(primeForm(s))

s = [8, 9, 2, 5, 0, 1]
print("\nTransformation levels:")
print(transformationLevels(s))

s = [9, 10, 0, 5, 3, 2]  # 6-Z25 (T8I, O2' D)
print("\nReferential collections:")
print(referentialCollections(s))
print("\nComplement (12T)")
print(complement(s))
print("\nModal complements:")
print(modalComplements(s))
print("\nSubsets of cardinality n")
print(subsets(s, 5))

s = [10, 0, 1, 3, 4]  # 5-10 (T4I, O0)
print("\nTranspositional invariants:")
print(transpositionalInvariants(s, 9))
print("\nInversional invariants:")
print(inversionalInvariants(s, 1))

s = [3, 4, 9, 10]  # 4-9 (T5, O2)
print("\nTranspositional symmetry:")
print(transpositionalSymmetry(s))
print("\nInversional symmetry:")
print(inversionalSymmetry(s))
