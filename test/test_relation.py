# relation_test.py
# Testing script for functions defined in the relation module

from pcsets.relation import *

a = {0, 1, 3, 5, 6}
b = {1, 3, 4}
print("\nBinary operations:")
print(union(a, b))
print(difference(a, b))
print(intersection(a, b))
print(symmetricDifference(a, b))

a = {11, 8, 7, 5}
b = {4, 7, 1, 3}
print("\nTranspositional equivalence:")
print(isTnEquivalent(a, b))

a = {8, 0, 5, 4}
b = {0, 5, 1, 9}
print("\nInversional equivalence:")
print(isTnIEquivalent(a, b))

a = {1, 2, 5}
b = {7, 4, 3}
print("\nPath to the same set:")
print(pathSame(a, b))

a = {7, 4, 3}
b = {1, 2, 5, 9}
print("\nPath to the literal subset:")
print(pathEmbed(a, b))

a = {1, 2, 5, 11}
b = {5, 8, 9}
print("\nPath to the literal superset:")
print(pathCover(a, b))

a, b = {2, 3, 7}, {0, 1, 4, 5, 8}  # SC 3-4 and 5-21
print("\nSubset relation:")
print(isSubset(a, b))

a, b = {11, 0, 3, 4, 7}, {0, 3, 7}
print("\nSuperset relation:")
print(isSuperset(a, b))

a, b = {2, 3, 6, 7, 9, 10}, {3, 6, 7}  # 6-Z19 and 3-3
print("\nInclusion")
print(inclusion(a, b))

a, b = [0, 1, 2, 5, 6, 8, 9], [3, 4, 7, 10, 11]
c, d = [0, 1, 2, 5, 6, 8, 9], [0, 1, 4, 7, 8]
print("\nComplement relation:")
print(isComplement(a, b))
print(isComplement(c, d))

a, b = [0, 1, 4], [9, 0, 1, 3, 4]  # f = 5-16 (T4I, O0)
print("\nComplementation:")
print(complementation(a, b))

a, b = [1, 2, 5, 7], [3, 4, 6, 10]  # SC 4-Z15 and 4-Z29
print("\nZ-relation:")
print(isZRelated(a, b))

# (4-4, 7-6) = Kh
a = {0, 1, 2, 5}
b = {0, 1, 2, 3, 4, 6, 7}
print("\nSet-complex relation (K and Kh)")
print(setComplexRelations(a, b))
# (4-4, 7-5) = K
b = {0, 1, 2, 3, 5, 6, 7}
print(setComplexRelations(a, b))
# (4-4, 7-28) = none
b = {0, 1, 3, 5, 6, 7, 9}
print(setComplexRelations(a, b))

a = [10, 0, 1, 3, 4]
b = [9, 10, 0, 2, 3]
print("\nPitch-class relation (Rp)")
print(simRp(a, b))

a = [3, 4, 5, 7, 9]
b = [11, 1, 3, 5, 6]
print("\nInterval-class relations (R1, R2, R0)")
print(simIC(a, b))
a = [10, 0, 1, 3, 4]
b = [9, 10, 0, 2, 3]
print(simIC(a, b))
a = [3, 4, 5, 7]
b = [8, 11, 1, 2]
print(simIC(a, b))
