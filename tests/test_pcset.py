# pcset_test.py
# Testing script for the class defined in the pcset module.

from pcpy.pcset import *
from pcpy.operation import *

set1 = {4, 0, 1}
set2 = [0, 3, 4, 5, 7, 8]

print("\nMethod chaining:")
s1 = Pcset(set1)
s2 = Pcset(set2)
print(s1.transpose(4).union(s2.transpose(2)))

print(s1.getSet(), "test")

print("\nObject copying:")
s1 = Pcset(set1)
s2 = Pcset(set2)
s = s1.copy().transpose(5)
print(s1, s)  # s1 is intact
print(s1.union(s2))  # Argument can be a Pcset object

print("\nInclusion:")
s2 = Pcset(set2)
for target, diff in s2.inclusion(set1):
    print(target, diff)
    print(primeForm(target), primeForm(diff))

print("\nComplementation:")
s1 = Pcset(set1)
for target, comp in s1.complementation(set2):
    print(target, comp)
    print(primeForm(target), primeForm(comp))
