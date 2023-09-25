"""
relation.py     Yota Kobayashi

A module for pitch-class set relations.

This module includes functions that take two pcsets to operate on or compute
a relation measure between them. For the functions of unary set operations,
use the operation module.

For the functions defined here, the parameters s1 and s2 are iterable of
type, that is, set, tuple, list, or Pcset object--an instance of Pcset class
from the module pcset--but not str, because each pc is represented by an int.

List of functions:

BINARY OPERATION
    union(s1, s2)
    difference(s1, s2)
    intersection(s1, s2)
    symmetricDifference(s1, s2)

TRANSFORMATION RELATION
    isTnEquivalent(s1, s2)
    isTnIEquivalent(s1, s2)
    pathSame(s1, s2)
    pathEmbed(s1, s2)
    pathCover(s1, s2)

INCLUSION RELATION
    isSubset(s1, s2)
    isSuperset(s1, s2)
    inclusion(s1, s2)

COMPLEMENT RELATION
    isComplement(s1, s2)
    complementation(s1, s2)

Z-RELATION
    isZRelated

SET-COMPLEX RELATION
    setComplexRelations(s1, s2)

SIMILARITY RELATION
    simRp(s1, s2)
    simIC(s1, s2)
"""

from itertools import combinations
from math import sqrt
from statistics import mean
from .pcset import Pcset
from .operation import interval, primeForm, subsets, complement, icv

__all__ = ['union',
           'difference',
           'intersection',
           'symmetricDifference',
           'isTnEquivalent',
           'isTnIEquivalent',
           'pathSame',
           'pathEmbed',
           'pathCover',
           'isSubset',
           'isSuperset',
           'inclusion',
           'isComplement',
           'complementation',
           'isZRelated',
           'setComplexRelations',
           'simRp',
           'simIC',
           'icvsim',
           'recrel']


# Binary operation functions ----------------------------------------------

def union(s1, s2):
    """
    Returns a new pcset containing pcs in either s1 or s2 with no duplicates.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a set resulted from the set union.
    """
    s1, s2 = set(s1), set(s2)
    return s1 | s2


def difference(s1, s2):
    """
    Returns a new pcset containing pcs in s1 that are not in s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a set resulted from the set difference.
    """
    s1, s2 = set(s1), set(s2)
    return s1 - s2


def intersection(s1, s2):
    """
    Returns a new pcset containing pcs in both s1 and s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a set resulted from the set intersection.
    """
    s1, s2 = set(s1), set(s2)
    return s1 & s2


def symmetricDifference(s1, s2):
    """
    Returns a new pcset containing pcs in either s1 or s2 but not both.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a set resulted from the set symmetrical difference.
    """
    s1, s2 = set(s1), set(s2)
    return s1 ^ s2


# Transformation relation functions ---------------------------------------

def isTnEquivalent(s1, s2):
    """
    A function to check whether the input two sets are transpositionally
    equivalent.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: an int for the transposition number with which s1 and s2
        are Tn equivalent, None if they are not Tn equivalent.
    """
    if len(s1) != len(s2):
        return None
    s1, s2 = Pcset(s1).normalForm(), Pcset(s2).normalForm()
    # When sets are transpositionally equivalent, they hold the same AIS.
    if interval(s1) == interval(s2):
        return (s2[0] - s1[0]) % 12
    else:
        return None


def isTnIEquivalent(s1, s2):
    """
    A function to check whether the input two pcsets are inversionally
    equivalent.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: an int for the index number with which s1 and s2 are TnI
        equivalent, None if they are not TnI equivalent.
    """
    if len(s1) != len(s2):
        return None
    s1, s2 = Pcset(s1).normalForm(), Pcset(s2).normalForm()
    ais1, ais2 = interval(s1), interval(s1)  # Create AIS
    # Sets are inversionally equivalent, if their AISs are mutually
    # retrogradable.
    if ais1 == list(reversed(ais2)):
        return (s2[-1] + s1[0]) % 12
    else:
        # Exceptions are some inversionally symmetrical set, in such case
        #   either set must be rotated to detect mutually retrogradable AIS.
        rots = [_rotate(s2, j) for j in range(len(s2))]
        for rot in rots:
            ais2 = interval(rot)
            if ais1 == list(reversed(ais2)):
                return (s1[0] + rot[-1]) % 12
        return None


def pathSame(s1, s2):
    """
    A function to find the operational paths to transform s1 into s2,
    that is, to make both sets comprise the same pc elements. If s1 can be
    transformed into s2 through Tn or TnI, this function will find the
    possible values for n.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a dict with two nested lists, {"Tn": [], "TnI": []}:
        the list Tn comprises possible values for n where Tn(s1) == s2
        the list TnI comprises possible values for n where TnI(s1) == s2.
    """
    return Pcset(s1).pathSame(s2)


def pathEmbed(s1, s2):
    """
    A function to find the operational paths to make s1 a literal subset of s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a dict with two nested lists, {"Tn": [], "TnI": []}:
        the list Tn comprises possible values for n where Tn(s1) < s2
        the list TnI comprises possible values for n where TnI(s1) < s2.
    """
    return Pcset(s1).pathEmbed(s2)


def pathCover(s1, s2):
    """
    A function to find the operational paths to make s1 a literal superset
    of the input set.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: a dict with two nested lists, {"Tn": [], "TnI": []}:
        the list Tn comprises possible values for n where Tn(s1) > s2
        the list TnI comprises possible values for n where TnI(s1) > s2.
    """
    return Pcset(s1).pathCover(s2)


# Inclusion relation functions --------------------------------------------

def isSubset(s1, s2):
    """
    A function to check the subset inclusion relation between s1 and s2,
    that is, whether s1 is a subset of s2 literally or abstractly.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: 0 (not a subset), 1 (literal subset), 2 (abstract subset).
    """
    if len(s1) >= len(s2):  # Pretest
        return 0
    s1, s2 = set(s1), set(s2)
    if s1 < s2:  # Literal subset
        return 1
    elif set(primeForm(s1)) in [set(primeForm(s)) for s in subsets(s2, len(s1))]:
        return 2
    else:
        return 0


def isSuperset(s1, s2):
    """
    A function to check the superset inclusion relation between s1 and s2,
    that is, whether s1 is a superset of s2 literally or abstractly.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: 0 (not a superset), 1 (literal superset), 2 (abstract superset).
    """
    if len(s1) <= len(s2):  # Pretest
        return 0
    s1, s2 = set(s1), set(s2)
    # s1 is a superset of s2, if s2 is a subset of s1.
    return isSubset(s2, s1)


def inclusion(s1, s2):
    """
    A function to find all the possible sets included in s1 that are the
    members of the set class of s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs representing the target set class.
    :return: a list of tuples, each tuple contains two sets (target, diff):
        target: a literal subset of the current set of the target set class
        diff: a set of difference pcs--pcs in the current set but not in target
    """
    return Pcset(s1).inclusion(s2)


# Complement relation functions -------------------------------------------

def isComplement(s1, s2):
    """
    A function to check the complement relation between s1 and s2, that is,
    whether s1 is a complement of s2 literally or abstractly.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: 0 (not a complement), 1 (literal complement), 2 (abstract complement)
    """
    if len(s1) + len(s2) != 12:  # Pretest
        return 0
    # Comparison is made in set data type
    s1 = set(s1)
    s2Comp = complement(s2)
    if s1 == s2Comp:  # Literal complement
        return 1
    elif primeForm(s1) == primeForm(s2Comp):  # Abstract complement
        return 2
    else:
        return 0


def complementation(s1, s2):
    """
    A function to find all the possible sets to complement s1 to form the
    members of the set class of s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs representing the target set class.
    :return: a list of tuples, each tuple contains two sets (target, comp):
        target: a complemented set that is a union of the current set and
            the set of complementing pcs s2
        comp: a set of complementing pcs for the current set to form
            the target set class
    """
    return Pcset(s1).complementation(s2)


# Z-relation function -----------------------------------------------------

def isZRelated(s1, s2):
    """
    A function to check the Z-relation between s1 and s2, that is, whether
    the input s1 and s2 are Z-related having the same ICV.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: True if Z-related, False otherwise.
    """
    return set(icv(s1)) == set(icv(s2))


# Set-complex relation ----------------------------------------------------

def setComplexRelations(s1, s2):
    """
    A function to check the set-complex relations between s1 and s2,
    that is, whether s1 is in a set-complex relation K or Kh about s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: 1 (set-complex relation K), 2 (set-complex relation Kh),
        0 (if none of set complex relations hold).
    """
    # Preliminary conditions of s1 and s2 are that their cardinalities
    #   are between 3 and 9 and not the same with each other.
    card1, card2 = len(s1), len(s2)
    if not(3 <= card1 <= 9) or not(3 <= card2 <= 9) or (card1 == card2):
        return 0
    s2Comp = complement(s2)
    bool1 = isSubset(s1, s2) or isSuperset(s1, s2)
    bool2 = isSubset(s1, s2Comp) or isSuperset(s1, s2Comp)
    if bool1 and bool2:
        return 2
    elif bool1 or bool2:
        return 1
    else:
        return 0


# Similarity relation -----------------------------------------------------

def simRp(s1, s2):
    """
    A function to check Forte's pitch-class similarity relation Rp between
    s1 and s2.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: True if Rp holds between s1 and s2, False otherwise.
    """
    if len(s1) != len(s2):  # s1 and s2 must have the same cardinality.
        return False
    s1 = set([tuple(primeForm(s)) for s in combinations(s1, len(s1) - 1)])
    s2 = set([tuple(primeForm(s)) for s in combinations(s2, len(s2) - 1)])
    return len(s1 & s2) != 0


def simIC(s1, s2):
    """
    A function to check Forte's interval-class relation between s1 and s2,
    that is, whether s1 and s2 are in a relation R1, R2, or R0.

        R1 (maximum similarity): same digits in ICV, 4 digits correspond.
        R2 (maximum similarity): same digits in 4 positions of ICV.
        R0 (minimum similarity): no ICV entries correspond.

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :return: 1 (R1), 2 (R2), 0 (R0), -1 (if none of R1, R2, or R0 holds).
    """
    if len(s1) != len(s2):  # s1 and s2 must have the same cardinality.
        return -1
    v1, v2 = icv(s1), icv(s2)
    nonEqual = []
    for i in range(6):
        if v1[i] != v2[i]:
            nonEqual.append((v1[i], v2[i]))
    if len(nonEqual) == 6:
        return 0
    elif len(nonEqual) == 2:
        if nonEqual[0][1] == nonEqual[1][0]:
            return 1
        else:
            return 2
    else:
        return -1


def icvsim(s1, s2, raw=False):
    """
    Return a set-class similarity value measured by Isaacson's ICVSIM
    (Isaacson, 1990).

    .. note::
        The original ICVSIM outputs a similarity value in the range [0.00,
        3.58] where a value of 0.00 indicates maximum similarity with
        respect to interval-class content. When raw=False, this method
        returns the raw similarity value; when raw=True, the similarity
        value is "flipped" and normalized so that it returns a float [0.00,
        1.00] where 1.00 indicates the maximum similarity.

    >>> s = {0,1,2}  # 3-1
    >>> icvsim(s, {0,1,2,3,6}, raw=True))    # ICVSIM(3-1, 5-4)
    0.3726779962499649
    >>> icvsim(s, {0,1,2,5,6,7}, raw=True))  # ICVSIM(3-1, 6-Z6)
    1.0
    >>> icvsim(s, {0,1,4,6,9}, raw=True))    # ICVSIM(3-1, 5-32)
    1.343709624716425
    >>> icvsim(s, {0,3,6,9}, raw=True))      # ICVSIM(3-1, 4-28)
    1.9790570145063195
    >>> icvsim({0,2,4,6,8,10}, {0,1,3,4,6,7,9,10}))  # ICVSIM(6-35, 8-28)
    0.0
    >>> icvsim({0,1,3,4,8}, {0,3,4,5,8}))  # ICVSIM(5-Z17, 5-Z37)
    1.0

    :param s1: an iterable with pcs.
    :param s2: an iterable with pcs.
    :param raw: If True, the returned ICVSIM value is in the range [0.0,
        3.578...], where 0.0 indicates the maximum similarity; if
        False (default), the value is flipped and normalized so
        that it is in the range [0.0, 1.0] where 1.0 indicates the maximum
        similarity.
    :return: ICVSIM value.
    :rtype: float
    """
    icv1, icv2 = icv(s1), icv(s2)
    idv = [icv2[i] - icv1[i] for i in range(6)]  # Interval-difference Vector
    idvMean = mean(idv)
    sim = sqrt(sum([(idv[i] - idvMean) ** 2 for i in range(6)]) / 6)  # ICVSIM

    if raw:
        return sim
    else:
        # Maximum dissimilarity returned from ICVSIM(SC6-35, SC8-28)
        maxDissim = 3.5784850922639815
        return (maxDissim - sim) / maxDissim


def recrel(s1, s2):
    """

    :param s1:
    :param s2:
    :return:
    """
    #TODO
    pass


# Private functions -------------------------------------------------------

def _rotate(lst, n):
    """A helper function to rotate a list by n items"""
    return lst[n:] + lst[:n]
