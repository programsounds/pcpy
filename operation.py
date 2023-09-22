"""
operation.py    Yota Kobayashi

A module for pitch-class set operations.

This module includes functions that take a single pcset to operate on,
while the relation module takes two pcsets and implements relational
functions between them.

For the functions defined here, the parameter, pcs, is an iterable type:
set, tuple, list, or Pcset object--an instance of Pcset class from the
module pcset--but not str because each pc is represented by an int.

List of functions:

INTERVAL CALCULATION
    pitchInterval(pseg, ordered=)
    interval(pcseg)
    intervalClass(pcseg)

SET TRANSFORMATION
    transpose(pcs, n)
    invert(pcs)
    invertXY(pcs, x, y)
    opT(pcs, n)
    opI(pcs)
    opTnI(pcs, n)
    opIxy(pcs, x, y)

SET PROFILE
    icv(pcs)
    indexVector(pcs)
    normalForm(pcs)
    primeForm(pcs)
    transformationLevels(pcs)
    referentialCollections(pcs)

SET ANALYSIS
    complement(pcs)
    modalComplements(pcs)
    subsets(pcs, n)
    transpositionalInvariants(pcs, n)
    inversionalInvariants(pcs, n)
    transpositionalSymmetry(pcs)
    inversionalSymmetry(pcs)
"""

from pcpy.pcset import Pcset

__all__ = ['pitchInterval',
           'interval',
           'intervalClass',
           'transpose',
           'invert',
           'invertXY',
           'opT',
           'opI',
           'opTnI',
           'opIxy',
           'icv',
           'indexVector',
           'normalForm',
           'primeForm',
           'transformationLevels',
           'referentialCollections',
           'complement',
           'modalComplements',
           'subsets',
           'transpositionalInvariants',
           'inversionalInvariants',
           'transpositionalSymmetry',
           'inversionalSymmetry']


# Interval calculation functions ----------------------------------------------

def pitchInterval(pseg, ordered=False):
    """
    A function to compute pitch intervals--ordered or unordere distance
    between a series of pitches.

    :param pseg: a list of ordered pitches, representing a series of notes
        in the pitch space.
    :param ordered: True for ordered (i.e., directed) pitch intervals,
        and False (default) for unordered pitch intervals.
    :return: a list of pitch intervals.
    """
    intervals = []
    for i in range(len(pseg) - 1):
        intvl = pseg[i + 1] - pseg[i]
        intervals.append(intvl)
    if ordered:
        return intervals
    else:
        return [abs(j) for j in intervals]


def interval(pcseg):
    """
    A function to compute ordered pitch-class intervals.

    :param pcseg: a list of ordered pcs, representing a series of notes in
        the pitch-class space.
    :return: a list of ordered pitch-class intervals.
    """
    intervals = []
    for i in range(len(pcseg) - 1):
        intvl = (pcseg[i + 1] - pcseg[i]) % 12
        intervals.append(intvl)
    return intervals


def intervalClass(pcseg):
    """
    A function to compute unordered pitch-class intervals (i.e., interval
    classes).

    :param pcseg: a list of unordered pcs, representing a series of notes in
        the pitch-class space with interval classes.
    :return: a list of integers in the range from 0 to 6.
    """
    ics = []
    for i in range(len(pcseg) - 1):
        intvl1 = (pcseg[i + 1] - pcseg[i]) % 12
        intvl2 = (pcseg[i] - pcseg[i + 1]) % 12
        if intvl1 == intvl2:
            ics.append(6)
        else:
            ics.append(min(intvl1, intvl2))
    return ics


# Set transformation functions ------------------------------------------------

def transpose(pcs, n):
    """
    A function to transpose the input pcset by the interval n.

    :param pcs: an iterable with pcs.
    :param n: an int (mod 12) for the transposition number.
    :return: a set of transposed pcs.
    """
    return {(pc + n) % 12 for pc in pcs}


def invert(pcs):
    """
    A function to invert the input set around pc 0.

    :param pcs: an iterable with pcs.
    :return: a set of inverted pcs.
    """
    return {(12 - pc) % 12 for pc in pcs}


def invertXY(pcs, x, y):
    """
    A function to invert the input pcset around an axis of symmetry
    specified by pcs x and y. As a result, pcs x and y map onto each other.
    This is an operation equivalent to TnI where n = x + y. The advantage of
    using I(x,y) over TnI is that it allows for specifying the axis of
    symmetry.

    :param pcs: an iterable with pcs.
    :param x: an int that inverts x onto y.
    :param y: an int that inverts y onto x.
    :return: a set of inverted pcs.
    """
    return opTnI(pcs, (x + y) % 12)


def opT(pcs, n):
    """
    Shorthand for transpose(n)
    """
    return transpose(pcs, n)


def opI(pcs):
    """
    Shorthand for invert()
    """
    return invert(pcs)


def opTnI(pcs, n):
    """
    Shorthand for TnI operation--inversion followed by transposition at n.
    """
    return transpose(invert(pcs), n)


def opIxy(pcs, x, y):
    """
    Shorthand for invertXY(x, y)
    """
    return invertXY(pcs, x, y)


# Set profile functions -------------------------------------------------------

def icv(pcs):
    """
    Return the interval-class vector (ICV) of the input set.

    ICV is an enumeration of unordered pc intervals (i.e., ics).

    # FIXME: icv method in pcpy.Pcset will have an option to use this
        function (recursive calculation of icv) as the mode 2.
        Pcset.icv() has its implementation already which is mode 1.
        As the module import hierarchy is that operation.py imports
        pcset.py, move the icv() function body below to pcset.py.
        Then, the function can be called as Pcset(pcs).icv(mode=2).

    :param pcs: an iterable with pcs.
    :return: a list for the ICV.
    """
    pcs = sorted(set(pcs))  # Remove duplicates and sort the pcs
    if len(pcs) == 1:  # Base case
        return [0] * 6
    else:
        vec = [0] * 6
        for pc in pcs[1:]:
            x = pc - pcs[0]
            ic = min(x, 12 - x)  # Convert ordered pc interval to ic
            vec[ic - 1] += 1     # Increment ICV digit
        return list(map(lambda a, b: a + b, vec, icv(pcs[1:])))


def indexVector(pcs):
    """
    A function to compute the index vector of the input pcset.

    :param pcs: an iterable with pcs.
    :return: a list of the index vector.
    """
    return Pcset(pcs).indexVector()


def normalForm(pcs):
    """
    A function to compute the normal form of the input set.

    :param pcs: an iterable with pcs.
    :return: a list of the normal form.
    """
    return Pcset(pcs).normalForm()


def primeForm(pcs):
    """
    A function to compute the prime form of the input pcset.

    :param pcs: an iterable with pcs.
    :return: a list of the prime form.
    """
    return Pcset(pcs).primeForm()


def transformationLevels(pcs):
    """
    A function to compute the Tn/TnI transformation level of the input pcset.

    :return: a dict with two nested lists, {"Tn": [], "TnI": []},
        where the lists associated with the keys Tn and TnI represent
        the transformation levels of the current set.
        For transpositionally and inversionally symmetrical sets,
        there would be multiple entries in the lists.
    """
    return Pcset(pcs).transformationLevels()


def referentialCollections(pcs):
    """
    A function to check the reference status to the modal collections,
    OCT, WT, HEX, and DT.

    :return: a dict with abbreviated keys for OCT, WT, and HEX as c+n,
        where c is the collection's initial and n is the transposition
        level (e.g., O0, W1, H3, etc.).
        These keys take an int of 3 possible values:
            0       no reference
            1       all but one pc are literal subset of the Tn level
                        of the collection
            2       all the pcs are literal subset of the Tn level of
                        the collection
        An additional key is for DT as D which takes a boolean value.
            True    if the current set is an abstract subset of DT
            False   otherwise
    """
    return Pcset(pcs).referentialCollections()


# Set analysis function -------------------------------------------------------

def complement(pcs):
    """
    Returns a set of the complement of the current pcset.

    :param pcs: an iterable with pcs.
    :return: a set of the complement of the current pcset.
    """
    return Pcset(pcs).complement()


def modalComplements(pcs):
    """
    Returns the complements to each of the modal collections OCT, WT, and HEX.

    :param pcs: an iterable with pcs.
    :return: a dict with abbreviated keys for OCT, WT, and HEX as c+n,
        where c is the collection's initial and n is the transposition
        level (e.g., O0, W1, H3, etc.).
        Each key has a value for the complementing pcs to the
        modal collection.
    """
    return Pcset(pcs).modalComplements()


def subsets(pcs, n):
    """
    Returns all the subsets of cardinality n of the input pcset.

    :param pcs: an iterable with pcs.
    :param n: an int for the cardinality of the subsets.
    :return: a list containing the subsets of cardinality n.
    """
    return Pcset(pcs).subsets(n)


def transpositionalInvariants(pcs, n):
    """
    Returns the invariants held under transposition of a pcset by interval n.
    The number of invariants is equal to the number of times the interval n
    occurs in the pcset. The digits of ICV represent intervals n and 12-n.
    Double the count of ic6.

    :param pcs: an iterable with pcs.
    :param n: an int for the ordered pc interval of transposition.
    :return: a set of invariant pcs.
    """
    return Pcset(pcs).transpositionalInvariants(n)


def inversionalInvariants(pcs, n):
    """
    Returns the invariants held under inverting a pc set by interval n. The
    number of invariants can be checked in the index vector of the set.

    :param pcs: an iterable with pcs.
    :param n: an int for the index number.
    :return: a set of invariant pcs.
    """
    return Pcset(pcs).inversionalInvariants(n)


def transpositionalSymmetry(pcs):
    """
    Returns the degree of transpositional symmetry.

    :param pcs: an iterable with pcs.
    :return: a list of ordered pc intervals--the input pcset maps entirely
        onto itself when transposing at the transpositional level(s).
        The number is at least 1, because T0 is an identity operator.
    """
    return Pcset(pcs).transpositinalSymmetry()


def inversionalSymmetry(pcs):
    """
    Returns the degree of inversional symmetry.

    :param pcs: an iterable with pcs.
    :return: a list of ordered pc intervals--the input pcset maps entirely
        onto itself when inverting with the index number(s).
        None is output when the set is not inversionally symmetrical.
    """
    return Pcset(pcs).inversionalSymmetry()
