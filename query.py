"""
query.py    Yota Kobayashi

This module includes functions that collectively generate a JSON file
for queries related to pc sets and modal set complexes.

A fresh copy of the JSON file can be made through directly running
this module script--a new catalog.json will be created in the same
directory where this program is in.

This module decodes the JSON file and assign the data to the
variable, catalog. So the queries can be made through importing the
variable in other scripts.

The variable, catalog, is a dict with nested dict structures.
The nested dicts are associated with four top-level keys, PFToSN, SC,
inclusionTable, and MSC. Each key corresponds with the subject of query.
Note that it is always str data type for the keys of the catalog dict.

PFToSN (key=PF, value=SN)

    A mapping for prime form and set name.

    PF: a str for prime form (e.g., "02468T")
    SN: a str for set name (e.g., "6-35")

    The key PF is a string with pc characters. To facilitate the query
    with list form of PF, use the function of this module, toPFStr().

SC (key=SN, value={PF, ICV, Z-corr, symmetry, MA, MSC})

    A mapping for SCs and their profiles. The top-level value is
    a nested dict.

    PF: a str with pc characters for prime form (10 is substituted by T)
    SN: a str for set name
    ICV: a list of ints representing an ICV
    Z-corr: a str for Z-correspondent represented by a str of SN if any,
        None otherwise
    symmetry: a list of two ints for the degrees of transpositional
        and inversional symmetries
    MA: a list with string elements for modal attribute symbols
    MSC: a list of MSC nexus sets, represented by SNs, to which the SC
        is a member.

    The top-level query is made by SN. If the current data structure
    of a set is an iterable of prime form, use toPFStr() and PFToSN to
    obtain the SN of the set.

inclusionTable (key=SN, value={{n: [inclVec]}, ..., {n+i: [inclVec]})

    A mapping for inclusion tables of SCs. The top-level value is a
    nested dict that represents the inclusion table of the SC.
    Each inclusion table is also a dict that is an inclusion vector
    of a given cardinality.

    SN: a str for set name
    n, ..., n+i: ints for the cardinalities between 3 and 9
        except for that of the query target
    inclVec: a list with ints representing the inclusion vector

    For a set-class S of cardinality j and inclusion vector of
    cardinality k, if j < k, the vector gives the number of
    abstractly included set classes; if j > k, the vector gives the
    number of distinct members of set-classes that include S.

    Note that it is important to distinguish the covering (j > k) of S
    by a set-class Q from the embedding (j < k) of S in Q: the numbers
    are not always the same. In some cases, counts for covering is
    smaller than those of embedding. This is due to the elimination of
    Tn/TnI equivalent sets in the case of covering. Therefore, to find
    how many instances of S are embedded in Q, one looks at the entry
    for Q; if one wants to know the number of members of Q that cover S,
    one consults the entry for S.

MSC ({key=SN, value={card: [{member1}...{memberN}]})

    A mapping for the modal set complexes and their profiles. The
    nested dict represents MSC members grouped by their cardinalities.

    SN: a str for set name of the nexus set (i.e., "5-10", "5-16",
        "5-19", "5-21", "5-25", "5-28", "5-31", "5-32", "5-33")
    card: an int for cardinalities of the MSC members (i.e., 3, 4, 6)

    The value of the card is a list of the MSC members, and each
    member is represented by a dict {SN, MA, symmetry, inclusion, Z-corr}:

    SN: a str for a set name
    MA: a list with string elements for modal attribute symbols
    symmetry: a list of two ints for the degrees of transpositional
        and inversional symmetries.
    inclusion: an int for the number of inclusion count
    Z-corr: a str for Z-correspondent represented by a str of SN if any,
        None otherwise

List of functions:

toPFStr(s)
    Converts a prime form in any iterable form into a textual representation:
    e.g., [0, 1, 3, 4, 6, 8, 10] to "012568T".

fromPFStr(sn)
    Inverse of toPFStr() function, and converts the prime form from
    string format to a list format: e.g., "012568T" to [0, 1, 3, 4, 6, 8, 10].
"""

import os.path
import json
from itertools import combinations
from pcsets.operation import (complement, primeForm, icv, referentialCollections,
                              opT, opTnI, transpositionalSymmetry, inversionalSymmetry)
from pcsets.relation import setComplexRelations, isSubset
from pcsets import constants as c

__all__ = ["toPFStr", "fromPFStr", "makePFLists", "catalog"]

filename = os.path.join(os.path.dirname(__file__), "catalog.json")
with open(filename, "r") as f:
    catalog = json.load(f)


# Public functions ------------------------------------------------------------

def toPFStr(s):
    """
    Converts a prime form in any iterable form into the textual representation:
    e.g., [0, 1, 3, 4, 6, 8, 10] to "012568T".

    :param s: an iterable of a prime form.
    :return: a str of characters representing pcs is a prime form.
        int 10 and 11 are substituted by T and E in the PF string.
    """
    # Substitute str T and E for int 10 and 11
    pf = ["T" if pc == 10 else "E" if pc == 11 else pc for pc in s]
    return "".join(str(pc) for pc in pf)  # Convert to string


def fromPFStr(sn):
    """
    Inverse of toPFStr() function, and converts the prime form from the
    string format to the list format: e.g., "012568T" to [0, 1, 3, 4, 6, 8, 10].

    :param sn: a str representing a prime form.
    :return: a list representing the same prime form as the input.
    """
    # Substitute int 10 and 11 for str T and E
    return [10 if ch == "T" else 11 if ch == "E" else int(ch) for ch in sn]


# Functions for creating the JSON file ----------------------------------------

def mapPFToSN(pfs):
    """
    Returns a mapping for prime form and set name.

    :param pfs: a sorted list of sets which correspond with all the PFs of
        cardinality from 3 to 9.
    :return: a dict (key=PF, value=SN). PF and SN are both str.
    """
    lst = []  # list for tuples (key=SN, value=PF)
    for card in range(3, 10):
        ord_ = 1
        sets = pfs[card-3]
        for s in sets:
            pf = toPFStr(s)
            sn = str(card) + "-" + c.ORDINAL_NUMS[card][ord_]  # Make SN str
            lst.append((pf, sn))
            ord_ += 1
    return dict(lst)


def mapSC(pfs, dct):
    """
    Returns a mapping for SC and its profile.

    :param pfs: a sorted list of sets which correspond with all the PFs of
        cardinality from 3 to 9.
    :param dct: a dict for PF to SN conversion.
    :return: a dict (key=SN, value={PF, ICV, Z-corr, symmetry, MA, MSC}),
        where the keys in the nested dictionary are:

        PF: a str with pc characters for prime form (10 is substituted by T).
        SN: a str for set name.
        ICV: a tuple for ICV.
        Z-corr: a str for Z-correspondent represented by a SN if any,
            None otherwise.
        symmetry: a list of two ints for the degrees of transpositional
            and inversional symmetries.
        MA: a tuple with string elements for modal attributes symbols.
        MSC: a list of MSC nexus sets, represented by SNs, to which the SC
            is a member.
    """
    lst = []  # list for (key=SN, value=[PF, ICV, Z-corr, MA])
    for card in range(3, 10):
        ord_ = 1
        sets = pfs[card-3]
        for s in sets:
            # SN
            sn = str(card) + "-" + c.ORDINAL_NUMS[card][ord_]
            # PF
            pf = toPFStr(s)
            # ICV
            icvec = icv(s)
            # Z-correspondent
            zcorr = None
            others = [t for t in sets if t != s]
            for zcand in others:
                if icvec == icv(zcand):
                    zcorr = dct[toPFStr(zcand)]
            # Degrees of symmetry
            degrees = [len(transpositionalSymmetry(s)),
                       len(inversionalSymmetry(s))]
            # Modal attributes
            matts = []
            ref = referentialCollections(s)
            for col in c.REF_COLS:
                sym = col[0]  # Extract only the modal collection symbol
                if ref[col] == 2:
                    if sym not in matts:
                        matts.append(sym)
                # Check prime reference only for hexachords
                elif card == 6 and ref[col] == 1:
                    sym += "'"
                    if sym not in matts:
                        matts.append(sym)
            if ref["D"]:
                matts.append("D")
            # MSC membership
            mscs = []
            for nexusSN, nexusPF in c.NEXUS_SETS:
                rels = setComplexRelations(s, nexusPF)
                if (int(card) <= 5 and rels == 2) \
                        or (int(card) == 6 and rels >= 1 and isSubset(nexusPF, s)):
                    mscs.append(nexusSN)
            # Consolidate the dict entry for the current SC
            lst.append((sn, {"PF": pf, "ICV": tuple(icvec),
                             "Z-corr": zcorr, "symmetry": degrees,
                             "MA": tuple(matts), "MSC": mscs}))
            ord_ += 1
    return dict(lst)


def countInclusions(s1, s2):
    """
    A function to compute the number of instances s1 includes s2 and vice versa:
    given that SC1 and SC2 are the set-classes of s1 and s2 respectively,
    if s1 is larger than s2, the function returns the number of instances of
    each included SC2 in SC1, namely the number of abstract subsets of SC2 in SC1;
    whereas, if s1 is smaller than s2, the function returns the number of
    unique members of SC2 that include SC1--"unique" denotes elimination
    of duplicates from equivalent Tn and TnI members of SC2.

    :param s1: a list for a prime form.
    :param s2: a list for a prime form.
    :return: an int for the number of inclusions.
    """
    card1, card2 = len(s1), len(s2)
    if card1 == card2:
        return 0
    count = 0  # Inclusion count
    if card1 > card2:
        for s in combinations(s1, card2):
            if primeForm(s) == s2:
                count += 1
    else:
        sets = []
        for n in range(12):
            sets.append(tuple(sorted(opT(s2, n))))
            sets.append(tuple(sorted(opTnI(s2, n))))
        sets = set(sets)  # Eliminate Tn and TnI equivalent sets
        s1 = set(s1)
        for s in sets:
            if set(s) > s1:
                count += 1
    return count


def mapInclusionTable(dct):
    """
    Returns a mapping for inclusion tables of each SC.

    :param dct: a dict of SC profiles for SN to PFStr conversion.
    :return: a dict {key=SN, value={{n: [inclVec]}, ..., {j: [inclVec]}}},
        where the nested dict represents the inclusion table of the SC.
        Each inclusion table is also a dict with key=card and value=inclVec.
    """
    tables = {}  # dict for all the inclusion tables
    for i in range(3, 10):
        for sn1 in c.SN_VECS[i][1:]:
            table = {}  # Inclusion table for the current sn1
            s1 = fromPFStr(dct[sn1]["PF"])  # list form of PF
            for j in [k for k in range(3, 10) if k != i]:
                card = str(j)
                vec = []  # Inclusion vector for the current card
                for sn2 in c.SN_VECS[j][1:]:
                    s2 = fromPFStr(dct[sn2]["PF"])  # list form of PF
                    vec.append(countInclusions(s1, s2))
                table[card] = vec
            tables[sn1] = table
    return dict(tables)


def mapMSC(dctSC, dctIncl):
    """
    Returns a mapping for the modal set complexes.

    :param dctSC: a dict of SC profiles.
    :param dctIncl: a dict of inclusion tables.
    :return: a dict {key=SN, value={card: [{member1}...{memberN}]}},
        where the nested dict represents MSC members grouped by
        their cardinalities:

        card: cardinalities of the MSC members (i.e., 3, 4, 6)

        The value of the card is a list of the MSC members, and each
        member is represented by a dict
        {SN, MA, symmetry, inclusionCount, Z-corr}:

        SN: a string for set name
        MA: a list with strings for modal attributes
        symmetry: a list of two ints for the degrees of transpositional
            and inversional symmetries.
        inclusionCount: an int for the number of inclusion count
        Z-corr: a str for the name of Z-correspondent represented by
            a SN if any, None otherwise
    """
    dct = {}
    for nexusSN, nexusPF in c.NEXUS_SETS:
        profile = {}  # dict for the profile of the current MSC
        for card in "346":
            snvec = c.SN_VECS[int(card)]  # SN vector of the current card
            members = []  # MSC members for the current card
            for ord_ in range(1, len(snvec)):
                sn = snvec[ord_]  # SN of the current set
                pf = fromPFStr(dctSC[sn]["PF"])  # PF of the current set
                # Add the current set to MSC members if it is Kh-related
                #   or K-related (in case of hexachord) to the nexus set.
                rels = setComplexRelations(pf, nexusPF)
                if rels == 2 or (int(card) == 6 and rels >= 1
                                 and isSubset(nexusPF, pf)):
                    # MA
                    matts = dctSC[sn]["MA"]
                    # Degrees of symmetry
                    degrees = dctSC[sn]["symmetry"]
                    # Inclusion count
                    count = dctIncl[nexusSN][card][ord_-1]
                    # Z-corr
                    zcorr = dctSC[sn]["Z-corr"]
                    # dict for this member
                    member = dict(SN=sn, MA=matts,
                                  symmetry=degrees, inclusion=count)
                    member["Z-corr"] = zcorr
                    members.append(member)
            profile[card] = members
        dct[nexusSN] = profile
    return dict(dct)


def makePFLists():
    """
    Return a list of PFs for the cardinalities between 3 and 9 inclusive.

    :return: PFs for the cardinalities between 3 and 9 inclusive, sorted into
        nested lists based on their cardinalities.
    :rtype: list
    """
    pfs = []
    for card in range(3, 10):
        if card <= 6:
            s1 = sorted({tuple(primeForm(s)) for s in combinations(c.PCS, card)})
            pfs.append(s1)
        else:
            # PFs of card 7, 8, 9 are complements of card 5, 4, 3 respectively
            s2 = [tuple(primeForm(complement(s))) for s in pfs[(12-card)-3]]
            pfs.append(s2)
    return pfs


def makeCatalog():
    """
    Make a JSON file for database queries.
    """
    # Create a list of PFs for the cardinalities between 3 and 9 inclusive.
    pfs = makePFLists()

    # Make dictionaries for each query category
    dctName = mapPFToSN(pfs)  # dict for PF to SN conversion
    dctSC = mapSC(pfs, dctName)  # dict for SC profiles
    dctIncl = mapInclusionTable(dctSC)  # dict for inclusion table
    dctMSC = mapMSC(dctSC, dctIncl)  # dict for MSC

    # Make a master dictionary for JSON catalog
    data = {
        "PFToSN": dctName,
        "SC": dctSC,
        "inclusionTable": dctIncl,
        "MSC": dctMSC,
    }

    # Write the master dict to a JSON file
    with open("catalog.json", "w") as outfile:
        json.dump(data, fp=outfile, indent=4, sort_keys=True)


if __name__ == "__main__":
    answer = input("\nGenerate a new catalog file (Y/N)?\n>>> ")
    if answer.lower()[0] == "y":
        makeCatalog()
        print("\nDone!")
