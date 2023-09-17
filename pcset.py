"""
pcset.py    Yota Kobayashi

This module includes a class pcset. The object of pcset class
can be instantiated as Pcset(pcs) where pcs is an iterable
(tuple, set, list, or Pcset object) and each pc is an int
(i.e., pcs is not a string).

The instance variable stores the pcs as a set representing an
unordered collection of unique elements. The class supports
set operational and relational methods. For the full list of the
methods, see the doc string of pcset class.

The class supports __iter__ and can return a generator object--
an iterable with the current pcs--therefore, method chaining is
possible as, for example:

a, b = {0,1,4}, {1,2}
s = Pcset(a)
s.transpose(4).union(b).normalForm()
"""

from itertools import combinations
from pcsets import constants as c

# __all__ = ["Pcset"]


class Pcset:  # TODO: will inherit MutableSet and define the abstract methods.
    """
    List of methods:

    BASIC
        __init__(pcs)
        __iter__()
        __str__()
        __len__()
        clone()

    ACCESSOR
        getSet()

    SET MEMBERSHIP
        union(pcs)
        difference(pcs)
        intersection(pcs)
        symmetricDifference(pcs)
        clear()

    SET OPERATION
        transpose(n)
        invert()
        invertXY(x, y)
        opT(n)
        opI()
        opTnI(n)
        opIxy(x, y)

    SET PROFILE
        icv()
        indexVector()
        normalForm()
        primeForm()
        transformationLevels()
        referentialCollections()

    SET ANALYSIS
        complement()
        modalComplements()
        subsets(n)
        transpositionalInvariants(n)
        inversionalInvariants(n)
        transpositionalSymmetry()
        inversionalSymmetry()

    SET RELATION
        pathSame(pcs)
        pathEmbed(pcs)
        pathCover(pcs)
        inclusion(pcs)
        complementation(pcs)
        icvsim(pcs)
        recrel(pcs)
    """

    # Basic methods -----------------------------------------------------------

    def __init__(self, pcs):
        """Constructor takes an iterable with pcs."""
        self.pcset = set(pcs)

    def __iter__(self):
        """Returns a generator object comprising the pcs."""
        for pc in self.pcset:
            yield pc

    def __repr__(self):
        """
        Return the official string representation of the object.

        >>> Pcset({0, 1, 4})
        PCSet({0, 1, 4})
        """
        return "Pcset({})".format(self.pcset)

    def __len__(self):
        """Returns the cardinality of the pc set."""
        return len(self.pcset)

    def copy(self):
        """Returns a deep copy of the current state of the object."""
        return Pcset(self.pcset)

    def clone(self):
        """Same as copy method--to be deprecated."""
        return self.copy()

    # Accessor methods --------------------------------------------------------

    def getSet(self):
        """
        Returns a shallow copy of the current pcset.

        To obtain a deep copy, use copy method.
        """
        return self.pcset

    # PC membership methods ---------------------------------------------------

    def union(self, pcs):
        """
        Add pcs to the current pc set. No duplicate counts.

        :param pcs: an iterable with pcs.
        """
        self.pcset |= set(pcs)
        return self

    def difference(self, pcs):
        """
        Mutates to a new set with elements in the current set and
        not in the input pcs.

        :param pcs: an iterable with pcs.
        """
        self.pcset -= set(pcs)
        return self

    def intersection(self, pcs):
        """
        Mutates to a new set with elements in both the current set and
        the input pcs.

        :param pcs: an iterable with pcs.
        """
        self.pcset &= set(pcs)
        return self

    def symmetricDifference(self, pcs):
        """
        Mutates to a new set with elements in either the current set
        or the input pcs but not both.

        :param pcs: an iterable with pcs.
        """
        self.pcset ^= set(pcs)
        return self

    def clear(self):
        """Removes all the pcs from the current set."""
        self.pcset.clear()
        return self

    # Set transformation methods ----------------------------------------------

    def transpose(self, n):
        """
        Transposes the current set by an ordered pc interval n.

        :param n: mod 12 integer for the transposition number.
        """
        self.pcset = {(pc + n) % 12 for pc in self.pcset}
        return self

    def invert(self):
        """
        Inverts the current set around pc 0.
        """
        self.pcset = {(12 - pc) % 12 for pc in self.pcset}
        return self

    def invertXY(self, x, y):
        """
        Inverts the current set around an axis of symmetry specified by pcs x and y:
        as a result, pcs x and y map onto each other.
        This is an operation equivalent to TnI where n = x + y. The advantage of
        using I(x,y) over TnI is that it allows for specifying the axis of symmetry.

        :param x: int that inverts x onto y.
        :param y: int that inverts y onto x.
        """
        self.opTnI((x + y) % 12)
        return self

    def opT(self, n):
        """
        Shorthand for transpose(n).
        """
        self.transpose(n)
        return self

    def opI(self):
        """
        Shorthand for invert().
        """
        self.invert()
        return self

    def opTnI(self, n):
        """
        Shorthand for TnI operation--inversion followed by transposition at n.
        """
        self.invert().transpose(n)
        return self

    def opIxy(self, x, y):
        """
        Shorthand for invertXY(x, y).
        """
        self.invertXY(x, y)
        return self

    # Set profile methods -----------------------------------------------------

    def icv(self):
        """
        A method to compute an interval-class vector (ICV) of the
        current set--a tabular representation for counting of unordered
        pc intervals (i.e., ics).

        :return: a list representing the ICV.
        """
        icv = [0] * 6
        pcList = sorted(self.pcset)
        while len(pcList) > 1:
            pc1 = pcList.pop()
            for pc in pcList:
                i = (pc - pc1) % 12  # Ordered interval class
                i = min(i, 12 - i)
                icv[i-1] += 1
        return icv

    def indexVector(self):
        """
        A method to compute an index vector of the input pc set.

        :return: a list representing the index vector.
        """
        pcs = list(self.pcset)
        ivec = [0] * 12
        # Index vector is computed using addition table.
        for i in pcs:
            for j in pcs:
                ivec[(i + j) % 12] += 1
        return ivec

    def normalForm(self):
        """
        A method to compute the normal form of the current set.

        :return: a list of pcs representing the normal form.
        """
        card = len(self.pcset)       # Set cardinality
        if card == 0:
            return []
        elif card == 12:
            return list(range(12))
        pcList = sorted(self.pcset)  # pc set to a sorted list
        rots = []                    # Array for rotations of pcList
        for r in range(card):        # Make all rotations of pcList
            rots.append(self.__rotate(pcList, r))
        k = card - 1
        while len(rots) > 1:
            if k == 0:  # If k=0 is reached, NF is the one with the smallest first pc.
                return rots[0]
            else:
                intervals = [(rot[k] - rot[0]) % 12 for rot in rots]
                smallest = min(intervals)  # Smallest intervals of each member of rots
                # Delete a set in rots if its interval > the smallest interval
                rots = [rots[i] for i in range(len(rots)) if intervals[i] == smallest]
                k -= 1
        return rots[0]

    def primeForm(self):
        """
        A method to compute the prime form of the current set.

        :return: a list of pcs representing the prime form.
        """
        card = len(self.pcset)
        if card == 0:
            return []
        elif card == 12:
            return list(range(12))
        s1 = self.normalForm()
        s2 = Pcset(s1).invert().normalForm()
        s1 = [(pc - s1[0]) % 12 for pc in s1]
        s2 = [(pc - s2[0]) % 12 for pc in s2]
        for i in range(len(s1) - 2, 0, -1):
            if s1[i] < s2[i]:
                return s1
            elif s1[i] > s2[i]:
                return s2
        return s1  # s1 = s2, a symmetrical set

    def transformationLevels(self):
        """
        Returns the Tn/TnI transformation level of the current set.

        :return: a dict with two nested lists, {"Tn": [], "TnI": []}, where
            the list Tn and list TnI represent the transformation levels of the
            current set. For transpositionally and inversionally symmetrical sets,
            there would be multiple entries in the lists.
        """
        nf = self.normalForm()
        pf = self.primeForm()
        return Pcset(pf).pathSame(nf)

    def referentialCollections(self):
        """
        A method to check the reference status to the modal collections,
        OCT, WT, HEX, and DT.

        :return: a dict with abbreviated keys for OCT, WT, and HEX as c+n where c is
            the collection's initial and n is the transposition level e.g., O0, W1, H3, etc.
            These keys take int of 3 possible values:
                0       no reference
                1       all but one pc are literal subset of the Tn level of the collection
                2       all the pcs are literal subset of the Tn level of the collection
            An additional key is for DT as D with a boolean value.
                True    if the current set is an abstract subset of DT
                False   otherwise
        """
        # Create referential collection dictionary with default val 0 to every key.
        refCols = {col: 0 for col in list(c.COL_DICT.keys())}
        # Check the subset status of the current set against OCT, WT, and HEX collections.
        for col in refCols:
            refCols[col] = self.__subsetStatus(c.COL_DICT[col])
        refCols["D"] = False  # Add abstract subset status for DT collection
        path = self.pathEmbed({0, 1, 3, 5, 6, 8, 10})  # Operational paths to DT
        if len(path["Tn"]) + len(path["TnI"]) > 0:
            refCols["D"] = True
        return refCols

    # Set analysis methods ----------------------------------------------------

    def complement(self):
        """
        Returns a set of the complement of the current set.
        """
        return c.TT - self.pcset

    def modalComplements(self):
        """
        Returns the complements to each of the modal collections OCT, WT, and HEX.

        :return: a dict with abbreviated keys for OCT, WT, and HEX as c+n where c is
            the collection's initial and n is the transposition level e.g., O0, W1, H3, etc.
            Each key has a value for the complementing pcs to the modal collection.
        """
        modalComps = {col: {} for col in list(c.COL_DICT.keys())}
        for col in modalComps:
            modalComps[col] = c.COL_DICT[col] - self.pcset
        return modalComps

    def subsets(self, n):
        """
        Returns all the subsets of cardinality n of the current set.

        :param n: the cardinality of the subsets.
        :return: a list of sets--the subsets of cardinality n.
        """
        if n >= len(self.pcset):
            return None
        return [set(s) for s in combinations(self.pcset, n)]

    def transpositionalInvariants(self, n):
        """
        Returns the invariants (common tones) held under transposition
        of a pc set by interval n. The number of invariants is equal to the
        number of times the interval n occurs in the set. The digits of
        ICV represent intervals n and 12-n. Double the count of ic 6.

        :param n: an int for the ordered pc interval of transposition.
        :return: a set of invariant pcs.
        """
        if n == 0:  # ICV has no entry for ic0
            return set(self.pcset)
        copy = self.clone()
        n = min(12-n, n)
        return self.pcset & copy.transpose(n).getSet()

    def inversionalInvariants(self, n):
        """
        Returns the invariants (common tones) held under inverting
        a pc set by interval n. The number of invariants can be checked
        in the index vector of the set.

        :param n: int for the index number.
        :return: a set of invariant pcs.
        """
        copy = self.clone()
        return self.pcset & copy.opTnI(n).getSet()

    def transpositinalSymmetry(self):
        """
        Returns the degree of transpositional symmetry.

        :return: a list of ordered pc intervals--the input set maps entirely
            onto itself when transposing at the transpositional level(s).
            The number is at least 1, because T0 is an identity operator.
        """
        icv = self.icv()
        icv[5] *= 2
        sym = {0}  # T0 is identity operation
        size = len(self.pcset)
        for i in range(6):
            if icv[i] == size:
                sym |= {(i + 1), 12 - (i + 1)}
        return sorted(sym)

    def inversionalSymmetry(self):
        """
        Returns the degree of inversional symmetry.

        :return: a list of ordered pc intervals--the input set maps entirely
            onto itself when inverting with the index number(s).
            None is output when the set is not inversionally symmetrical.
        """
        ivec = self.indexVector()
        return [i for i in range(12) if ivec[i] == len(self.pcset)]

    # Set relation methods ----------------------------------------------------

    def pathSame(self, pcs):
        """
        A method to find the operational paths to transform the current set
        into the input set, that is, to make both sets comprise the same pcs.
        If set1 can be transformed into set2 through Tn or TnI, this method
        will find the possible values for n.

        :param pcs: an iterable with pcs.
        :return: a dict with two nested lists, {"Tn": [], "TnI": []}:
            The list Tn comprises possible values for n where Tn(set1) == set2.
            The list TnI comprises possible values for n where TnI(set1) == set2.
        """
        path = {"Tn": [], "TnI": []}
        pcs = set(pcs)
        if self.primeForm() != Pcset(pcs).primeForm():
            return path
        for n in range(12):
            tn, tni = Pcset(self.pcset), Pcset(self.pcset)
            if tn.opT(n).getSet() == pcs:
                path["Tn"].append(n)
            if tni.opTnI(n).getSet() == pcs:
                path["TnI"].append(n)
        return path

    def pathEmbed(self, pcs):
        """
        A method to find the operational path to make the current set
        a literal subset of the input set.

        :param pcs: an iterable with pcs.
        :return: a dict with two nested lists, {"Tn": [], "TnI": []}:
            The list Tn comprises possible values for n where Tn(set1) < set2.
            The list TnI comprises possible values for n where TnI(set1) < set2.
        """
        path = {"Tn": [], "TnI": []}
        pcs = set(pcs)
        if len(self.pcset) >= len(pcs):
            return path
        for n in range(12):
            tn, tni = Pcset(self.pcset), Pcset(self.pcset)
            if tn.opT(n).getSet() < pcs:
                path["Tn"].append(n)
            if tni.opTnI(n).getSet() < pcs:
                path["TnI"].append(n)
        return path

    def pathCover(self, pcs):
        """
        A method to find the operational path to make the current set
        a literal superset of the input set.

        :param pcs: an iterable with pcs.
        :return: a dict with two nested lists, {"Tn": [], "TnI": []}:
            The list Tn comprises possible values for n where Tn(set1) > set2.
            The list TnI comprises possible values for n where TnI(set1) > set2.
        """
        path = {"Tn": [], "TnI": []}
        pcs = set(pcs)
        if len(self.pcset) <= len(pcs):
            return path
        for n in range(12):
            tn, tni = Pcset(self.pcset), Pcset(self.pcset)
            if tn.opT(n).getSet() > pcs:
                path["Tn"].append(n)
            if tni.opTnI(n).getSet() > pcs:
                path["TnI"].append(n)
        return path

    def inclusion(self, pcs):
        """
        A method to find all the possible sets included in the current set
        that are the members of the set class of pcs.

        :param pcs: an iterable with pcs representing the target set class.
        :return: a list of tuples, each tuple contains two sets (target, diff):
            target: a literal subset of the current set of the target set class
            diff: a set of difference pcs--pcs in the current set but not in the
                target set.
        """
        target = set(Pcset(pcs).primeForm())
        card = len(target)
        targets = []
        # Pretest--target must be smaller than the current set
        if len(self.pcset) - card <= 0:
            return targets
        for s in list(self.subsets(card)):
            if set(Pcset(s).primeForm()) == target:
                targets.append((s, self.pcset - s))
        return targets

    def complementation(self, pcs):
        """
        A method to find all the possible sets to complement the current
        set to form the members of the set class of pcs.

        :param pcs: an iterable with pcs representing the target set class.
        :return: a list of tuples, each tuple contains two sets (target, diff):
            target: a complemented set that is a union of the current set and
                the set of complementing pcs s2
            diff: a set of difference pcs--pcs in the target set but not in the
                current set, that is, the complementing pcs for the current set
                to form the target set class.
        """
        target = set(Pcset(pcs).primeForm())
        targets = []
        gap = len(target) - len(self.pcset)
        if gap <= 0:  # Pretest--target must be larger than the current set
            return targets
        candidates = [set(s) for s in combinations(c.TT - self.pcset, gap)]
        for diff in candidates:
            s = Pcset(self.pcset | diff)
            if set(s.primeForm()) == target:
                targets.append((set(s), diff))
        return targets

    def icvsim(self, pcs, raw=False):
        """
        Return a set-class similarity value measured by Isaacson's ICVSIM
        (Isaacson, 1990).

        .. note::
            The original ICVSIM outputs a similarity value in the range [0.00,
            3.578...] where a value of 0.00 indicates maximum similarity with
            respect to interval-class content. When raw=False, this method
            returns the raw similarity value; when raw=True, the similarity
            value is "flipped" and normalized so that it returns a float [0.00,
            1.00] where 1.00 indicates the maximum similarity.

        >>> s = {0, 1, 2}  # 3-1
        >>> Pcset(s).icvsim({0,1,2,3,6}, raw=True)    # ICVSIM(3-1, 5-4)
        0.3726779962499649
        >>> Pcset(s).icvsim({0,1,2,5,6,7}, raw=True)  # ICVSIM(3-1, 6-Z6)
        1.0
        >>> Pcset(s).icvsim({0,1,4,6,9}, raw=True)    # ICVSIM(3-1, 5-32)
        1.343709624716425
        >>> Pcset(s).icvsim({0,3,6,9}, raw=True)      # ICVSIM(3-1, 4-28)
        1.9790570145063195
        >>> Pcset({0,2,4,6,8,10}).icvsim({0,1,3,4,6,7,9,10})  # ICVSIM(6-35, 8-28)
        0.0
        >>> Pcset({0,1,3,4,8}).icvsim({0,3,4,5,8})  # ICVSIM(5-Z17, 5-Z37)
        1.0

        :param pcs: an iterable with pcs.
        :param raw: If True, the returned ICVSIM value is in the range [0.0,
            3.578...], where 0.0 indicates the maximum similarity; if False
            (default), the value is flipped and normalized so that it is in
            the range [0.0, 1.0] where 1.0 indicates the maximum similarity.
        :return: ICVSIM value.
        :rtype: float
        """
        from pcsets.relation import icvsim  # FIXME: import will be at top.
        return icvsim(self.pcset, pcs, raw=raw)

    # Private methods ---------------------------------------------------------

    def __rotate(self, lst, n):
        """A helper function to rotate a list by n items."""
        return lst[n:] + lst[:n]

    def __subsetStatus(self, pcs):
        """
        A helper method that returns the literal subset status
        of the current set against the input set.

        :param pcs: an iterable with pcs.
        :return:
            0: more than one pc are not the elements of the input set
            1: all but one pc are the elements of the input set
            2: literal subset in, or the same as, the input set
        """
        pcs = set(pcs)
        if self.pcset > pcs:  # Pretest
            return 0
        n = len(self.pcset - pcs)
        if n == 0:
            return 2
        elif n == 1:
            return 1
        else:
            return 0
