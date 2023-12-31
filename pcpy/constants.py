# Package-wide constants

# Pitch-classes
PCS = set(range(12))

# Collection aggregates
TT = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
OCT0 = {0, 1, 3, 4, 6, 7, 9, 10}
OCT1 = {1, 2, 4, 5, 7, 8, 10, 11}
OCT2 = {2, 3, 5, 6, 8, 9, 11, 0}
WT0 = {0, 2, 4, 6, 8, 10}
WT1 = {1, 3, 5, 7, 9, 11}
HEX0 = {0, 1, 4, 5, 8, 9}
HEX1 = {1, 2, 5, 6, 9, 10}
HEX2 = {2, 3, 6, 7, 10, 11}
HEX3 = {3, 4, 7, 8, 11, 0}

# Abbreviated names of the referential collections and transposition levels
REF_COLS = ["O0", "O1", "O2", "W0", "W1", "H0", "H1", "H2", "H3"]

# Collection dictionary
COL_DICT = {"O0": OCT0, "O1": OCT1, "O2": OCT2,
            "W0": WT0, "W1": WT1,
            "H0": HEX0, "H1": HEX1, "H2": HEX2, "H3": HEX3}

# Ordinal numbers in the order of sequentially output prime forms
ORDINAL_NUMS = [
    ["dummy"],
    # Cardinality 1
    ["dummy", "1"],
    # Cardinality 2
    ["dummy", "1", "2", "3", "4", "5", "6"],
    # Cardinality 3
    ["dummy", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"],
    # Cardinality 4
    ["dummy", "1", "2", "4", "5", "6", "3", "11", "13", "Z29", "7", "Z15",
     "18", "19", "8", "16", "20", "9", "10", "12", "14", "21", "22", "24",
     "23", "27", "25", "17", "26", "28"],
    # Cardinality 5
    ["dummy", "1", "2", "4", "5", "3", "9", "Z36", "13", "6", "14", "Z38",
     "7", "15", "10", "16", "Z17", "Z12", "24", "27", "19", "29", "31",
     "Z18", "21", "30", "32", "22", "20", "8", "11", "23", "25", "28", "26",
     "33", "34", "35", "Z37"],
    # Cardinality 6
    ["dummy", "1", "2", "Z36", "Z37", "Z3", "9", "Z40", "5", "Z41", "Z42",
     "Z38", "Z4", "Z11", "15", "Z12", "22", "Z46", "Z17", "Z47", "Z6",
     "Z43", "Z44", "18", "Z48", "7", "Z10", "14", "Z13", "Z24", "27", "Z19",
     "Z49", "Z25", "Z28", "Z26", "34", "30", "16", "31", "20", "Z50", "8",
     "Z39", "21", "Z45", "Z23", "33", "Z29", "32", "35"],
    # Cardinality 7
    ["dummy", "1", "2", "4", "5", "3", "9", "Z36", "13", "6", "14", "Z38",
     "7", "15", "10", "16", "Z17", "Z12", "24", "27", "19", "29", "31",
     "Z18", "21", "30", "32", "22", "20", "8", "11", "23", "25", "28", "26",
     "33", "34", "35", "Z37"],
    # Cardinality 8
    ["dummy", "1", "2", "4", "5", "6", "3", "11", "13", "Z29", "7", "Z15",
     "18", "19", "8", "16", "20", "9", "10", "12", "14", "21", "22", "24",
     "23", "27", "25", "17", "26", "28"],
    # Cardinality 9
    ["dummy", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"],
    # Cardinality 10
    ["dummy", "1", "2", "3", "4", "5", "6"],
    # Cardinality 11
    ["dummy", "1"],
    # Cardinality 12
    ["dummy", "1"]
]

# Set-name vectors
SN_VECS = [
    ["dummy"],
    ["dummy"],
    ["dummy"],
    ["dummy", "3-1", "3-2", "3-3", "3-4", "3-5", "3-6", "3-7", "3-8", "3-9",
     "3-10", "3-11", "3-12"],
    ["dummy", "4-1", "4-2", "4-3", "4-4", "4-5", "4-6", "4-7", "4-8", "4-9",
     "4-10", "4-11", "4-12", "4-13", "4-14", "4-Z15", "4-16", "4-17",
     "4-18", "4-19", "4-20", "4-21", "4-22", "4-23", "4-24", "4-25", "4-26",
     "4-27", "4-28", "4-Z29"],
    ["dummy", "5-1", "5-2", "5-3", "5-4", "5-5", "5-6", "5-7", "5-8", "5-9",
     "5-10", "5-11", "5-Z12", "5-13", "5-14", "5-15", "5-16", "5-Z17",
     "5-Z18", "5-19", "5-20", "5-21", "5-22", "5-23", "5-24", "5-25",
     "5-26", "5-27", "5-28", "5-29", "5-30", "5-31", "5-32", "5-33", "5-34",
     "5-35", "5-Z36", "5-Z37", "5-Z38"],
    ["dummy", "6-1", "6-2", "6-Z3", "6-Z4", "6-5", "6-Z6", "6-7", "6-8", "6-9",
     "6-Z10", "6-Z11", "6-Z12", "6-Z13", "6-14", "6-15", "6-16", "6-Z17",
     "6-18", "6-Z19", "6-20", "6-21", "6-22", "6-Z23", "6-Z24", "6-Z25",
     "6-Z26", "6-27", "6-Z28", "6-Z29", "6-30", "6-31", "6-32", "6-33",
     "6-34", "6-35", "6-Z36", "6-Z37", "6-Z38", "6-Z39", "6-Z40", "6-Z41",
     "6-Z42", "6-Z43", "6-Z44", "6-Z45", "6-Z46", "6-Z47", "6-Z48", "6-Z49",
     "6-Z50"],
    ["dummy", "7-1", "7-2", "7-3", "7-4", "7-5", "7-6", "7-7", "7-8", "7-9",
     "7-10", "7-11", "7-Z12", "7-13", "7-14", "7-15", "7-16", "7-Z17",
     "7-Z18", "7-19", "7-20", "7-21", "7-22", "7-23", "7-24", "7-25",
     "7-26", "7-27", "7-28", "7-29", "7-30", "7-31", "7-32", "7-33", "7-34",
     "7-35", "7-Z36", "7-Z37", "7-Z38"],
    ["dummy", "8-1", "8-2", "8-3", "8-4", "8-5", "8-6", "8-7", "8-8", "8-9",
     "8-10", "8-11", "8-12", "8-13", "8-14", "8-Z15", "8-16", "8-17",
     "8-18", "8-19", "8-20", "8-21", "8-22", "8-23", "8-24", "8-25", "8-26",
     "8-27", "8-28", "8-Z29"],
    ["dummy", "9-1", "9-2", "9-3", "9-4", "9-5", "9-6", "9-7", "9-8", "9-9",
     "9-10", "9-11", "9-12"]
]

# Set names and prime forms of the MSC nexus sets
NEXUS_SETS = [
    ("5-10", (0, 1, 3, 4, 6)),
    ("5-16", (0, 1, 3, 4, 7)),
    ("5-19", (0, 1, 3, 6, 7)),
    ("5-21", (0, 1, 4, 5, 8)),
    ("5-25", (0, 2, 3, 5, 8)),
    ("5-28", (0, 2, 3, 6, 8)),
    ("5-31", (0, 1, 3, 6, 9)),
    ("5-32", (0, 1, 4, 6, 9)),
    ("5-33", (0, 2, 4, 6, 8))
]
