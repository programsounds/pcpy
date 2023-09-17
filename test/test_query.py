# query_test.py
# Testing script for querying the catalog dict variable defined in the query
#  module.

from pprint import pprint
from pcsets.query import catalog

print(catalog["SC"]["3-1"])
print(catalog["inclusionTable"]["4-12"]["3"])
pprint(catalog["MSC"]["5-10"])
