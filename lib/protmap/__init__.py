"""
Library of common functions and other definitions
"""

import sys
from collections import defaultdict

def prMsg(msg):
    print(msg, file=sys.stderr, flush=True)

def dropVersion(ident):
    return ident.split('.')[0]

def buildDfUniqueIndex(df, column):
    """build a dictionary to serve as a unique index for a DataFrame column.
    Used because this can be faster and you can have multiple indexes.
    """
    idx = {}
    for _, row in df.iterrows():
        assert row[column] not in idx
        idx[row[column]] = row
    return idx

def buildDfMultiIndex(df, column):
    """build a dictionary to serve as a multi-index for a DataFrame column
    Used because this can be faster and you can have multiple and non-unique indexes.
    """
    idx = defaultdict(list)
    for _, row in df.iterrows():
        idx[row[column]].append(row)
    idx.default_factory = None
    return idx
