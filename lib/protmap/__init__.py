"""
Library of common functions and other definitions
"""

import sys

def prMsg(msg):
    print(msg, file=sys.stderr, flush=True)

def dropVersion(ident):
    return ident.split('.')[0]
