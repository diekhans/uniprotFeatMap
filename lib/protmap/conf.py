"""
configuration for decorator build.
"""

# must match runBlast
blastDir = "/cluster/bin/blast/x86_64/blast-2.2.16/bin"

paraHost = "ku"

##
# initial filtering of alignments
##

# want to allow for alignments to other isoforms, so let a lot pass
protAlnMinCover = 0.25

# stricter filter keeps if from being too far off
protAlnMinId = 0.80
