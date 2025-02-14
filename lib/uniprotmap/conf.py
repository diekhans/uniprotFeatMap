"""
configuration for decorator build.
"""

# must match runBlast
blastDir = "/cluster/bin/blast/x86_64/blast-2.2.16/bin"

paraHost = None

##
# initial filtering of protein/transcript alignments
##

# want to allow for alignments to other isoforms, so let a lot pass
protTransAlnMinCover = 0.0

# stricter filter keeps if from being too far off
protTransAlnMinId = 0.80
