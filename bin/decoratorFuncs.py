"""
Library of common functions and other definitions for building domain decorators.
"""

from os import path as osp

def _isOutOfDate(doneFile, dependencies):
    modTime = osp.getmtime(doneFile)
    for depend in dependencies:
        if osp.getmtime(depend) > modTime:
            return True
    return False

def getDoneFile(target):
    "get path to done file corresponding to target"
    return osp.normpath(target) + ".done"

class runIfNotDone:
    """Cheap make-like hack. checks if a file target.done exists.  If
    dependencies are listed, the check that target.done is newer.  If target
    is a directory, create a file target.done beside it, not inside of it.
    """

    def _isOutOfDate(self, doneFile, dependencies):
        if isinstance(dependencies, str):
            dependencies = [dependencies]
        modTime = osp.getmtime(doneFile)
        for depend in dependencies:
            if osp.getmtime(depend) > modTime:
                return True
        return False

    def __init__(self, target, dependencies=None):
        "dependencies can be a string filename or list"
        self.doneFile = getDoneFile(target)
        self.outOfDate = ((not osp.exists(self.doneFile)) or
                          ((dependencies is not None) and _isOutOfDate(dependencies)))

    def __enter__(self):
        return self if self.outOfDate else None

    def __exit__(self, exc_type, exc_value, traceback):
        if (exc_type is None) and self.outOfDate:
            open(self.doneFile, 'w').close()
