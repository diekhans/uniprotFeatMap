"""
check and dirty file generation dependency checking
"""
from os import path as osp

class DependsError(Exception):
    pass

def getDoneFile(target):
    "get path to done file corresponding to target"
    return osp.normpath(target) + ".done"

def _normalizeDepends(depends):
    "make list if str or None"
    if depends is None:
        depends = ()
    elif isinstance(depends, str):
        depends = (depends,)
    return depends

def _checkDepends(depends):
    for depend in depends:
        if not osp.exists(depend):
            raise DependsError(f"dependency does not exist: {depend}")

def _isOutOfDate(doneFile, depends):
    modTime = osp.getmtime(doneFile)
    for depend in depends:
        if osp.getmtime(depend) > modTime:
            return True
    return False

def _checkDoneDepends(doneDepends):
    for doneDepend in doneDepends:
        depend = getDoneFile(doneDepend)
        if not osp.exists(depend):
            raise DependsError(f"dependency does not exist: {depend}")

def _isDoneOutOfDate(doneFile, doneDepends):
    return _isOutOfDate(doneFile, [getDoneFile(d) for d in doneDepends])

class runIfNotDone:
    """Context manager for a Cheap make-like hack. Checks if a flag file to
    see if something is done. The with context returns a True if target needs
    build, false otherwise.  Checks if target exists. If there are
    dependencies, the check that target.done is newer.  depends are check for
    actually files, doneDepends are check for done files

    If the body succeeds, a file target.done is created. If target is a
    directory, create a file target.done beside it, not inside of it.

    """

    def __init__(self, target, *, depends=None, doneDepends=None):
        "depends can be a string filename or list"
        depends = _normalizeDepends(depends)
        _checkDepends(depends)
        doneDepends = _normalizeDepends(doneDepends)
        _checkDoneDepends(doneDepends)

        self.doneFile = getDoneFile(target)
        self.outOfDate = ((not osp.exists(self.doneFile)) or
                          _isOutOfDate(self.doneFile, depends) or
                          _isDoneOutOfDate(self.doneFile, doneDepends))

    def __enter__(self):
        return self.outOfDate

    def __exit__(self, exc_type, exc_value, traceback):
        if (exc_type is None) and self.outOfDate:
            open(self.doneFile, 'w').close()

class runIfOutOfDate:
    """Context manager for a Cheap make-like hack. Checks if an atomically
    create does not exist or is out-of-date with optional dependencies
    or done flag files for dependenices.
    The with context returns a True if target needs built, false otherwise.
    """

    def __init__(self, target, *, depends=None, doneDepends=None):
        "depends can be a string filename or list"
        if isinstance(depends, str):
            depends = [depends]
        if isinstance(doneDepends, str):
            doneDepends = [doneDepends]

        self.outOfDate = ((not osp.exists(target)) or
                          ((depends is not None) and _isOutOfDate(target, depends)) or
                          ((doneDepends is not None) and _isDoneOutOfDate(target, doneDepends)))

    def __enter__(self):
        return self.outOfDate

    def __exit__(self, exc_type, exc_value, traceback):
        pass
