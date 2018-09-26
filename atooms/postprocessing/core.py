from ._version import __version__ as __bare_version

try:
    from ._commit import __commit__, __date__
    __version__ = '%s+%s (%s)' % (__bare_version, __commit__, __date__)
except ImportError:
    __commit__ = ""
    __date__ = ""
    __version__ = __bare_version
