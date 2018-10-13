# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""
Fake decorator to compute partial correlation functions.

It uses filters internally.
"""

import logging

from .helpers import filter_species

_log = logging.getLogger(__name__)


class Partial(object):

    def __init__(self, corr_cls, species, *args, **kwargs):
        """The first positional argument must be the trajectory instance."""
        # Instantiate correlation objects
        # with args passed upon construction
        self.partial = {}
        self.nbodies = corr_cls.nbodies
        self.species = species

        if self.nbodies == 1:
            for i in range(len(self.species)):
                isp = self.species[i]
                self.partial[isp] = corr_cls(*args, **kwargs)
                self.partial[isp].add_filter(filter_species, isp)
                self.partial[isp].tag = str(isp)
                self.partial[isp].tag_description = 'species %s' % isp

        elif self.nbodies == 2:
            for i in range(len(self.species)):
                for j in range(len(self.species)):
                    if j < i:
                        continue
                    isp = self.species[i]
                    jsp = self.species[j]
                    self.partial[(isp, jsp)] = corr_cls(*args, **kwargs)
                    self.partial[(isp, jsp)].add_filter(filter_species, isp)
                    # Slight optimization: avoid filtering twice when isp==jsp
                    if isp != jsp:
                        self.partial[(isp, jsp)].add_filter(filter_species, jsp)
                    self.partial[(isp, jsp)].tag = '%s-%s' % (isp, jsp)
                    self.partial[(isp, jsp)].tag_description = 'species pair %s-%s' % (isp, jsp)

    def need_update(self):
        need = False
        for partial in self.partial.values():
            if partial.need_update():
                need = True
                break
        return need

    def compute(self):
        if self.nbodies == 1:
            for i in range(len(self.species)):
                isp = self.species[i]
                self.partial[isp].compute()

        elif self.nbodies == 2:
            for i in range(len(self.species)):
                for j in range(len(self.species)):
                    isp = self.species[i]
                    jsp = self.species[j]
                    if j < i:
                        continue
                    self.partial[(isp, jsp)].compute()

    def do(self, update=False):
        if update and not self.need_update():
            return

        self.compute()

        for partial in self.partial.values():
            try:
                partial.analyze()
            except ImportError as e:
                _log.warn('Could not analyze due to missing modules, continuing...')
                _log.warn(e.message)

            partial.write()
