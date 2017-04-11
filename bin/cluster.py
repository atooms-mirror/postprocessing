#!/usr/bin/env python

"""
Cluster analysis
"""

from postprocessing.cluster import main

if __name__ == '__main__':

    import argh
    argh.dispatch_command(main)
    
