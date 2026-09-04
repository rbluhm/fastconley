#!/usr/bin/env python3
"""Regenerate stata/upstream/Conley.mata from stata/src/fastconley.mata.

The engine functions are copied with a reghdfe_conley_ prefix and the
reghdfe_vce_conley() front-end (modelled on reghdfe_vce_dkraay) is prepended.
Run from the repository root after any change to the Stata engine, then
refresh reghdfe-vce-conley.patch (see README.md)."""
import re, sys
src = open('stata/src/fastconley.mata').read()
def grab(start_marker):
    i = src.index(start_marker)
    hdr = src.rfind('// ---------------------------------------------------------------------------\n', 0, i)
    j = src.index('\n}\n', i) + 3
    return src[hdr:j]
parts = [grab('real colvector fastconley_group('), grab('void fastconley_aggregate('),
         grab('real matrix fastconley_spatial_meat(')]
i = src.index('// One cell pair (rows ia..ib against ja..jb), tiled.')
parts.append(src[i:src.index('\n}\n', i) + 3])
parts.append(grab('real matrix fastconley_serial_meat('))
engine = '\n'.join(parts).replace('fastconley_', 'reghdfe_conley_')
front = open('stata/upstream/Conley.mata').read()
cut = front.index('// ---------------------------------------------------------------------------\n// Engine (derived from')
front = front[:cut]
out = front + '''// ---------------------------------------------------------------------------
// Engine (derived from fastconley.mata of the fastconley Stata package by
// make_conley_mata.py; do not edit here). Plain Mata, no reghdfe types.
// ---------------------------------------------------------------------------
mata:

''' + engine + '''
end
'''
open('stata/upstream/Conley.mata', 'w').write(out)
print('Conley.mata regenerated:', out.count('\n'), 'lines')
