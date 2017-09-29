#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
if state.method_style=='basic': get_angle_torsion_distn_from_pentamer()
elif state.method_style=='advanced': dextran_backmapper()
else: raise Exception('invalid method_style %s'%method_style)
