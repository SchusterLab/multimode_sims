# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 09:54:15 2021

@author: ziqian
"""

from slab.circuits import *
from slab.circuits.mp_components import *

import os, time, sys
import subprocess
from time import sleep
import numpy as np
from matplotlib import pyplot as plt
import csv
import matplotlib as mpl
import datetime as dt

plt.close('all')

# CHIPNAME = "UC"  # dt.datetime.today().strftime("%Y%m%d") #+ "ParaAmp"
today = time.strftime('%Y%m%d')
author = ''

# Only change these two variables!!
pad_layer = 1
junction_layer = 1

# I don't know why it works...
if pad_layer == 1 and junction_layer == 0:
    draw_layer1 = 1
    draw_layer2 = 0
    draw_layer3 = 1

if junction_layer == 1 and pad_layer == 0:
    draw_layer1 = 0
    draw_layer2 = 1
    draw_layer3 = 1

if pad_layer == 1 and junction_layer == 1:
    draw_layer1 = 1
    draw_layer2 = 1
    draw_layer3 = 1

show_structure = 0
show_wafer = 1
two_layer = 1
perf = False
solid = 0
square = 0
etching = True  # False if a bare wafer is used
open_dwgviewer = False
open_klayout = True

if draw_layer1 == 1:
    MaskName = ""
else:
    MaskName = " "  # MaskName not needed for junction mask

if draw_layer1 == 1 and draw_layer2 == 0:
    two_layer = 0
elif draw_layer1 == 0 and draw_layer3 == 1:
    two_layer = 0
else:
    two_layer = 1

if junction_layer == 1:
    two_layer = 1
# Don't change anything above!

### CPW Parameters
cpw_length = 10
cpw_pinw = 10.0
cpw_rad = 50
taperl = 50
eps_eff = (1. + 10.4) / 2.
cpw_gapw = calculate_gap_width(eps_eff, 50, cpw_pinw)
flux_pinw = 3.5
flux_gapw = calculate_gap_width(eps_eff, 50, flux_pinw) + 0.21

### JPA capacitor parameters
cap_botm_y = 40
cap_botm_x = 140
cap_top_y = 120
cap_top_x = 80
cap_sep = 11
con_pin_l = 50.

junc_pad_y = 6.5

print(cpw_gapw)

chip_names = list()
chanstarts = list()

### Close DWG Viewer & KLayout
if open_dwgviewer:
    subprocess.Popen(r'taskkill /F /im "dwgviewr.exe"')
if open_klayout:
    subprocess.Popen(r'taskkill /F /im "klayout_app.exe"')


def set_mask_init():
    ### Setting defaults
    d = ChipDefaults()
    d.Q = 1000
    d.radius = 50
    d.segments = 6
    d.pinw_rsn = 2.  # this is the resonator pinwitdth that we are goign to use.
    d.gapw_rsn = 8.5

    d.pinw = 1.5  # d.pinw
    d.gapw = 1.
    d.center_gapw = 1
    ### Now calculate impedance
    d.imp_rsn = 80.  # calculate_impedance(d.pinw_rsn,d.gapw_rsn,d.eps_eff)
    d.solid = True
    return d


def chipInit(c, defaults):
    """
    This makes the launch pads on the chip. Input is an object c, which is the chip.
    From the launch pads, we can make connections on the chip.
    """
    # The following creates 8 launch pads on the chip. There are 3 pads per side, two
    # in each of the corners and then one in the middle of each side.

    setattr(c, 's0',
            Structure(c, start=(c.size[0] / 2 - 1900 - 800-100+100+100, c.size[1] / 2 - 800 - 5.), direction=-90, defaults=defaults,
                      layer='0'))
    setattr(c, 's1',
            Structure(c, start=(c.size[0] / 2 - 1900 - 800-100+100+100, c.size[1] / 2 + 1000 + 450), direction=-90, defaults=defaults,
                      layer='0'))
    setattr(c, 's2',
            Structure(c, start=(c.size[0] / 2 - 1900, c.size[1] - 600+100-100), direction=270, defaults=defaults, layer='0'))
    setattr(c, 's3',
            Structure(c, start=(c.size[0] / 2 - 9, c.size[1] - 600+100-100), direction=270, defaults=defaults, layer='0'))
    setattr(c, 's4',
            Structure(c, start=(c.size[0] / 2 + 1900+200, c.size[1] - 600+100-100), direction=180, defaults=defaults, layer='0'))
    setattr(c, 's5',
            Structure(c, start=(c.size[0] / 2 + 1900 + 800+100-100-100, c.size[1] / 2 + 1000 + 450), direction=-90, defaults=defaults,
                      layer='0'))
    setattr(c, 's6',
            Structure(c, start=(c.size[0] / 2 + 1900 + 800+100-100-100, c.size[1] / 2 - 800), direction=-90, defaults=defaults,
                      layer='0'))
    setattr(c, 's7', Structure(c, start=(c.size[0] / 2 + 1900+200, 800-100+100), direction=180, defaults=defaults, layer='0'))
    setattr(c, 's8', Structure(c, start=(c.size[0] / 2, 800-100+100), direction=90, defaults=defaults, layer='0'))
    setattr(c, 's9', Structure(c, start=(c.size[0] / 2 - 1900-200, 800-100+100), direction=0, defaults=defaults, layer='0'))


def draw_launchers(c, d, exclude=[]):
    chipInit(c, defaults=d)

    for k in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:

        if k in exclude:
            pass
        else:
            if k==k:
                Launcher(vars(vars()['c'])['s%d' % k], pinw=cpw_pinw, gapw=cpw_gapw, launcher_pinw=150, launcher_gapw=75)
            else:
                Launcher(vars(vars()['c'])['s%d' % k], pinw=cpw_pinw, gapw=cpw_gapw, launcher_pinw=300,
                         launcher_gapw=150)


def cover_launchers(c, d, exclude=[], h=500, w=400):
    """
    Cover the launchers with a square so that we can wirebond to the pads on the chip.
    """

    # taper_length = 250
    # taper_to_width = 2 * 50 + 20
    #
    # s = Structure(c, start=c.top_midpt, direction=270, defaults=d)
    #
    # if not (1 in exclude):  # middle, top
    #     lo_left = (c.top_midpt[0] - w / 2., c.top_midpt[1] - (h - taper_length))
    #     lo_right = (c.top_midpt[0] + w / 2., c.top_midpt[1] - (h - taper_length))
    #     tp_left = (c.top_midpt[0] - w / 2., c.top_midpt[1])
    #     tp_right = (c.top_midpt[0] + w / 2., c.top_midpt[1])
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    #     tp_left = lo_left
    #     tp_right = lo_right
    #     midx = (lo_left[0] + lo_right[0]) / 2.
    #     midy = tp_left[1] - taper_length
    #     lo_left = (midx - taper_to_width / 2., midy)
    #     lo_right = (midx + taper_to_width / 2., midy)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (2 in exclude):  # middle, bottom
    #     lo_left = (c.bottom_midpt[0] - w / 2., c.bottom_midpt[1])
    #     lo_right = (c.bottom_midpt[0] + w / 2., c.bottom_midpt[1])
    #     tp_left = (c.bottom_midpt[0] - w / 2., c.bottom_midpt[1] + (h - taper_length))
    #     tp_right = (c.bottom_midpt[0] + w / 2., c.bottom_midpt[1] + (h - taper_length))
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     lo_left = tp_left
    #     lo_right = tp_right
    #     midx = (lo_left[0] + lo_right[0]) / 2.
    #     midy = tp_left[1] + taper_length
    #     tp_left = (midx - taper_to_width / 2., midy)
    #     tp_right = (midx + taper_to_width / 2., midy)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (3 in exclude):  # Left, middle
    #     lo_left = (c.left_midpt[0], c.left_midpt[1] - w / 2.)
    #     lo_right = (c.left_midpt[0] + (h - taper_length), c.left_midpt[1] - w / 2.)
    #     tp_left = (c.left_midpt[0], c.left_midpt[1] + w / 2.)
    #     tp_right = (c.left_midpt[0] + (h - taper_length), c.left_midpt[1] + w / 2.)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     lo_left = lo_right
    #     tp_left = tp_right
    #     midx = tp_left[0] + taper_length
    #     midy = (lo_left[1] + tp_left[1]) / 2.
    #     tp_right = (midx, midy + taper_to_width / 2.)
    #     lo_right = (midx, midy - taper_to_width / 2.)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (4 in exclude):  # Right, middle
    #     lo_left = (c.right_midpt[0] - (h - taper_length), c.right_midpt[1] - w / 2.)
    #     lo_right = (c.right_midpt[0], c.right_midpt[1] - w / 2.)
    #     tp_left = (c.right_midpt[0] - (h - taper_length), c.right_midpt[1] + w / 2.)
    #     tp_right = (c.right_midpt[0], c.right_midpt[1] + w / 2.)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     lo_right = lo_left
    #     tp_right = tp_left
    #     midx = tp_left[0] - taper_length
    #     midy = (lo_left[1] + tp_left[1]) / 2.
    #     tp_left = (midx, midy + taper_to_width / 2.)
    #     lo_left = (midx, midy - taper_to_width / 2.)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (5 in exclude):  # Top, left
    #     lo_left = (c.top_left[0] - w / 2., c.top_left[1] - (h - taper_length))
    #     lo_right = (c.top_left[0] + w / 2., c.top_left[1] - (h - taper_length))
    #     tp_left = (c.top_left[0] - w / 2., c.top_left[1])
    #     tp_right = (c.top_left[0] + w / 2., c.top_left[1])
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     tp_left = lo_left
    #     tp_right = lo_right
    #     midx = (lo_left[0] + lo_right[0]) / 2.
    #     midy = tp_left[1] - taper_length
    #     lo_left = (midx - taper_to_width / 2., midy)
    #     lo_right = (midx + taper_to_width / 2., midy)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (6 in exclude):  # Top, right
    #     lo_left = (c.top_right[0] - w / 2., c.top_right[1] - (h - taper_length))
    #     lo_right = (c.top_right[0] + w / 2., c.top_right[1] - (h - taper_length))
    #     tp_left = (c.top_right[0] - w / 2., c.top_right[1])
    #     tp_right = (c.top_right[0] + w / 2., c.top_right[1])
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     tp_left = lo_left
    #     tp_right = lo_right
    #     midx = (lo_left[0] + lo_right[0]) / 2.
    #     midy = tp_left[1] - taper_length
    #     lo_left = (midx - taper_to_width / 2., midy)
    #     lo_right = (midx + taper_to_width / 2., midy)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (7 in exclude):  # Bottom, left
    #     lo_left = (c.bottom_left[0] - w / 2., c.bottom_left[1])
    #     lo_right = (c.bottom_left[0] + w / 2., c.bottom_left[1])
    #     tp_left = (c.bottom_left[0] - w / 2., c.bottom_left[1] + (h - taper_length))
    #     tp_right = (c.bottom_left[0] + w / 2., c.bottom_left[1] + (h - taper_length))
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     lo_left = tp_left
    #     lo_right = tp_right
    #     midx = (lo_left[0] + lo_right[0]) / 2.
    #     midy = tp_left[1] + taper_length
    #     tp_left = (midx - taper_to_width / 2., midy)
    #     tp_right = (midx + taper_to_width / 2., midy)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))
    #
    # if not (8 in exclude):  # Bottom, right
    #     lo_left = (c.bottom_right[0] - w / 2., c.bottom_right[1])
    #     lo_right = (c.bottom_right[0] + w / 2., c.bottom_right[1])
    #     tp_left = (c.bottom_right[0] - w / 2., c.bottom_right[1] + (h - taper_length))
    #     tp_right = (c.bottom_right[0] + w / 2., c.bottom_right[1] + (h - taper_length))
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, lo_right, tp_right, tp_left, lo_left]))
    #
    #     lo_left = tp_left
    #     lo_right = tp_right
    #     midx = (lo_left[0] + lo_right[0]) / 2.
    #     midy = tp_left[1] + taper_length
    #     tp_left = (midx - taper_to_width / 2., midy)
    #     tp_right = (midx + taper_to_width / 2., midy)
    #
    #     if solid:
    #         s.append(sdxf.Solid([lo_left, lo_right, tp_right, tp_left]))
    #     else:
    #         s.append(sdxf.PolyLine([lo_left, tp_left, tp_right, lo_right, lo_left]))


def draw_chip_alignment_marks(solid, d, c):
    """
    Draw the alignment marks on the chip.
    """
    CrossShapeAlignmentMarks(Structure(c, start=(125, 125), direction=90, defaults=d), width=2,
                             armlength=120, solid=solid, layer='0')
    CrossShapeAlignmentMarks(Structure(c, start=(c.size[0] - 125, c.size[1] - 125), direction=90, defaults=d),
                             width=2, armlength=120, solid=solid, layer='0')
    CrossShapeAlignmentMarks(Structure(c, start=(c.size[0] - 125, 125), direction=90, defaults=d),
                             width=2, armlength=120, solid=solid, layer='0')
    CrossShapeAlignmentMarks(Structure(c, start=(125, c.size[1] - 125), direction=90, defaults=d),
                             width=2, armlength=120, solid=solid, layer='0')


def draw_square_alignment_marks(structure, flip):
    # alignment boxes for e-beam
    c = structure

    if flip:
        direc = c.s2.last_direction

        c.s2.last_direction = 180+c.s2.last_direction

        align_box = 490  # center is the reference
        c.s2.last = (150, 18100-(8300 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 18100-(8300 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)

        c.s2.last = (150, 18100-(700 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 18100-(700 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        #  draw negative marker shadow
        align_box = 500  # center is the reference
        c.s2.last = (150, 18100-(8300 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 18100-(8300 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

        c.s2.last = (150, 18100-(700 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 18100-(700 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)


        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        c.s2.last_direction = direc

    else:

        align_box = 80  # center is the reference
        c.s2.last = (150, 28000 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 28000 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

        c.s2.last = (150, 700 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 700 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        #  draw negative marker shadow
        align_box = 500  # center is the reference
        c.s2.last = (150, 28000 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 28000 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

        c.s2.last = (150, 700 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1450, 700 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1000, 9800 - 50)
        CPWStraight(c.s2, 100, pinw=0, gapw=1500 / 2)
        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)




## Bridge Junction created by Ziqian Li
def RoundLCorner(structure, square, cornerl=10.0, gapw=2.0, whichedge=0):
    # starting point is at the outer corner of a L corner (L corner must have an outer corner)
    # whichedge = 0 or 1 : direction is along one edge or the other

    if whichedge == 0:
        angle = -90
    else:
        angle = 90
    structure.last_direction += angle
    structure.last_direction += angle

    structure.last_direction += angle
    structure.move(cornerl)
    structure.last_direction -= angle
    structure.last_direction += 2 * angle
    structure.move(cornerl * 2)
    structure.last_direction -= 2 * angle
    CPWBendNew(structure, angle=angle, pinw=0, gapw=cornerl, radius=cornerl, polyarc=1, segments=60,
               square=square)
    structure.move(-cornerl * 2)
    CPWStraight(structure, cornerl * 2, pinw=0, gapw=cornerl)


def RoundZCorner(structure, square, cornerl=10.0, gapw=2.0, cornertype=0):
    # starting point is at the outer corner of a Z corner
    # direction is pointing from the out corner to the inner corner
    # cornerl is the distance between the two corners of the z corner
    # cornertype = 0 or 1 : z or mirrored z

    if cornertype == 0:
        angle = -90
    else:
        angle = 90
    structure.move(0.75 * cornerl)
    structure.last_direction += angle
    CPWBendNew(structure, angle=angle, pinw=0, gapw=0.25 * cornerl, radius=0.25 * cornerl, polyarc=1, segments=60,
               square=square)
    ref = structure.last
    CPWStraight(structure, 0.5 * cornerl, 0, 0.25 * cornerl)
    structure.last = ref
    structure.last_direction -= angle
    structure.move(0.25 * cornerl)
    structure.last_direction += angle
    structure.move(0.25 * cornerl)
    structure.last_direction -= angle
    CPWStraight(structure, 0.5 * cornerl, 0, 0.25 * cornerl)
    structure.move(-0.25 * cornerl, structure.last_direction)
    structure.last_direction += angle
    structure.move(-0.25 * cornerl)

    CPWBend(structure, -angle, pinw=0, gapw=0.25 * cornerl, radius=0.25 * cornerl, polyarc=1, segments=60)
    # CPWBendNew(structure, angle=-angle, pinw=0, gapw=gapw, radius=0.5 * cornerl + gapw, polyarc=1, segments=60,
    #            square=square)

def CPWBendNew(structure, angle=90, pinw=10.0, gapw=20.0, radius=50.0, polyarc=1, segments=60, square=0):
    if square == 0:
        CPWBend(structure, angle, pinw=pinw, gapw=gapw, radius=radius, polyarc=polyarc, segments=segments)

    elif abs(angle) == 90:
        start_point = structure.last
        start_dir = structure.last_direction
        structure.move((gapw + pinw) / 2., structure.last_direction - (abs(angle) / angle) * 90)
        CPWStraight(structure, radius + gapw + pinw / 2., pinw=0, gapw=gapw / 2.)
        structure.last = start_point
        structure.move(radius + (gapw + pinw) / 2.)
        structure.move(pinw / 2., structure.last_direction - (abs(angle) / angle) * 90)
        structure.last_direction = structure.last_direction + (abs(angle) / angle) * 90
        CPWStraight(structure, pinw / 2. + radius, pinw=0, gapw=gapw / 2.)
        structure.last = start_point
        structure.last_direction = start_dir
        structure.move((gapw + pinw) / 2., structure.last_direction + (abs(angle) / angle) * 90)
        CPWStraight(structure, radius - pinw / 2., pinw=0, gapw=gapw / 2.)
        structure.move(-gapw / 2.)
        structure.move(gapw / 2., structure.last_direction - (abs(angle) / angle) * 90)
        structure.last_direction = structure.last_direction + (abs(angle) / angle) * 90
        CPWStraight(structure, radius - pinw / 2., pinw=0, gapw=gapw / 2.)
        structure.move((gapw + pinw) / 2., structure.last_direction - (abs(angle) / angle) * 90)

    elif abs(angle) == 180:
        start_point = structure.last
        start_dir = structure.last_direction
        structure.move((gapw + pinw) / 2., structure.last_direction - (abs(angle) / angle) * 90)
        CPWStraight(structure, radius + gapw + pinw / 2., pinw=0, gapw=gapw / 2.)
        structure.last = start_point
        structure.move(radius + (gapw + pinw) / 2.)
        structure.move(pinw / 2., structure.last_direction - (abs(angle) / angle) * 90)
        structure.last_direction = structure.last_direction + (abs(angle) / angle) * 90
        CPWStraight(structure, pinw + 2. * radius, pinw=0, gapw=gapw / 2.)
        structure.last = start_point
        structure.last_direction = start_dir
        structure.move((gapw + pinw) / 2, structure.last_direction + (abs(angle) / angle) * 90)
        CPWStraight(structure, radius - pinw / 2, pinw=0, gapw=gapw / 2)
        structure.move(-gapw / 2.)
        structure.move(gapw / 2., structure.last_direction - (abs(angle) / angle) * 90)
        structure.last_direction = structure.last_direction + (abs(angle) / angle) * 90
        CPWStraight(structure, 2. * radius - pinw, pinw=0, gapw=gapw / 2.)
        structure.move(-gapw / 2.)
        structure.move(gapw / 2., structure.last_direction - (abs(angle) / angle) * 90)
        structure.last_direction = structure.last_direction + (abs(angle) / angle) * 90
        CPWStraight(structure, radius - pinw / 2., pinw=0, gapw=gapw / 2.)
        structure.move(pinw + gapw, structure.last_direction - (abs(angle) / angle) * 90)
        structure.last_direction += 180
        CPWStraight(structure, radius + gapw + pinw / 2., pinw=0, gapw=gapw / 2.)
        structure.last_direction += 180
        structure.move(radius + gapw + pinw / 2.)
        structure.move((pinw + gapw) / 2., structure.last_direction + (abs(angle) / angle) * 90)

    else:
        print("this angle is not supported")


def curve_corner(structure, xpos, ypos, radius, type):
    c2 = structure
    ref_pos = c2.last
    ref_direction = c2.last_direction
    square = 0
    ########################################
    # Corner type: 0:left top, 1:right top, 2:right bottom, 3:left bottom
    #
    #
    #  1#0
    #################
    #  2#3
    #
    #
    # starting position at corner
    c2.last = (xpos, ypos)
    if type == 0:
        c2.last_direction = 0
        c2.move(radius / 2)
        c2.last_direction = 90
        CPWStraight(c2, radius, pinw=0, gapw=radius / 2)
        c2.last = (xpos, ypos)
        c2.last_direction = 0
        c2.move(radius / 2)
        c2.last_direction = 90
        c2.move(radius)
        c2.last_direction = 270
        CPWBendNew(c2, angle=90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=60,
                   square=square)
    if type == 1:
        c2.last_direction = 180
        c2.move(radius / 2)
        c2.last_direction = 90
        CPWStraight(c2, radius, pinw=0, gapw=radius / 2)
        c2.last = (xpos, ypos)
        c2.last_direction = 180
        c2.move(radius / 2)
        c2.last_direction = 90
        c2.move(radius)
        c2.last_direction = 270
        CPWBendNew(c2, angle=-90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=60,
                   square=square)
    if type == 2:
        c2.last_direction = 180
        c2.move(radius / 2)
        c2.last_direction = 270
        CPWStraight(c2, radius, pinw=0, gapw=radius / 2)
        c2.last = (xpos, ypos)
        c2.last_direction = 180
        c2.move(radius / 2)
        c2.last_direction = 270
        c2.move(radius)
        c2.last_direction = 90
        CPWBendNew(c2, angle=90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=60,
                   square=square)
    if type == 3:
        c2.last_direction = 0
        c2.move(radius / 2)
        c2.last_direction = 270
        CPWStraight(c2, radius, pinw=0, gapw=radius / 2)
        c2.last = (xpos, ypos)
        c2.last_direction = 0
        c2.move(radius / 2)
        c2.last_direction = 270
        c2.move(radius)
        c2.last_direction = 90
        CPWBendNew(c2, angle=-90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=60,
                   square=square)
    ########################################
    c2.last_direction = ref_direction
    c2.last = ref_pos


def perforate(chip, grid_x, grid_y):
    nx, ny = map(int, [chip.size[0] / grid_x, chip.size[1] / grid_y])
    occupied = [[False] * ny for i in range(nx)]
    for i in range(nx):
        occupied[i][0] = True
        occupied[i][-1] = True
    for i in range(ny):
        occupied[0][i] = True
        occupied[-1][i] = True

    for e in chip.entities:
        o_x_list = []
        o_y_list = []
        for p in e.points:
            o_x, o_y = map(int, (p[0] / grid_x, p[1] / grid_y))
            if 0 <= o_x < nx and 0 <= o_y < ny:
                o_x_list.append(o_x)
                o_y_list.append(o_y)
        if o_x_list:
            for x in range(min(o_x_list), max(o_x_list) + 1):
                for y in range(min(o_y_list), max(o_y_list) + 1):
                    occupied[x][y] = True

    second_pass = deepcopy(occupied)
    for i in range(nx):
        for j in range(ny):
            if occupied[i][j]:
                for ip, jp in [(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)]:
                    try:
                        second_pass[ip][jp] = True
                    except IndexError:
                        pass

    for i in range(nx):
        for j in range(ny):
            if not second_pass[i][j]:
                size = 2.5
                pos = i * grid_x + grid_x / 2., j * grid_y + grid_y / 2.
                p0 = vadd(pos, (-size, -size))
                p1 = vadd(pos, (size, size))
                abs_rect(chip, p0, p1)


def Draw_test_large(structure, xpos, ypos):
    c2 = structure
    ##########################################
    ref_pos1 = c2.s5.last
    ref_direction = c2.s5.last_direction
    pad_x = 700
    pad_y = 350
    gap = 50
    length = 200
    cpw_pinw = 10.0
    cpw_gapw = calculate_gap_width(eps_eff, 50, cpw_pinw)
    # Default settings for junction box (together with top flux line)
    xpos1 = xpos
    ypos1 = ypos
    shape3 = {}
    shape3['curve'] = 5
    shape3['x'] = 60
    shape3['y1'] = 30
    shape3['y2'] = 30
    shape3['flux_dis'] = 2
    shape3['flux_shift'] = 0  # x position shift relative to the center
    shape3['extension'] = 100  # extending lengtth of CPW
    shift_x = [-300,-200, -100, 100, 200, 300, 0]
    for sx in shift_x:
        junctionbox_test(c2, xpos1+sx, ypos1, shape3)

    c2.s5.last_direction = 90

    for sx in shift_x:
        c2.s5.last = (xpos1+sx, ypos1 + shape3['x'] / 2 + shape3['extension'])
        CPWStraight(c2.s5, length, pinw=cpw_pinw, gapw=cpw_gapw)
        CPWStraight(c2.s5, gap, pinw=0, gapw=(cpw_pinw)/2)
    c2.s5.move(-gap)
    CPWStraight(c2.s5, gap, pinw=0, gapw=(gap * 2 + pad_x) / 2)
    CPWStraight(c2.s5, pad_y, pinw=pad_x, gapw=gap)
    CPWStraight(c2.s5, gap, pinw=0, gapw=gap + pad_x / 2)

    c2.s5.last_direction = 270
    for sx in shift_x:
        c2.s5.last = (xpos1+sx, ypos1 - shape3['x'] / 2 - shape3['extension'])
        CPWStraight(c2.s5, length, pinw=cpw_pinw, gapw=cpw_gapw)
        CPWStraight(c2.s5, gap, pinw=0, gapw=(cpw_pinw) / 2)
    c2.s5.move(-gap)
    CPWStraight(c2.s5, gap, pinw=0, gapw=(gap * 2 + pad_x ) / 2)
    CPWStraight(c2.s5, pad_y, pinw=pad_x, gapw=gap)
    CPWStraight(c2.s5, gap, pinw=0, gapw=gap + pad_x / 2)

    ########################################
    c2.s5.last_direction = ref_direction
    c2.s5.last = ref_pos1


def draw_junction_frame(structure, xpos, ypos, shape, flag):
    c = structure
    pin_ref = c.s2.pin_layer.last
    pin_direction = c.s2.pin_layer.last_direction
    gap_ref = c.s2.gap_layer.last
    gap_direction = c.s2.gap_layer.last_direction

    gap_new = 0.7
    gap_thin = 0.1

    if flag == 0:

        c.s2.pin_layer.last_direction = 90
        c.s2.pin_layer.last = (xpos, ypos)
        CPWStraight(c.s2.pin_layer, shape['overlap'], pinw=0, gapw=shape['overlap_width']/2)
        CPWLinearTaper(c.s2.pin_layer, shape['shrink_length'], start_pinw=0, stop_pinw=0, start_gapw=shape['overlap_width']/2,
                       stop_gapw=shape['mjj_bar_width']/2)
        CPWStraight(c.s2.pin_layer, shape['mjj_x_bar_length'], pinw=0, gapw=shape['mjj_bar_width'] / 2)

        pin_store1 = c.s2.pin_layer.last

        CPWStraight(c.s2.pin_layer, shape['mjj_x_length'], pinw=0, gapw=shape['mjj_x_width'] / 2)

        c.s2.pin_layer.last = pin_store1
        c.s2.pin_layer.move(3+shape['mjj_y_width'] / 2)

        c.s2.pin_layer.last_direction = 0
        c.s2.pin_layer.move(shape['offset']+shape['mjj_bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, -shape['offset']-shape['mjj_bar_width'] / 2-shape['mjj_x_width']/2, pinw=0, gapw=shape['mjj_y_width'] / 2)
        c.s2.pin_layer.move(shape['mjj_x_width'])
        CPWStraight(c.s2.pin_layer, shape['mjj_y_length']+shape['offset']+shape['mjj_bar_width'] / 2-shape['mjj_x_width']/2, pinw=0, gapw=shape['mjj_y_width'] / 2)

        c.s2.pin_layer.last_direction = 90
        c.s2.pin_layer.last = pin_store1
        c.s2.pin_layer.last = (xpos+shape['offset'], ypos+shape['overlap']*2+30.0)
        c.s2.pin_layer.last_direction = 270
        CPWStraight(c.s2.pin_layer, shape['overlap'], pinw=0, gapw=shape['overlap_width'] / 2)
        CPWLinearTaper(c.s2.pin_layer, shape['shrink_length'], start_pinw=0, stop_pinw=0,
                       start_gapw=shape['overlap_width'] / 2,
                       stop_gapw=shape['mjj_bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, shape['mjj_y_bar_length'], pinw=0, gapw=shape['mjj_bar_width'] / 2)

        c.s2.gap_layer.last_direction = 270
        c.s2.gap_layer.last = c.s2.pin_layer.last
        c.s2.gap_layer.last_direction = 180
        c.s2.gap_layer.move(gap_new/2)
        c.s2.gap_layer.last_direction = 270
        CPWStraight(c.s2.gap_layer, gap_new, pinw=0, gapw=shape['mjj_bar_width'] / 2+gap_new/2)

        c.s2.gap_layer.last_direction = 180
        c.s2.gap_layer.last = (c.s2.pin_layer.last[0]-shape['mjj_bar_width'] / 2, c.s2.pin_layer.last[1]+shape['mjj_y_bar_length'] / 2)
        CPWStraight(c.s2.gap_layer, gap_new, pinw=0, gapw=shape['mjj_y_bar_length'] / 2)

        c.s2.gap_layer.last_direction = 90
        c.s2.gap_layer.last = (
        c.s2.pin_layer.last[0], c.s2.pin_layer.last[1] + shape['mjj_y_bar_length'])
        CPWLinearTaper(c.s2.gap_layer, shape['shrink_length'], start_pinw=shape['mjj_bar_width'], stop_pinw=shape['overlap_width'],
                       start_gapw=gap_new,
                       stop_gapw=gap_new)
        CPWStraight(c.s2.gap_layer, shape['overlap'], pinw=shape['overlap_width'], gapw=gap_new)
        CPWStraight(c.s2.gap_layer, gap_new, pinw=0, gapw=gap_new+shape['overlap_width']/2)

        c.s2.gap_layer.last_direction = 180
        c.s2.gap_layer.last = (
        pin_store1[0] - shape['mjj_bar_width'] / 2, pin_store1[1] - shape['mjj_y_bar_length'] / 2)
        CPWStraight(c.s2.gap_layer, gap_new, pinw=0, gapw=shape['mjj_y_bar_length'] / 2)
        c.s2.gap_layer.last_direction = 0
        c.s2.gap_layer.last = (
            pin_store1[0] + shape['mjj_bar_width'] / 2, pin_store1[1] - shape['mjj_y_bar_length'] / 2)
        CPWStraight(c.s2.gap_layer, gap_new, pinw=0, gapw=shape['mjj_y_bar_length'] / 2)

        c.s2.gap_layer.last_direction = 90
        c.s2.gap_layer.last = (
            pin_store1[0], pin_store1[1])
        CPWStraight(c.s2.gap_layer, gap_new, pinw=shape['mjj_y_width'], gapw=gap_new -shape['mjj_y_width']/2+ shape['mjj_bar_width'] / 2)

        c.s2.gap_layer.last_direction = 270
        c.s2.gap_layer.last = (
        pin_store1[0], pin_store1[1] - shape['mjj_y_bar_length'])
        CPWLinearTaper(c.s2.gap_layer, shape['shrink_length'], start_pinw=shape['mjj_bar_width'],
                       stop_pinw=shape['overlap_width'],
                       start_gapw=gap_new,
                       stop_gapw=gap_new)
        CPWStraight(c.s2.gap_layer, shape['overlap'], pinw=shape['overlap_width'], gapw=gap_new)
        CPWStraight(c.s2.gap_layer, gap_new, pinw=0, gapw=gap_new + shape['overlap_width'] / 2)

    else:
        c.s2.pin_layer.last_direction = 90
        c.s2.pin_layer.last = (xpos, ypos)
        CPWStraight(c.s2.pin_layer, shape['overlap'], pinw=0, gapw=shape['overlap_width'] / 2)
        CPWLinearTaper(c.s2.pin_layer, shape['shrink_length'], start_pinw=0, stop_pinw=0,
                       start_gapw=shape['overlap_width'] / 2,
                       stop_gapw=shape['mjj_bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, shape['mjj_x_bar_length'], pinw=0, gapw=shape['mjj_bar_width'] / 2)

        pin_store = c.s2.pin_layer.last

        CPWStraight(c.s2.pin_layer, shape['mjj_x_length'], pinw=0, gapw=shape['mjj_x_width'] / 2)

        c.s2.pin_layer.last = pin_store
        c.s2.pin_layer.move(3 + shape['mjj_y_width'] / 2)

        c.s2.pin_layer.last_direction = 0
        c.s2.pin_layer.move(shape['offset'] + shape['mjj_bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, -shape['offset'] - shape['mjj_bar_width'] / 2 - shape['mjj_x_width'] / 2, pinw=0,
                    gapw=shape['mjj_y_width'] / 2)
        c.s2.pin_layer.move(shape['mjj_x_width'])
        CPWStraight(c.s2.pin_layer,
                    shape['mjj_y_length'] + shape['offset'] + shape['mjj_bar_width'] / 2 - shape['mjj_x_width'] / 2,
                    pinw=0, gapw=shape['mjj_y_width'] / 2)

        c.s2.pin_layer.last_direction = 90
        c.s2.pin_layer.last = pin_store
        c.s2.pin_layer.last = (xpos + shape['offset'], ypos + shape['overlap'] * 2 + 30.0)
        c.s2.pin_layer.last_direction = 270
        CPWStraight(c.s2.pin_layer, shape['overlap'], pinw=0, gapw=shape['overlap_width'] / 2)
        CPWLinearTaper(c.s2.pin_layer, shape['shrink_length'], start_pinw=0, stop_pinw=0,
                       start_gapw=shape['overlap_width'] / 2,
                       stop_gapw=shape['mjj_bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, shape['mjj_y_bar_length'], pinw=0, gapw=shape['mjj_bar_width'] / 2)


    ################################################
    c.s2.pin_layer.last = pin_ref
    c.s2.gap_layer.last = gap_ref
    c.s2.pin_layer.last_direction = pin_direction
    c.s2.gap_layer.last_direction = gap_direction


def draw_test_frame(structure, xpos, ypos, junc_para, loop):
    c = structure
    pin_ref = c.s2.pin_layer.last
    pin_direction = c.s2.pin_layer.last_direction
    gap_ref = c.s2.gap_layer.last
    gap_direction = c.s2.gap_layer.last_direction

    c.s2.pin_layer.last = (xpos, ypos)
    arm = loop['arm']
    boxx = 10
    boxy = 7
    c.s2.pin_layer.last_direction = 90
    CPWStraight(c.s2.pin_layer, 30, pinw=0, gapw=arm / 2)
    CPWStraight(c.s2.pin_layer, boxy, pinw=0, gapw=boxx / 2)

    c.s2.pin_layer.last = (xpos, ypos)
    c.s2.pin_layer.last_direction = 270
    CPWStraight(c.s2.pin_layer, 30, pinw=0, gapw=arm / 2)
    CPWStraight(c.s2.pin_layer, boxy, pinw=0, gapw=boxx / 2)

    ################################################
    c.s2.pin_layer.last = pin_ref
    c.s2.gap_layer.last = gap_ref
    c.s2.pin_layer.last_direction = pin_direction
    c.s2.gap_layer.last_direction = gap_direction


def draw_junction_small(structure, xpos, ypos, junc_correction, loop, junction_number, flag):
    c = structure
    pin_ref = c.s2.pin_layer.last
    pin_direction = c.s2.pin_layer.last_direction
    gap_ref = c.s2.gap_layer.last
    gap_direction = c.s2.gap_layer.last_direction
    linear = junc_correction['junction_linear']
    gap_protect = junc_correction['gap_protect']
    junc_arm = junc_correction['junc_arm']
    looparm = loop['arm']
    if junction_number == 3:
        junction_width = junc_correction['junction3_length']
        junction_gap = junc_correction['junction3_gap']
    if junction_number == 4:
        junction_width = junc_correction['junction4_length']
        junction_gap = junc_correction['junction4_gap']
    if junction_number == 5:
        junction_width = junc_correction['junction5_length']
        junction_gap = junc_correction['junction5_gap']
    if junction_number == 6:
        junction_width = junc_correction['junction6_length']
        junction_gap = junc_correction['junction6_gap']
    if junction_number == 7:
        junction_width = junc_correction['junction7_length']
        junction_gap = junc_correction['junction7_gap']
    if junction_number == 8:
        junction_width = junc_correction['junction8_length']
        junction_gap = junc_correction['junction8_gap']

    block_length = junction_gap + junc_arm + linear
    if (flag == 0) or (flag == -1):
        # Clear the junction block
        c.s2.pin_layer.last = (xpos, ypos)
        c.s2.pin_layer.last_direction = 90
        CPWStraight(c.s2.pin_layer, block_length, pinw=0, gapw=looparm / 2)
    if flag == -1:
        # Draw junction
        c.s2.pin_layer.last = (xpos, ypos + junction_gap)
        c.s2.pin_layer.last_direction = 90
        CPWStraight(c.s2.pin_layer, junc_arm, pinw=0, gapw=junction_width / 2)
        CPWLinearTaper(c.s2.pin_layer, linear, 0, 0, junction_width / 2, looparm / 2)
        # Draw gap layer
        gap_wrap = 3
        c.s2.gap_layer.last = (xpos, ypos - gap_wrap)
        c.s2.gap_layer.last_direction = 90
        CPWStraight(c.s2.gap_layer, gap_wrap, pinw=looparm, gapw=gap_protect)
        CPWStraight(c.s2.gap_layer, junction_gap, pinw=0, gapw=looparm / 2 + gap_protect)
        CPWStraight(c.s2.gap_layer, junc_arm, pinw=junction_width, gapw=gap_protect)
        CPWLinearTaper(c.s2.gap_layer, linear, junction_width, looparm, gap_protect, gap_protect)
        CPWStraight(c.s2.gap_layer, gap_wrap, pinw=looparm, gapw=gap_protect)

    ################################################
    c.s2.pin_layer.last = pin_ref
    c.s2.gap_layer.last = gap_ref
    c.s2.pin_layer.last_direction = pin_direction
    c.s2.gap_layer.last_direction = gap_direction


def multimode_qubit(structure, xpos, ypos, shape, flag):
    c = structure
    pin_ref = c.s2.last
    pin_direction = c.s2.last_direction
    ################################################
    c.s2.last = (xpos, ypos)
    c.s2.last_direction = 90+180*flag

    CPWStraight(c.s2, shape['y_left_thick'], pinw=0, gapw=shape['x_left_thick'] / 2)
    CPWStraight(c.s2, shape['y_left_thin'], pinw=0, gapw=shape['x_thin'] / 2)

    #  draw overlapping T shape
    x_last = c.s2.last[0]
    y_last = c.s2.last[1]
    c.s2.last_direction = 180 + c.s2.last_direction
    T_block = 5
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block / 2)
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block * 3 / 2)
    c.s2.last_direction = 180 + c.s2.last_direction
    c.s2.last = (x_last, y_last)

    c.s2.move(shape['jjgap'])

    #  draw overlapping T shape
    x_last = c.s2.last[0]
    y_last = c.s2.last[1]
    T_block = 5
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block / 2)
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block * 3 / 2)
    c.s2.last = (x_last, y_last)

    CPWStraight(c.s2, shape['y_right_thin'], pinw=0, gapw=shape['x_thin'] / 2)
    CPWStraight(c.s2, shape['y_right_thick'], pinw=0, gapw=shape['x_right_thick'] / 2)
    x_last = c.s2.last[0]
    y_last = c.s2.last[1]


    # curve all corners
    curve_rad_large = 40
    curve_rad_small = 10
    if flag==0:
        curve_corner(c.s2, xpos + shape['x_left_thick'] / 2, ypos, curve_rad_large, 1)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos, curve_rad_large, 0)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos + shape['y_left_thick'], curve_rad_large, 3)
        curve_corner(c.s2, xpos + shape['x_left_thick'] / 2, ypos + shape['y_left_thick'], curve_rad_large, 2)

        curve_corner(c.s2, xpos - shape['x_thin'] / 2, ypos + shape['y_left_thick'], curve_rad_small, 1)
        curve_corner(c.s2, xpos + shape['x_thin'] / 2, ypos + shape['y_left_thick'], curve_rad_small, 0)
        #
        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last, curve_rad_large, 2)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last, curve_rad_large, 3)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last - shape['y_right_thick'], curve_rad_large, 0)
        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last - shape['y_right_thick'], curve_rad_large, 1)
        #
        curve_corner(c.s2, x_last - shape['x_thin'] / 2, y_last - shape['y_right_thick'], curve_rad_small, 2)
        curve_corner(c.s2, x_last + shape['x_thin'] / 2, y_last - shape['y_right_thick'], curve_rad_small, 3)
    else:
        curve_corner(c.s2, xpos+shape['x_left_thick'] / 2, ypos, curve_rad_large, 2)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos, curve_rad_large, 3)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos-shape['y_left_thick'], curve_rad_large, 0)
        curve_corner(c.s2, xpos + shape['x_left_thick'] / 2, ypos-shape['y_left_thick'], curve_rad_large, 1)

        curve_corner(c.s2, xpos - shape['x_thin'] / 2, ypos - shape['y_left_thick'], curve_rad_small, 2)
        curve_corner(c.s2, xpos + shape['x_thin'] / 2, ypos - shape['y_left_thick'], curve_rad_small, 3)


        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last, curve_rad_large, 1)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last, curve_rad_large, 0)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last + shape['y_right_thick'], curve_rad_large, 3)
        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last + shape['y_right_thick'], curve_rad_large, 2)

        curve_corner(c.s2, x_last - shape['x_thin'] / 2, y_last + shape['y_right_thick'], curve_rad_small, 1)
        curve_corner(c.s2, x_last + shape['x_thin'] / 2, y_last + shape['y_right_thick'], curve_rad_small, 0)
    ###############################################
    # draw surface top marker
    c.s2.move(8000+1900)
    x_mark = c.s2.last[0]
    y_mark = c.s2.last[1]
    if flag == 0:

        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)
        CPWStraight(c.s2, 500, pinw=200, gapw=100)
        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)

        c.s2.last = (x_mark + 600, y_mark)
        CPWStraight(c.s2, 100, pinw=0, gapw=500 / 2)
        CPWStraight(c.s2, 600, pinw=0, gapw=100 / 2)

        c.s2.last = (x_mark - 600, y_mark)
        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)
        CPWStraight(c.s2, 150, pinw=200, gapw=100)
        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)
        c.s2.last = (x_mark - 450, y_mark + 350)
        CPWStraight(c.s2, 350, pinw=0, gapw=100 / 2)

    else:

        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)
        CPWStraight(c.s2, 500, pinw=200, gapw=100)
        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)

        c.s2.last = (x_mark-600, y_mark)
        CPWStraight(c.s2, 100, pinw=0, gapw=500 / 2)
        CPWStraight(c.s2, 600, pinw=0, gapw=100/2)

        c.s2.last = (x_mark + 600, y_mark)
        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)
        CPWStraight(c.s2, 150, pinw=200, gapw=100)
        CPWStraight(c.s2, 100, pinw=0, gapw=400 / 2)
        c.s2.last = (x_mark + 450, y_mark - 350)
        CPWStraight(c.s2, 350, pinw=0, gapw=100 / 2)


    ################################################
    c.s2.last = pin_ref
    c.s2.last_direction = pin_direction

def draw_test_optical(structure, xpos, ypos, shape, flag):
    c = structure
    pin_ref = c.s2.last
    pin_direction = c.s2.last_direction
    ################################################
    c.s2.last = (xpos, ypos)
    c.s2.last_direction = 90+180*flag

    CPWStraight(c.s2, shape['y_left_thick'], pinw=0, gapw=shape['x_left_thick'] / 2)
    CPWStraight(c.s2, shape['y_left_thin'], pinw=0, gapw=shape['x_thin'] / 2)

    #  draw overlapping T shape
    x_last = c.s2.last[0]
    y_last = c.s2.last[1]
    c.s2.last_direction = 180 + c.s2.last_direction
    T_block = 5
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block / 2)
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block * 3 / 2)
    c.s2.last_direction = 180 + c.s2.last_direction
    c.s2.last = (x_last, y_last)

    c.s2.move(shape['jjgap'])

    #  draw overlapping T shape
    x_last = c.s2.last[0]
    y_last = c.s2.last[1]
    T_block = 5
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block / 2)
    CPWStraight(c.s2, T_block, pinw=0, gapw=T_block * 3 / 2)
    c.s2.last = (x_last, y_last)

    CPWStraight(c.s2, shape['y_right_thin'], pinw=0, gapw=shape['x_thin'] / 2)
    CPWStraight(c.s2, shape['y_right_thick'], pinw=0, gapw=shape['x_right_thick'] / 2)
    x_last = c.s2.last[0]
    y_last = c.s2.last[1]


    # curve all corners
    curve_rad_large = 40
    curve_rad_small = 10
    if flag==0:
        curve_corner(c.s2, xpos + shape['x_left_thick'] / 2, ypos, curve_rad_large, 1)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos, curve_rad_large, 0)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos + shape['y_left_thick'], curve_rad_large, 3)
        curve_corner(c.s2, xpos + shape['x_left_thick'] / 2, ypos + shape['y_left_thick'], curve_rad_large, 2)

        curve_corner(c.s2, xpos - shape['x_thin'] / 2, ypos + shape['y_left_thick'], curve_rad_small, 1)
        curve_corner(c.s2, xpos + shape['x_thin'] / 2, ypos + shape['y_left_thick'], curve_rad_small, 0)
        #
        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last, curve_rad_large, 2)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last, curve_rad_large, 3)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last - shape['y_right_thick'], curve_rad_large, 0)
        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last - shape['y_right_thick'], curve_rad_large, 1)
        #
        curve_corner(c.s2, x_last - shape['x_thin'] / 2, y_last - shape['y_right_thick'], curve_rad_small, 2)
        curve_corner(c.s2, x_last + shape['x_thin'] / 2, y_last - shape['y_right_thick'], curve_rad_small, 3)
    else:
        curve_corner(c.s2, xpos+shape['x_left_thick'] / 2, ypos, curve_rad_large, 2)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos, curve_rad_large, 3)
        curve_corner(c.s2, xpos - shape['x_left_thick'] / 2, ypos-shape['y_left_thick'], curve_rad_large, 0)
        curve_corner(c.s2, xpos + shape['x_left_thick'] / 2, ypos-shape['y_left_thick'], curve_rad_large, 1)

        curve_corner(c.s2, xpos - shape['x_thin'] / 2, ypos - shape['y_left_thick'], curve_rad_small, 2)
        curve_corner(c.s2, xpos + shape['x_thin'] / 2, ypos - shape['y_left_thick'], curve_rad_small, 3)


        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last, curve_rad_large, 1)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last, curve_rad_large, 0)
        curve_corner(c.s2, x_last - shape['x_left_thick'] / 2, y_last + shape['y_right_thick'], curve_rad_large, 3)
        curve_corner(c.s2, x_last + shape['x_left_thick'] / 2, y_last + shape['y_right_thick'], curve_rad_large, 2)

        curve_corner(c.s2, x_last - shape['x_thin'] / 2, y_last + shape['y_right_thick'], curve_rad_small, 1)
        curve_corner(c.s2, x_last + shape['x_thin'] / 2, y_last + shape['y_right_thick'], curve_rad_small, 0)
    ###############################################
    # draw surface top marker
    c.s2.move(7000)
    x_mark = c.s2.last[0]
    y_mark = c.s2.last[1]


    ################################################
    c.s2.last = pin_ref
    c.s2.last_direction = pin_direction


def draw_bridge_junction(structure, d, junc_correction, flip):
    c = structure
    ## alignment boxes for e-beam
    draw_launchers(c, d, exclude=[2, 3, 8])
    draw_launchers(c, d, exclude=[2, 3, 8])
    draw_square_alignment_marks(c, flip)
    ## Junction parameters

    xpos = 1000.
    ypos = 17700.
    if flip:
        ypos = 18500-17700.

    shape = {}
    shape['x_thin'] = 30
    shape['y_right_thin'] = 350
    shape['x_right_thick'] = 150
    shape['y_right_thick'] = 3100

    shape['y_left_thin'] = 370
    shape['x_left_thick'] = 150
    shape['y_left_thick'] = 3600+100

    shape['jjgap'] = 30

    shape1 = {}
    shape1['x_thin'] = 30
    shape1['y_right_thin'] = 100
    shape1['x_right_thick'] = 700
    shape1['y_right_thick'] = 700

    shape1['y_left_thin'] = 100
    shape1['x_left_thick'] = 700
    shape1['y_left_thick'] = 700

    shape1['jjgap'] = 30
    if draw_layer1 == True:
        # Draw Multimode star pad
        flag = 1
        if flip:
            flag = 0
        xpos = 1000.-200.
        ypos = 17700.-200.+7780
        if flip:
            ypos = 18100 - ypos
        multimode_qubit(c, xpos, ypos-400, shape, flag)
        # draw a few test frames

        draw_test_optical(c, xpos-500, ypos-20000+300-0, shape1, flag)
        draw_test_optical(c, xpos + 500, ypos - 20000+300-0, shape1, flag)

        draw_test_optical(c, xpos - 500, ypos - 22500 + 300-0, shape1, flag)
        draw_test_optical(c, xpos + 500, ypos - 22500 + 300-0, shape1, flag)

    pin_length = 0
    if draw_layer2 == True:

        # in total four junctions
        junc_para = {}
        junc_para['shrink_length'] = 7.0
        junc_para['overlap'] = 20.0
        junc_para['overlap_width'] = 25.0
        junc_para['offset'] = -2.5
        junc_para['mjj_x_width'] = junc_correction['mjjx']
        junc_para['mjj_y_width'] = junc_correction['mjjy']
        junc_para['mjj_x_length'] = 4.0+1.3
        junc_para['mjj_y_length'] = 3.0+1.3
        junc_para['mjj_x_bar_length'] = 6.5
        junc_para['mjj_y_bar_length'] = 6.5
        junc_para['mjj_bar_width'] = 1.0
        junc_para['gap'] = 0.2

        # flag: direction of the junction
        xpos = 1000.-200.
        ypos = 14200.-720+7780-500
        flag1 = 0
        if flip:
            flag1 = 1
            xpos = 1000. - 200.
            ypos = 18030-(14200. - 720)
        draw_junction_frame(c, xpos, ypos, junc_para, flag1)

        # draw test structures
        # junc_para['mjj_x_width'] = junc_correction['mjjx1']
        # junc_para['mjj_y_width'] = junc_correction['mjjy1']
        draw_junction_frame(c, xpos-500, ypos-20000+300+3170+500, junc_para, flag1)
        # junc_para['mjj_x_width'] = junc_correction['mjjx2']
        # junc_para['mjj_y_width'] = junc_correction['mjjy2']
        draw_junction_frame(c, xpos+500, ypos-20000+300+3170+500, junc_para, flag1)
        junc_para['mjj_x_width'] = junc_correction['mjjx3']
        junc_para['mjj_y_width'] = junc_correction['mjjy3']
        draw_junction_frame(c, xpos-500, ypos- 22500 + 300+3170+500, junc_para, flag1)
        junc_para['mjj_x_width'] = junc_correction['mjjx3']
        junc_para['mjj_y_width'] = junc_correction['mjjx3']
        draw_junction_frame(c, xpos+500, ypos- 22500 + 300+3170+500, junc_para, flag1)




def draw_large_junction_test(structure, d, junc_correction):
    c = structure
    ## alignment boxes for e-beam
    if (draw_layer1 == True) and (draw_layer2 == False):
        draw_launchers(c, d, exclude=[2, 4, 8])
        draw_launchers(c, d, exclude=[2, 4, 8])
        draw_square_alignment_marks(c)
    else:
        draw_launchers(c, d, exclude=[2, 4, 8])
        draw_launchers(c, d, exclude=[2, 4, 8])
        draw_square_alignment_marks(c)
    ## Junction parameters
    junc_length = junc_correction['s_length']
    junc_gap = junc_correction['s_gap']
    large_length = junc_correction['l_length']
    large_gap = junc_correction['l_gap']
    length_step = junc_correction['length_step']
    length_step_large = junc_correction['length_step_large']

    xpos = (c.s1.last[0] + c.s5.last[0]) / 2.
    ypos = (c.s3.last[1] + c.s8.last[1]) / 2.
    # Draw VSLQ test structure
    ylist = [1500, 3400, 5300]
    xlist = [770, 1910, 3050, 4190,  5330]
    if draw_layer1 == True:

        for i in range(len(xlist)):
            for j in range(len(ylist)):
                Draw_test_large(c, xlist[i], ylist[j])
    if draw_layer2 == True:
        # in total four junctions
        junc_para = {}

        junc_para['junction_linear'] = 3
        junc_para['junction_linear_large'] = 2

        junc_para['gap_protect'] = 0.2
        junc_para['junc_arm'] = 0.5

        # define loop parameter
        loop = {}
        loop['x'] = 15
        loop['y'] = 15
        # position of the SQUID junctions (ratio)
        loop['position'] = 0.5
        loop['arm'] = 2

        # draw test part
        single = [-300, -200, -100, 100, 200, 300, 0]
        for i in range(len(xlist)):
            junc_para['junction1_length'] = large_length + i * length_step_large
            junc_para['junction3_length'] = junc_length + i * length_step

            junc_para['junction1_gap'] = large_gap
            junc_para['junction3_gap'] = junc_gap
            if i < len(xlist) - 1:
                ff = -1
            else:
                ff = 1
            for j in range(len(ylist)):
                for k in single:
                    draw_test_frame(structure, xlist[i]+k, ylist[j], junc_para, loop)
                    draw_junction_large(structure, xlist[i]+k, ylist[j], junc_para, loop, 1, ff)


def output_wafer():
    # Meant for writing junctions on a wafer with metallic base layer

    m = WaferMask(MaskName, diameter=50800., flat_angle=270., flat_distance=24100., wafer_padding=2500,
                  chip_size=(1600, 25950),
                  dicing_border=400, etchtype=False, wafer_edge=False,
                  dashed_dicing_border=80, ndashes=3, dice_corner=True, square_arr=False)

    # 4in wafer -> diameter=101600, flat_distance=48200
    # 3in wafer -> diameter=76200., flat_distance=37100.
    # 2in wafer -> diameter=50800., flat_distance=24100.

    points = [(-24000., -4000.), (-24000., 4000.), (24000., -4000.), (24000., 4000.)]
    points_medium = [(-24000. + 1200, -4000.), (-24000. + 1200, 4000.), (24000. - 1200, -4000.), (24000. - 1200, 4000.)]
    points_small = [(-24000. + 2000, -4000.), (-24000. + 2000, 4000.), (24000. - 2000, -4000.), (24000. - 2000, 4000.)]

    # Create the alignment crosses on the wafer. NOTE: currently only works for solid = True setting.
    # Something goes wrong with the layers cause that's built in otherwise. Hacky solution:
    if solid:
        if draw_layer1 == 2:
            AlignmentCross(m, linewidth=50, size=(1000, 1000), solid=solid, points=points, layer=None, name='cross')
            AlignmentCross(m, linewidth=25, size=(200, 200), solid=solid, points=points_medium, layer=None,
                           name='cross_medium')
            AlignmentCross(m, linewidth=5, size=(20, 20), solid=solid, points=points_small, layer=None,
                           name='cross_small')


    else:
        if draw_layer1 == 2:
            AlignmentCross(m, linewidth=50, size=(1000, 1000), solid=False, points=points, layer='0', name='cross')
            AlignmentCross(m, linewidth=25, size=(200, 200), solid=False, points=points_medium, layer='0',
                           name='cross_medium')
            AlignmentCross(m, linewidth=5, size=(20, 20), solid=False, points=points_small, layer='0',
                           name='cross_small')

    d = set_mask_init()

    # --------------------------#

    junc_correction = {'mjjx': 0.14, 'mjjy': 0.14,'mjjx3':0.14, 'mjjy3':0.14}
    CHIPNAME = 'MM2A'
    flip = False
    c0 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(100, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c0, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.14, 'mjjy': 0.14,'mjjx3':0.14, 'mjjy3':0.14}
    flip = False
    CHIPNAME = 'MM2B'
    c1 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(100, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c1, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.155, 'mjjy': 0.155,'mjjx3':0.155, 'mjjy3':0.155}
    flip = False
    CHIPNAME = 'MM2C'
    c2 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    # draw_bridge_junction(c2, d, junc_correction)
    draw_bridge_junction(c2, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.155, 'mjjy': 0.155,'mjjx3':0.155, 'mjjy3':0.155}
    flip = False
    CHIPNAME = 'MM2D'
    c3 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c3, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.17, 'mjjy': 0.17,'mjjx3':0.17, 'mjjy3':0.17}
    flip = False
    CHIPNAME = 'MM2E'
    c4 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c4, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.17, 'mjjy': 0.17,'mjjx3':0.17, 'mjjy3':0.17}
    flip = False
    CHIPNAME = 'MM2F'
    c5 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c5, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.185, 'mjjy': 0.185,'mjjx3':0.185, 'mjjy3':0.185}
    flip = False
    CHIPNAME = 'MM2G'
    c6 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c6, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.185, 'mjjy': 0.185,'mjjx3':0.185, 'mjjy3':0.185}
    flip = False
    CHIPNAME = 'MM2H'
    c7 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c7, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.2, 'mjjy': 0.2,'mjjx3':0.2, 'mjjy3':0.2}
    flip = False
    CHIPNAME = 'MM2I'
    c8 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c8, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.2, 'mjjy': 0.2,'mjjx3':0.2, 'mjjy3':0.2}
    flip = False
    CHIPNAME = 'MM2J'
    c9 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c9, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.22, 'mjjy': 0.22,'mjjx3':0.22, 'mjjy3':0.22}
    flip = False
    CHIPNAME = 'MM2K'
    SIPF_correction = 8100
    c10 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c10, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.22, 'mjjy': 0.22,'mjjx3':0.22, 'mjjy3':0.22}
    flip = False
    CHIPNAME = 'MM2L'
    c11 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c11, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.14, 'mjjy': 0.14,'mjjx3':0.14, 'mjjy3':0.14}
    flip = False
    CHIPNAME = 'MM2M'
    c12 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c12, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.155, 'mjjy': 0.155,'mjjx3':0.155, 'mjjy3':0.155}
    flip = False
    CHIPNAME = 'MM2N'
    c13 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c13, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.17, 'mjjy': 0.17,'mjjx3':0.17, 'mjjy3':0.17}
    flip = False
    CHIPNAME = 'MM2O'
    SIPF_correction = 8300
    c14 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c14, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.185, 'mjjy': 0.185,'mjjx3':0.185, 'mjjy3':0.185}
    flip = False
    CHIPNAME = 'MM2P'
    c15 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c15, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.2, 'mjjy': 0.2,'mjjx3':0.2, 'mjjy3':0.2}
    flip = False
    CHIPNAME = 'MM2Q'
    c16 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c16, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.22, 'mjjy': 0.22,'mjjx3':0.22, 'mjjy3':0.22}
    flip = False
    CHIPNAME = 'MM2R'
    c17 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c17, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = True
    CHIPNAME = 'MM2S'
    c18 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c18, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = True
    CHIPNAME = 'MM2T'
    c19 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c19, d, junc_correction, flip)

    # --------------------------#

    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = True
    CHIPNAME = 'MM2U'
    c20 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c20, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = True
    CHIPNAME = 'MM2V'
    c21 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c21, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = True
    CHIPNAME = 'MM2W'
    SIPF_correction = 8150
    c22 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c22, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = True
    CHIPNAME = 'MM2X'
    c23 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c23, d, junc_correction, flip)

    # --------------------------#
    junc_correction = {'mjjx': 0.15, 'mjjy': 0.16,'mjjx3':0.31, 'mjjy3':0.31}
    flip = False
    CHIPNAME = 'MM2Y'
    c24 = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
               chip_id_loc=(200, 100), textsize=(70, 70), two_layer=two_layer, solid=solid)
    draw_bridge_junction(c24, d, junc_correction, flip)

    # --------------------------#

    ### ===================== ###
    # Last row: JPA transformers
    ### --------------------- ###
    #####################################
    global sign
    global draw_len
    # common features for eps = 10.4, 7 GHz transformers
    cpw_rad1 = 400
    cpw_rad2 = 1.0 * cpw_rad1

    cpw_pinw = 10.0
    eps_eff = (1. + 10.4) / 2.
    cpw_gapw = calculate_gap_width(eps_eff, 50, cpw_pinw)
    cpw_pinw1 = 54.2
    cpw_gapw1 = 10.2
    cpw_pinw2 = 10.0
    cpw_gapw2 = 10.2
    max_len = 3200
    lambda2_len = 9000
    lambda4_len = 4500

    # --------------------------------#
    # 7GHz, eps = 10.4
    junc_l = 2.1
    CHIPNAME = 'ZQL'
    ct = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
              chip_id_loc=(300, 100), textsize=(70, 70), two_layer=two_layer, solid=solid, layer='0')
    # draw_launchers(ct, d, exclude=[0, 1, 2, 4, 5, 6, 7, 9])
    # draw_square_alignment_marks(ct, d)
    #
    # draw_test_full_transformer(ct, cpw_pinw, cpw_gapw, cpw_pinw1, cpw_pinw2, cpw_gapw1, cpw_gapw2, cpw_rad1, cpw_rad2,
    #                       lambda4_len, lambda2_len, max_len)
    # --------------------------------#

    if draw_layer1 == 0 and draw_layer2 == 1:
        label = False
    else:
        label = True

    # m.add_chip(ct, 5, label=label)
    m.add_chip(c0, 1, label=label)
    m.add_chip(c1, 1, label=label)
    m.add_chip(c2, 1, label=label)
    m.add_chip(c3, 1, label=label)
    m.add_chip(c4, 1, label=label)
    m.add_chip(c5, 1, label=label)
    m.add_chip(c6, 1, label=label)
    m.add_chip(c7, 1, label=label)
    m.add_chip(c8, 1, label=label)
    m.add_chip(c9, 1, label=label)
    m.add_chip(c10, 1, label=label)
    m.add_chip(c11, 1, label=label)
    m.add_chip(c12, 1, label=label)
    m.add_chip(c13, 1, label=label)
    m.add_chip(c14, 1, label=label)
    m.add_chip(c15, 1, label=label)
    m.add_chip(c16, 1, label=label)
    m.add_chip(c17, 1, label=label)
    # m.add_chip(c18, 1, label=label)
    # m.add_chip(c19, 1, label=label)
    # m.add_chip(c20, 1, label=label)
    # m.add_chip(c21, 1, label=label)
    # m.add_chip(c22, 1, label=label)
    # m.add_chip(c23, 1, label=label)
    # m.add_chip(c24, 1, label=label)
    # m.add_chip(c25, 1, label=label)
    # m.add_chip(c26, 1, label=label)
    # m.add_chip(c27, 1, label=label)

    return m


if __name__ == "__main__":

    m = output_wafer()
    m.save()

    print("\n\n Chip names are:")
    print("_____________________")
    for name in chip_names:
        print(name)
    print("_____________________\n\n")

    sleep(.1)

    if show_structure:
        if open_dwgviewer:
            subprocess.Popen(
                r'"C:\Program Files\Autodesk\DWG TrueView 2019 - English\dwgviewr.exe" "' + os.getcwd() + '\\' + MaskName + '-' + '6000' + '.dxf" ')
        if open_klayout:
            subprocess.Popen(
                r'"C:\Users\slab\AppData\Roaming\KLayout\klayout_app.exe" "' + os.getcwd() + '\\' + MaskName + '-' + '6000' + '.dxf" -e')  # editor mode
    elif show_wafer:
        if open_dwgviewer:
            subprocess.Popen(
                r'"C:\Program Files\Autodesk\DWG TrueView 2019 - English\dwgviewr.exe" "' + os.getcwd() + '\\' + MaskName + '.dxf" ')
        if open_klayout:
            subprocess.Popen(
                r'"C:\Users\slab\AppData\Roaming\KLayout\klayout_app.exe" "' + os.getcwd() + '\\' + MaskName + '.dxf" -e')

    try:
        # save a text file with the dimensions
        text = "Specs for %s created at %s at %s\n\nTapered: \t\t%s\nTotal length: \t\t%s um\nTotal width: \t\t%s um\nFinger length: \t\t%s um\nFinger width: \t\t%s um\nFinger separation: \t%s um\nNumber of fingers: \t%s\nTaper length: \t\t%s um\nCPW length: \t\t%s um\nCPW center pin: \t%s um\nCPW gap width: \t\t%s um" \
               % (MaskName + '-' + c.name, today, time.strftime("%H:%M:%S"), str(taper), c.total_size, c.total_width,
                  finger_length, finger_width, finger_spacing, noof_fingers, taper_length, cpw_length, cpw_pinw,
                  cpw_pinw)

        text_file = open(MaskName + '-' + c.name + '-specs.txt', "w")
        text_file.write(text)
        text_file.close()

        print("\nFinished creating capacitor with following parameters:")
        print(text)
    except:
        pass