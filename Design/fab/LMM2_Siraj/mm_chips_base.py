# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16, 2025 

@author: eesh
"""

# from slab.circuits import *
# from slab.circuits.mp_components import *
from mask_maker import *

import os, time
import subprocess
from time import sleep
import numpy as np
from copy import deepcopy

square =0


class MMChipsBase(): 
    """
    author: Eesh Gupta
    This class is base class for designing 3D chips. 
    """
    def __init__(self, c = None):
        
        self.c = c
        # self.mask_defaults = self.set_mask_init()
        # self.chipInit(self.c)
        # mm_chips_base = MMChipsBase(c)
        # mm_chips_base.chipInit(c, defaults=d)

        pass

    def draw_alignment_marks(self, draw_marks = False):
        c = self.structure
        ## alignment boxes for e-beam
        self.draw_launchers(c, self.d, exclude=[2, 3, 8])
        self.draw_launchers(c, self.d, exclude=[2, 3, 8])
        if draw_marks:
            self.draw_square_alignment_marks(c, self.flip)

    def curve_corner(self, structure, xpos, ypos, radius, type):
        '''
        Corner type: 0:left top, 1:right top, 2:right bottom, 3:left bottom
        xpos, ypos: starting  position of the corner
        '''
        c2 = structure
        ref_pos = c2.last
        ref_direction = c2.last_direction
        square = 0
       

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
            self.CPWBendNew(c2, angle=90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
            self.CPWBendNew(c2, angle=-90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
            # print('radius', radius)
            # print('c2.last_direction', c2.last_direction)
            self.CPWBendNew(c2, angle=90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
            self.CPWBendNew(c2, angle=-90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
                    square=square)
        ########################################
        c2.last_direction = ref_direction
        c2.last = ref_pos

    def set_mask_init(self):
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
        d.solid = False
        return d


    def chipInit(self, c, defaults = None):
        """
        This makes the launch pads on the chip. Input is an object c, which is the chip.
        From the launch pads, we can make connections on the chip.
        """
        # The following creates 8 launch pads on the chip. There are 3 pads per side, two
        # in each of the corners and then one in the middle of each side.
        if defaults is None:
            defaults = self.mask_defaults

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
        # self.draw_launchers(self.structure, self.d, exclude=[2, 3, 8])


    def CPWBendNew(self, structure, angle=90, pinw=10.0, gapw=20.0, 
                   radius=50.0, polyarc=1, segments=4, square=0):
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

    @staticmethod
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

class MMChipsBaseUtils:
    @staticmethod
    def draw_launchers(c, d, exclude=[]):
        cpw_pinw = 10.0
        eps_eff = (1. + 10.4) / 2.
        cpw_gapw = calculate_gap_width(eps_eff, 50, cpw_pinw)


        for k in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
            if k in exclude:
                pass
            else:
                if k == k:
                    Launcher(vars(vars()['c'])['s%d' % k], pinw=cpw_pinw, gapw=cpw_gapw, launcher_pinw=150, launcher_gapw=75)
                else:
                    Launcher(vars(vars()['c'])['s%d' % k], pinw=cpw_pinw, gapw=cpw_gapw, launcher_pinw=300, launcher_gapw=150)

    @staticmethod
    def draw_chip_alignment_marks(solid, d, c):
        CrossShapeAlignmentMarks(Structure(c, start=(125, 125), direction=90, defaults=d), width=2,
                                 armlength=120, solid=solid, layer='0')
        CrossShapeAlignmentMarks(Structure(c, start=(c.size[0] - 125, c.size[1] - 125), direction=90, defaults=d),
                                 width=2, armlength=120, solid=solid, layer='0')
        CrossShapeAlignmentMarks(Structure(c, start=(c.size[0] - 125, 125), direction=90, defaults=d),
                                 width=2, armlength=120, solid=solid, layer='0')
        CrossShapeAlignmentMarks(Structure(c, start=(125, c.size[1] - 125), direction=90, defaults=d),
                                 width=2, armlength=120, solid=solid, layer='0')

    @staticmethod
    def draw_square_alignment_marks(structure, flip):
        c = structure

        if flip:
            direc = c.s2.last_direction
            c.s2.last_direction = 180 + c.s2.last_direction

            align_box = 80
            c.s2.last = (150, 2100 + 18100 - (8300 + align_box / 2))
            CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
            c.s2.last = (5350, 2100 + 18100 - (8300 + align_box / 2))
            CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)

            c.s2.last = (150, 2100 + 18100 - (700 + align_box / 2))
            CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
            c.s2.last = (5350, 2100 + 18100 - (700 + align_box / 2))
            CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
            c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

            align_box = 500
            c.s2.last = (150, 2100 + 18100 - (8300 + align_box / 2))
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (5350, 2100 + 18100 - (8300 + align_box / 2))
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

            c.s2.last = (150, 2100 + 18100 - (700 + align_box / 2))
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (5350, 2100 + 18100 - (700 + align_box / 2))
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

            c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)
            c.s2.last_direction = direc

        else:
            shift_y = 0
            align_box = 80
            c.s2.last = (200, 2100 + 50500 + align_box / 2 + shift_y)
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (5250, 2100 + 50500 + align_box / 2 + shift_y)
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

            c.s2.last = (200, 700 + align_box / 2)
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (5250, 700 + align_box / 2)
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

            align_box = 500
            c.s2.last = (200, 2100 + 50500 + align_box / 2 + shift_y)
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (5250, 2100 + 50500 + align_box / 2 + shift_y)
            CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
            c.s2.last = (1000, 9800 - 50)
            c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)



            

# Attach these as static methods to MMChipsBase for convenience
MMChipsBase.draw_launchers = staticmethod(MMChipsBaseUtils.draw_launchers)
MMChipsBase.draw_chip_alignment_marks = staticmethod(MMChipsBaseUtils.draw_chip_alignment_marks)
MMChipsBase.draw_square_alignment_marks = staticmethod(MMChipsBaseUtils.draw_square_alignment_marks)