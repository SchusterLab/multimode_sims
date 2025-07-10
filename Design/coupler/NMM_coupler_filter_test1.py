# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 09:54:15 2021

@author: ziqian
"""

from slab.circuits import *
from slab.circuits.mp_components import *

import os, time
import subprocess
from time import sleep
import numpy as np

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

        pass

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
            CPWBendNew(c2, angle=90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
            CPWBendNew(c2, angle=-90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
            CPWBendNew(c2, angle=90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
            CPWBendNew(c2, angle=-90, pinw=0, gapw=radius / 2, radius=radius / 2, polyarc=1, segments=4,
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
        d.solid = True
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


def draw_launchers(c, d, exclude=[]):
    cpw_pinw = 10.0
    cpw_gapw = calculate_gap_width(eps_eff, 50, cpw_pinw)

    mm_chips_base = MMChipsBase(c)
    mm_chips_base.chipInit(c, defaults=d)

    for k in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:

        if k in exclude:
            pass
        else:
            if k==k:
                Launcher(vars(vars()['c'])['s%d' % k], pinw=cpw_pinw, gapw=cpw_gapw, launcher_pinw=150, launcher_gapw=75)
            else:
                Launcher(vars(vars()['c'])['s%d' % k], pinw=cpw_pinw, gapw=cpw_gapw, launcher_pinw=300,
                         launcher_gapw=150)


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

        align_box = 80  # center is the reference
        c.s2.last = (150, 2100+18100-(8300 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
        c.s2.last = (5350, 2100+18100-(8300 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)

        c.s2.last = (150, 2100+18100-(700 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
        c.s2.last = (5350, 2100+18100-(700 + align_box / 2))
        CPWStraight(c.s2, 80, pinw=0, gapw=align_box / 2)
        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        #  draw negative marker shadow
        align_box = 500  # center is the reference
        c.s2.last = (150, 2100+18100-(8300 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (5350, 2100+18100-(8300 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

        c.s2.last = (150, 2100+18100-(700 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (5350, 2100+18100-(700 + align_box / 2))
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)


        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        c.s2.last_direction = direc

    else:
        shift_y = 5000*0
        align_box = 80  # center is the reference
        c.s2.last = (200, 2100+50500 + align_box / 2+shift_y)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (5250, 2100+50500 + align_box / 2+shift_y)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

        c.s2.last = (200, 700 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (5250, 700 + align_box / 2)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        #  draw negative marker shadow
        align_box = 500  # center is the reference
        c.s2.last = (200, 2100+50500 + align_box / 2+shift_y)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (5250, 2100+50500 + align_box / 2+shift_y)
        CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)

        # c.s2.last = (200, 700 + align_box / 2)
        # CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        # c.s2.last = (5250, 700 + align_box / 2)
        # CPWStraight(c.s2, align_box, pinw=0, gapw=align_box / 2)
        c.s2.last = (1000, 9800 - 50)
        # CPWStraight(c.s2, 100, pinw=0, gapw=1500 / 2)
        c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)




def CPWBendNew(structure, angle=90, pinw=10.0, gapw=20.0, radius=50.0, polyarc=1, segments=4, square=0):
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


class MultimodeBalun():
    def __init__(self, structure, xpos, ypos, shape):
        self.structure = structure
        self.xpos = xpos
        self.ypos = ypos
        self.shape = shape
        self.coords = {}  # key: coordinate identifier label, value: (x, y) tuple
    
    def draw_cpw(self, pin_width, gap_width, ground_width, length):
        """
        Draws a coplanar waveguide (CPW) structure with specified dimensions.
        Parameters:
            pin_width (float): Width of the center conductor (pin) of the CPW.
            gap_width (float): Width of the gap between the center conductor and the ground planes.
            ground_width (float): Width of the ground conductor regions.
            length (float): Length of the CPW section to be drawn.
        Notes:
            - The function draws the center conductor and ground conductors as separate CPW sections.
            - The CPWStraight function is used to create the straight segments of the CPW.
            - The gaps are drawn on either side of the pin as per the CPW geometry.
        """
        c = self.structure
        shape = self.shape
        ref4 = c.s2.last

        # center conductor 
        CPWStraight(c.s2, length,
                    pinw=0,
                    gapw=pin_width/2)
        c.s2.last = ref4
        # ground conductor
        CPWStraight(c.s2,length,
                    pinw=pin_width+2*gap_width,
                    gapw=ground_width)

    def draw_linear_tapered_cpw(self, start_pin_width, start_gap_width,
                                 start_ground_width, stop_pin_width,
                                   stop_gap_width, stop_ground_width, length):
        """
        Draws a linear tapered coplanar waveguide (CPW) section with specified start and stop dimensions.
        This method creates a CPW structure where the pin (center conductor), gap, and ground widths
        linearly taper from their starting values to their stopping values over the given length.
        Parameters:
            start_pin_width (float): Width of the center conductor at the start of the taper.
            start_gap_width (float): Width of the gap between the center conductor and ground at the start.
            start_ground_width (float): Width of the ground plane at the start of the taper.
            stop_pin_width (float): Width of the center conductor at the end of the taper.
            stop_gap_width (float): Width of the gap between the center conductor and ground at the end.
            stop_ground_width (float): Width of the ground plane at the end of the taper.
            length (float): Length of the tapered CPW section.
        Returns:
            None
        """
        c = self.structure
        ref = c.s2.last
        # center conductor 
        CPWLinearTaper(c.s2, length, 
                       0, 0, 
                       start_pin_width/2, stop_pin_width / 2) # center conductor
        c.s2.last = ref
        # ground 
        CPWLinearTaper(c.s2, length,
                        start_pin_width+2*start_gap_width, stop_pin_width+2*stop_gap_width,
                         start_ground_width, stop_ground_width)
    
    def draw_port3(self):
        """
        Draws the geometry for port 3 of the balun.

        This method constructs the CPW and ground structures for port 3, starting from the balun center.
        It draws a straight CPW section, a linear tapered CPW, a CPW with ground, and a stripline section,
        updating the structure's state as it proceeds.

        Args:

        Side Effects:
            - Modifies the structure by adding CPW and ground segments for port 3.
            - Updates the internal coordinates dictionary with relevant reference points.
        """
        # port 3
        shape = self.shape
        port3_ref = self.coords['balun_center']
        l = shape['cpw_pin']
        c = self.structure
        c.s2.last = (port3_ref[0], port3_ref[1] - shape['cpw_pin'] / np.sqrt(3))
        c.s2.last_direction = 150
        c.s2.move(shape['cpw_pin'] / 2 / np.sqrt(3))
        # from center to CPS section (just the pin)
        CPWStraight(c.s2, l / 2 * np.sqrt(3), pinw=0, gapw=l / 2)
        # draws tapered CPW
        ref3 = c.s2.last
        self.draw_linear_tapered_cpw(
            start_pin_width=shape['cpw_pin'],
            start_gap_width=shape['cpw_gap'],
            start_ground_width=0,
            stop_pin_width=shape['cpw_pin'],
            stop_gap_width=shape['cpw_gap'],
            stop_ground_width=shape['cpw_ground_width'],
            length=shape['cpw_ground_width'] / np.sqrt(3)
        )
        # Draws CPW which consists of one CPWStraight for center conductor and one CPWStraight for ground
        self.draw_cpw(
            length=shape['cps_port2_length'] + shape['cpw_ground_width'] / np.sqrt(3) * 2,
            pin_width=l,
            gap_width=shape['cpw_gap'],
            ground_width=shape['cpw_ground_width']
        )
        # stripline
        CPWStraight(
            c.s2,
            shape['cps_w'] / 1.5,
            pinw=l + shape['cpw_gap'] * 2,
            gapw=shape['cpw_ground_width']
        )
        # closing the stripline
        CPWStraight(
            c.s2,
            shape['cps_w'] / 1.5,
            pinw=0,
            gapw=shape['cpw_ground_width'] + (l + shape['cpw_gap'] * 2) / 2
        )

    def draw_port5(self):
        """
        Draws port 5 of the coupler structure using a series of CPW (coplanar waveguide) and tapering operations.
        This method performs the following steps:
        1. Sets the starting position and direction for port 5 based on the reference from port 3 and the CPW pin width.
        2. Moves the drawing cursor to the correct location for the start of port 5.
        3. Draws a straight CPW segment with specified pin and gap widths.
        4. Draws a linear tapered CPW section, transitioning ground width from 0 to the specified value.
        5. Draws a standard CPW section with given pin, gap, and ground widths for the port.
        6. Closes the CPW by drawing a final straight segment with adjusted gap width.
        Assumes the following instance variables and arguments are available:
            - c: An object containing the drawing state, including 's2' (the drawing cursor).
            - port3_ref: Tuple (x, y) reference position for port 3.
            - shape: Dictionary containing geometric parameters such as 'cpw_pin', 'cpw_gap', 'cpw_ground_width', and 'cps_port2_length'.
            - l: Length parameter used for gap width calculations.
            - CPWStraight: Function or class for drawing straight CPW segments.
            - self.draw_linear_tapered_cpw: Method for drawing a tapered CPW section.
            - self.draw_cpw: Method for drawing a standard CPW section.
        Returns:
            None
        """
        shape = self.shape
        port3_ref = self.coords['balun_center']
        l = shape['cpw_pin']
        c = self.structure

        # port 5
        c.s2.last = (port3_ref[0], port3_ref[1] - shape['cpw_pin'] / np.sqrt(3))
        # c.s2.last = port3_ref
        c.s2.last_direction = 30
        c.s2.move(shape['cpw_pin'] / 2 / np.sqrt(3))
        CPWStraight(c.s2, shape['cpw_pin'] / 2 * np.sqrt(3), pinw=0,
                    gapw=l / 2)

        ref3 = c.s2.last
        self.draw_linear_tapered_cpw( start_pin_width=shape['cpw_pin'],
                                     stop_pin_width=shape['cpw_pin'],
                                     start_gap_width=shape['cpw_gap'],
                                     stop_gap_width=shape['cpw_gap'],
                                     start_ground_width=0,
                                     stop_ground_width=shape['cpw_ground_width'],
                                     length=shape['cpw_ground_width'] / np.sqrt(3))
        ## cpw 
        self.draw_cpw(pin_width=shape['cpw_pin'],
                      gap_width=shape['cpw_gap'],
                      ground_width=shape['cpw_ground_width'],
                      length=shape['cps_port2_length'] + shape['cpw_ground_width'] / np.sqrt(3) * 2
                      )
        ## closing the cpw
        CPWStraight(c.s2, shape['cpw_ground_width']/1.5, pinw=0,
                    gapw=shape['cpw_ground_width'] + (l + shape['cpw_gap'] * 2) / 2)
        
    def draw_port1(self):
        """
        Draws the geometry for port 1 of the coupler, including the CPW sections, a tapered section, 
        and the balun center pin connection.
        The method performs the following steps:
            1. Draws the initial CPW section using parameters from the `shape` dictionary.
            2. Draws a second CPW section with different ground width and length.
            3. Stores the end coordinate of the second CPW section as 'balun_start_port1'.
            4. Draws a linear tapered CPW section starting from 'balun_start_port1', tapering the ground width to zero.
            5. Draws the balun center pin connection from port 1 using a straight CPW and a linear taper.
        Assumes that `self.shape` contains all necessary geometric parameters, and that `self.structure` 
        and `self.coords` are properly initialized.
        Side Effects:
            - Modifies `self.coords` by adding 'balun_start_port1'.
            - Draws on the structure `self.structure` using various CPW drawing methods.
        """
        
        shape = self.shape
        c = self.structure

        self.draw_cpw(pin_width=shape['cpw_pin'],
                      gap_width=shape['cpw_gap'],
                      ground_width=shape['cpw_ground_width'],
                      length=shape['scpw_length'] )
        
        self.draw_cpw(pin_width=shape['cpw_pin'],
                      gap_width=shape['cpw_gap'],
                      ground_width=shape['cpw_ground4'],
                      length= shape['cpw_ground4_length'])
        self.coords['balun_start_port1'] = c.s2.last


        ## tapered section
        ref3 = self.coords['balun_start_port1']
        self.draw_linear_tapered_cpw(start_pin_width=shape['cpw_pin'],
                                     start_gap_width=shape['cpw_gap'],
                                     start_ground_width=shape['cpw_ground4'],
                                     stop_pin_width=shape['cpw_pin'],
                                     stop_gap_width=shape['cpw_gap'],
                                     stop_ground_width=0,
                                     length=shape['cpw_ground4_length']/np.sqrt(3))

        ## balun center pin connection from port 1  
        cpw_gap = shape['cpw_gap']            
        CPWStraight(c.s2, (shape['cps_gap2']+cpw_gap) / np.sqrt(3), pinw=0,
                    gapw=shape['cpw_pin'] / 2)
        CPWLinearTaper(c.s2, shape['cpw_pin'] / 2 * np.sqrt(3), 0,
                    0,
                    shape['cpw_pin'] / 2, 0)
       
    def draw(self):
        c = self.structure
        shape = self.shape
        xpos, ypos = self.xpos, self.ypos

        pin_ref = c.s2.last
        pin_direction = c.s2.last_direction
        ################################################
        c.s2.last = (xpos, ypos)
        eps_eff = (1. + 10.4) / 2.
        cpw_gap = shape['cpw_gap']

        initial_length = 600-200
        c.s2.last_direction = 90

        # draw launcher part
        CPWStraight(c.s2, initial_length, pinw=0, gapw=shape['ground']) # ground
        CPWStraight(c.s2, shape['gap2ground'], pinw=shape['launcher_width']*2,  # 2 because total gap and launcher are together twice launcher width
                    gapw=shape['ground']-shape['launcher_width'])
        
        self.coords['begin_launcher'] = c.s2.last
        ## The launcher 
        self.draw_cpw(length = shape['launcher_length'],
            pin_width = shape['launcher_width'],
            gap_width = shape['launcher_width']/2,
            ground_width = shape['ground']-shape['launcher_width'])
        # c.s2.last = self.coords['begin_launcher']
        self.draw_linear_tapered_cpw(start_pin_width=shape['launcher_width'],
                                     start_gap_width=shape['launcher_width']/2,
                                     start_ground_width=shape['ground']-shape['launcher_width'],
                                     stop_pin_width=shape['cpw_pin'],
                                     stop_gap_width=shape['cpw_gap'],
                                     stop_ground_width=shape['cpw_ground_width'],
                                     length=shape['launcher_taper'])
        self.coords['end_launcher'] = c.s2.last
        

        # draw port 1 section
        self.draw_port1()
        port3_ref = c.s2.last
        self.coords['balun_center'] = port3_ref


        # port 3 
        self.draw_port3()

        # port 5
        self.draw_port5()

        # port 2
        l = shape['cps_gap2']
        # c.s2.last = (port2_ref[0] - (cpw_gap * 2 + shape['cpw_pin']) / 2 - shape['cpw_ground_width'] - l / 2 / 2,
        #              port2_ref[1] + l / 2 / 2 * np.sqrt(3))
        c.s2.last = port3_ref
        c.s2.last_direction = 270
        c.s2.move(np.sqrt(3) / 3 * shape['cpw_pin'])
        c.s2.last_direction = 210
        c.s2.move(shape['cpw_ground_width'] * 2 / np.sqrt(3) + np.sqrt(3) / 3 * shape['cpw_pin']+np.sqrt(3)/2*shape['cpw_pin'])
        CPWLinearTaper(c.s2, shape['cps_w'] / np.sqrt(3)-0.0001*0, l,
                    l, 0, shape['cps_w'])
        CPWStraight(c.s2, shape['cps_port2_length'], pinw=l,
                    gapw=shape['cps_w'])

        CPWStraight(c.s2, shape['cps_w'], pinw=0,
                    gapw=shape['cps_w']+l/2)

        # port 6
        l = shape['cps_gap2']
        # c.s2.last = (port2_ref[0] - (cpw_gap * 2 + shape['cpw_pin']) / 2 - shape['cpw_ground_width'] - l / 2 / 2,
        #              port2_ref[1] + l / 2 / 2 * np.sqrt(3))
        c.s2.last = port3_ref
        c.s2.last_direction = 270
        c.s2.move(np.sqrt(3) / 3 * shape['cpw_pin'])
        c.s2.last_direction = 330
        c.s2.move(shape['cpw_ground_width'] * 2 / np.sqrt(3) + np.sqrt(3) / 3 * shape['cpw_pin'] + np.sqrt(3) / 2 * shape[
            'cpw_pin'])
        CPWLinearTaper(c.s2, shape['cps_w'] / np.sqrt(3)-0.0001*0, l,
                    l, 0, shape['cps_w'])
        CPWStraight(c.s2, shape['cps_port2_length'], pinw=l,
                    gapw=shape['cps_w'])

        # port 4 CPS out
        c.s2.last = port3_ref
        c.s2.last_direction = 90
        c.s2.move(shape['cpw_ground_width']*2/np.sqrt(3)+np.sqrt(3)/2*shape['cpw_pin'])
        CPWLinearTaper(c.s2, shape['cps_w'] / np.sqrt(3), shape['cps_gap2'],
                    shape['cps_gap2'], 0, shape['cps_w'])
        CPWStraight(c.s2, 500, pinw=shape['cps_gap2']+0.0015*0,
                    gapw=shape['cps_w'])


        # linear taper
        CPWLinearTaper(c.s2, shape['cps_pin_taper'], shape['cps_gap2'],
                    shape['cps_pin_target'], shape['cps_w'], shape['cps_gap_target'])





        reff = c.s2.last


    ####################################################
        c.s2.last = pin_ref
        c.s2.last_direction = pin_direction

        return reff

class MultimodeCouplerPad(MMChipsBase):
    def __init__(self, structure, xpos, ypos, shape):
        super().__init__()
        self.structure = structure
        self.xpos = xpos
        self.ypos = ypos
        self.shape = shape
        self.coords = {} # key: coordinate identifier label, value: (x, y) tuple
        
    
    def draw_bar(self, coupler_center_y, coupler_gap,  thin_bar_length): 

        movey = coupler_center_y
        gap_squid = coupler_gap

        c = self.structure
        xpos = self.xpos
        ypos = self.ypos
        shape = self.shape
        thinner = thin_bar_length
        fat_bar_length = shape['bar_length'] - gap_squid / 2 - thinner

        # Left Bar 
        c.s2.last_direction = 180
        c.s2.last = (xpos - gap_squid / 2, ypos+ movey)
        self.coords['left_bar_start'] = c.s2.last # save 7

        # thinner bar ; here the 2 gaps in CPW is the metal pad; that's why pin width is 0 and gap width is dquid size/2
        CPWStraight(c.s2, thinner, pinw=0, gapw=shape['squid_y'] / 2)
        self.coords['left_bar_thin_end'] = c.s2.last # save 3

        # fatter bar
        CPWStraight(c.s2, fat_bar_length, pinw=0, gapw=shape['bar_thick']/2)
        self.coords['left_bar_thick_end'] = c.s2.last 

        # bar towards right
        c.s2.last_direction = 0 
        c.s2.last = (xpos + gap_squid / 2, ypos+ movey)
        self.coords['right_bar_start'] = c.s2.last  # save 5
        CPWStraight(c.s2, thinner, pinw=0, gapw=shape['squid_y'] / 2)
        
        self.coords['right_bar_thin_end'] = c.s2.last # save 4
        CPWStraight(c.s2, fat_bar_length, pinw=0, gapw=shape['bar_thick'] / 2)
        self.coords['right_bar_thick_end'] = c.s2.last
 

    def draw_left_pad(self, width_of_long_left_pad): 
        """
        Draws the left pad of the structure by creating two consecutive CPW (coplanar waveguide) straight sections:
        a long section and a short section. The method updates the structure's coordinates accordingly.
        Args:
            width_of_long_left_pad (float): The length of the long section of the left pad.
        Side Effects:
            - Modifies the structure's state by drawing two CPW straight sections.
            - Updates the 'left_pad_start' and 'left_pad_end' entries in the self.coords dictionary.
        Notes:
            - The width and gap of each section are determined by the 'shape' dictionary.
            - The starting position and direction are set based on the structure's current state.
        """
        
        shape = self.shape
        c = self.structure
        xpos = self.xpos
        ypos = self.ypos


        for_diff = width_of_long_left_pad
        c.s2.last_direction = 180
        c.s2.last = (self.coords['left_bar_thick_end'][0], ypos )
        self.coords['left_pad_start'] = c.s2.last # save 6
        # long section of left pad
        CPWStraight(c.s2, for_diff, pinw=0, gapw=shape['right_pady'] / 2)
        self.coords['left_pad_long_end'] = c.s2.last 
        
        # short section of left pad
        CPWStraight(c.s2, shape['left_padx']-for_diff, pinw=0, gapw=shape['left_pady'] / 2)
        self.coords['left_pad_end'] = c.s2.last # save 1
    
    def draw_right_pad(self):
        """
        Adds a right pad to the structure at the specified position.
        This method updates the structure by extending it with a right pad, starting from the end of the right bar.
        It uses the specified pad dimensions from the `shape` attribute and updates the coordinates dictionary
        with the end position of the right pad.
        Side Effects:
            - Modifies `self.structure` by adding a CPWStraight segment.
            - Updates `self.coords['right_pad_end']` with the new end position.
        Attributes Used:
            - self.structure: The structure object to which the pad is added.
            - self.xpos, self.ypos: The current x and y positions.
            - self.shape: Dictionary containing pad dimensions.
            - self.coords: Dictionary for storing coordinate points.
        """

        c = self.structure
        xpos = self.xpos
        ypos = self.ypos
        shape = self.shape
        
        c.s2.last_direction = 0
        c.s2.last = (self.coords['right_bar_thick_end'][0], ypos)
        CPWStraight(c.s2, shape['right_padx'], pinw=0, gapw=shape['right_pady'] / 2)
        self.coords['right_pad_end'] = c.s2.last # save2 = c.s2.last

    def curve_bar_corners(self, curve=5):
        """
        Draws curved corners for the bars of the multimode coupler pad.

        Args:
            curve (float): The radius of curvature for the corners.

        This method draws curved connections between the thin and thick bar sections,
        as well as at the ends of the bars, using the `curve_corner` function.
        """
        c = self.structure
        shape = self.shape

        # left bar thick
        save3 = self.coords['left_bar_thin_end']
        self.curve_corner(c.s2, save3[0], save3[1] + shape['bar_thick'] / 2, curve * 2, 2)
        self.curve_corner(c.s2, save3[0], save3[1] - shape['bar_thick'] / 2, curve * 2, 1)

        # right bar thick
        save4 = self.coords['right_bar_thin_end']
        self.curve_corner(c.s2, save4[0], save4[1] + shape['bar_thick'] / 2, curve, 3)
        self.curve_corner(c.s2, save4[0], save4[1] - shape['bar_thick'] / 2, curve, 0)

        # left bar thin to thick connection
        save7 = save3
        self.curve_corner(c.s2, save7[0], save7[1] + shape['squid_y'] / 2, curve, 0)
        self.curve_corner(c.s2, save7[0], save7[1] - shape['squid_y'] / 2, curve, 3)

        # right bar thin to thick connection
        save5 = self.coords['right_bar_thin_end']
        self.curve_corner(c.s2, save5[0], save5[1] + shape['squid_y'] / 2, curve, 1)
        self.curve_corner(c.s2, save5[0], save5[1] - shape['squid_y'] / 2, curve, 2)

        # left bar thick to pad connection
        save = self.coords['left_bar_thick_end']
        self.curve_corner(c.s2, save[0], save[1] + shape['bar_thick'] / 2, curve, 0)
        self.curve_corner(c.s2, save[0], save[1] - shape['bar_thick'] / 2, curve, 3)

        # right bar thick to pad connection
        save6 = self.coords['right_bar_thick_end']
        self.curve_corner(c.s2, save6[0], save6[1] + shape['bar_thick'] / 2, curve, 1)
        self.curve_corner(c.s2, save6[0], save6[1] - shape['bar_thick'] / 2, curve, 2)

    def curve_left_pad_corners(self, curve=5):
        """
        Draws curved corners for the left pad of the multimode coupler pad.

        Args:
            curve (float): The radius of curvature for the corners.

        This method draws curved connections at the ends and junctions of the left pad,
        using the `curve_corner` function.
        """
        c = self.structure
        shape = self.shape

        # left pad leftmost corners
        save1 = self.coords['left_pad_end']
        self.curve_corner(c.s2, save1[0], save1[1] + shape['left_pady'] / 2, curve * 2, 3)
        self.curve_corner(c.s2, save1[0], save1[1] - shape['left_pady'] / 2, curve * 2, 0)

        # left pad short to long connection
        save = self.coords['left_pad_long_end']
        self.curve_corner(c.s2, save[0], save[1] + shape['left_pady'] / 2, curve * 20, 1)
        self.curve_corner(c.s2, save[0], save[1] - shape['left_pady'] / 2, curve * 20, 2)

        # long left pad corners near connection to short left pad
        self.curve_corner(c.s2, save[0], save[1] + shape['right_pady'] / 2, curve * 2, 3)
        self.curve_corner(c.s2, save[0], save[1] - shape['right_pady'] / 2, curve * 2, 0)

        # long left pad near left bar thick end connection
        save6 = self.coords['left_pad_start']
        self.curve_corner(c.s2, save6[0], save6[1] + shape['right_pady'] / 2, curve * 2, 2)
        self.curve_corner(c.s2, save6[0], save6[1] - shape['right_pady'] / 2, curve * 2, 1)

    def curve_right_pad_corners(self, curve=5):
        """
        Draws curved corners for the right pad of the multimode coupler pad.

        Args:
            curve (float): The radius of curvature for the corners.

        This method draws curved connections at the ends and junctions of the right pad,
        using the `curve_corner` function.
        """
        c = self.structure
        shape = self.shape

        save2 = self.coords['right_pad_end']
        self.curve_corner(c.s2, save2[0], save2[1] + shape['right_pady'] / 2, curve * 2, 2)
        self.curve_corner(c.s2, save2[0], save2[1] - shape['right_pady'] / 2, curve * 2, 1)

        self.curve_corner(c.s2, save2[0] - shape['right_padx'], save2[1] + shape['right_pady'] / 2, curve * 2, 3)
        self.curve_corner(c.s2, save2[0] - shape['right_padx'], save2[1] - shape['right_pady'] / 2, curve * 2, 0)


    def draw(self):
        c = self.structure
        shape = self.shape
        xpos = self.xpos
        ypos = self.ypos

        pin_ref = c.s2.last
        pin_direction = c.s2.last_direction
        ################################################
        c.s2.last = (xpos, ypos)
        movey = 813-250-1400 # coupler center position
        curve = 5

        if shape['tripplr_jj']:
            gap_squid = shape['squid_x']*3
        else:
            gap_squid = shape['squid_x']
        
    
        thin_bar_length =200
        fat_bar_length = shape['bar_length'] - gap_squid / 2 - thin_bar_length
        width_long_left_pad = 120 # this is the width of long section of left pad (length in terms of drawn CPW)
        
        self.draw_bar(coupler_center_y=movey, coupler_gap=gap_squid, thin_bar_length=thin_bar_length)
        
        self.draw_left_pad(width_of_long_left_pad=width_long_left_pad)
        self.draw_right_pad()

        # curve corners
        self.curve_bar_corners(curve=curve)
        self.curve_left_pad_corners(curve=curve)
        self.curve_right_pad_corners(curve=curve)

        
    ####################################################
        c.s2.last = pin_ref
        c.s2.last_direction = pin_direction


class MultimodeFluxLine(MMChipsBase):
    # In this class, the CPW pin is the gap between striplines and CPW gap is the striplines
    def __init__(self, structure, xpos, ypos, shape):
        """
        Initialize the object with the given structure, position, and shape.
        Args:
            structure: The structure object to which this instance is associated.
            xpos (float): The x-coordinate position.
            ypos (float): The y-coordinate position.
            shape: The geometric shape or profile for the object. Contains parameters like lengths, widths, and gaps for the multimode coupler.
        Attributes:
            structure: Stores the provided structure object.
            xpos (float): Stores the x-coordinate.
            ypos (float): Stores the y-coordinate.
            shape: Stores the shape information.
            pin_ref: Reference to the last pin in structure.s2.
            pin_direction: Direction of the last pin in structure.s2.
            round_c (int): Radius for rounding corners of high and low impedance sections.
        """
        super().__init__()
        self.structure = structure
        self.xpos = xpos
        self.ypos = ypos
        self.shape = shape

        c = self.structure
        self.pin_ref = c.s2.last
        self.pin_direction = c.s2.last_direction
        ################################################
        self.round_c = 25+25*2 # radius for rounding corners of high and low impedance sections
        

    def draw_low_impedance_pads(self, last_point_coord):
        """
        Draws low impedance pads starting from the specified coordinate.
        This method updates the structure's last point and direction, then creates a 
        low impedance coplanar waveguide (CPW) section with specified length, width, 
        and gap parameters from the shape configuration. It returns the coordinate 
        of the last point after drawing the pads.
        Args:
            last_point_coord (tuple or list): The (x, y) coordinate to start drawing the pads from.
        Returns:
            tuple: The (x, y) coordinate of the last point after drawing the low impedance pads.
        """

        self.structure.s2.last = last_point_coord
        self.structure.s2.last_direction = 90
        target = self.shape['length_lo']
        width = self.shape['w_lo']
        sgap = self.shape['s_lo']
        CPWStraight(self.structure.s2, target, pinw=sgap, gapw=width)
        return self.structure.s2.last
    
    def draw_high_impedance_section(self, last_point_coord):
        '''
        Goal: Draw the S shape of the high impedance section of a multimode coupler.
        Draws the high impedance section of a multimode coupler starting from the specified coordinate.
        This method constructs a complex high impedance coplanar waveguide (CPW) section by sequentially adding
        straight and bent segments, using parameters defined in the object's shape configuration. The sequence 
        includes multiple straight sections and bends, with careful length accounting to match the target 
        high impedance section length. The method updates the structure's state as it draws, and returns the 
        final coordinate after completing the section.
            last_point_coord (tuple or list): The (x, y) coordinate from which to start drawing the high impedance section.
        Notes:
            - The method relies on several external functions/classes (e.g., CPWStraight, CPWBendNew) and 
              attributes (e.g., self.structure, self.shape).
            - The lengths and radii for the segments are calculated based on both hardcoded values and 
              configuration parameters.
        '''
        c = self.structure
        shape = self.shape

        section_length = 940 - 85 + 320 - 120 + 100+25+25+40 # large length/2 (should be roughly the size of low impedance section / size of 1 pad)
        section1_length = 180 - 70-15
        c.s2.last = last_point_coord
        c.s2.last_direction = 90
        radius = 100-10-10
        target = shape['length_hi'] # total length of the high impedance section
        sec1 = section1_length
        width = shape['w_hi']
        sgap = shape['s_hi']

        # 1. Towards the high impedance section
        CPWStraight(c.s2, sec1, pinw=sgap, gapw=width)
        remain = target - sec1 # remaining length after the first straight section

        # 2. First bend into the high impedance section
        CPWBendNew(c.s2, angle=-90, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi/2*radius # subtract the length of the bend

        # 3. First Straight section in the high impedance section
        sec2 = section_length
        CPWStraight(c.s2, sec2, pinw=sgap, gapw=width)
        remain = remain - sec2

        # 4. Second bend in the high impedance section
        CPWBendNew(c.s2, angle=180, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi * radius

        # 5. Second Straight section in the high impedance section (the largest)
        sec3 = section_length*2+radius*2 # larg length (2x low impedance pad length)
        CPWStraight(c.s2, sec3, pinw=sgap, gapw=width)
        remain = remain - sec3

        # 6. Third bend in the high impedance section
        CPWBendNew(c.s2, angle=-180, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi * radius

        # 7. Third Straight section in the high impedance section
        CPWStraight(c.s2, sec2, pinw=sgap, gapw=width)
        remain = remain - sec2
        CPWBendNew(c.s2, angle=90, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi / 2 * radius
        # 8. Final Straight section in the high impedance section
        CPWStraight(c.s2, remain, pinw=sgap, gapw=width)
        return c.s2.last

        """def curve_corner(self, structure, xpos, ypos, radius, type):
        '''
        Corner type: 0:left top, 1:right top, 2:right bottom, 3:left bottom
        xpos, ypos: starting  position of the corner
        """

    def draw_flux_loop(self, centeral_coords, loop_curve_radius, square = 0):
        """
        Draws a flux loop structure for the multimode coupler.

        Args:
            centeral_coords (tuple): The (x, y) coordinates at the center of the flux loop.
            loop_curve_radius (float): The radius for the curved corners of the flux loop.
            square (int, optional): If set, draws square bends instead of rounded. Default is 0.
        Notes:
            - Use the same sgap and width as the high impedance section.
        This method draws the right and left halves of the flux loop, each consisting of a sequence of bends and straight sections,
        using the provided geometric parameters. The structure's state is updated as the loop is drawn.
        """
        c = self.structure
        shape = self.shape
        loop_curve = loop_curve_radius
        width = shape['w_hi']
        sgap = shape['s_hi']

        # Save the starting point for the right half of the loop
        right_start = (centeral_coords[0] - sgap / 2 - width / 2, centeral_coords[1])
        c.s2.last_direction = 90
        c.s2.last = right_start

        # Right half of the flux loop
        CPWBendNew(c.s2, angle=90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_x'] / 2 - loop_curve * 2 - sgap / 2, pinw=0, gapw=width / 2)
        CPWBendNew(c.s2, angle=-90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_y'] - loop_curve * 2, pinw=0, gapw=width / 2)
        CPWBendNew(c.s2, angle=-90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_x'] - loop_curve * 2 + width, pinw=0, gapw=width / 2)

        # Save the starting point for the left half of the loop
        left_start = (centeral_coords[0] + sgap / 2 + width / 2, centeral_coords[1])
        c.s2.last_direction = 90
        c.s2.last = left_start

        # Left half of the flux loop
        CPWBendNew(c.s2, angle=-90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_x'] / 2 - loop_curve * 2 - sgap / 2, pinw=0, gapw=width / 2)
        CPWBendNew(c.s2, angle=90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_y'] - loop_curve * 2, pinw=0, gapw=width / 2)
        CPWBendNew(c.s2, angle=90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        

    def round_corners_low_to_high(self, central_coords, radius):
        """
        Rounds the corners between low and high impedance sections of the multimode coupler.

        Args:
            central_coords (tuple): The (x, y) coordinates at the center of the transition.
            radius (float): The radius for the rounded corners.

        This method draws the rounded corners connecting the low and high impedance sections
        of the multimode coupler, as well as the outer corners of the low impedance section.
        It uses the `curve_corner` function to create the appropriate corner geometry.
        """
        c = self.structure
        shape = self.shape
        sgap = shape['s_lo']
        # save2 is the central coordinate for the transition
        save2 = central_coords
        round_c = radius

        # Connection between low and high impedance sections
        self.curve_corner(c.s2, save2[0] - sgap / 2 - shape['w_hi'], save2[1], round_c, 1)  # left corner
        self.curve_corner(c.s2, save2[0] + sgap / 2 + shape['w_hi'], save2[1], round_c, 0)  # right corner

        # Outer corners of low impedance section
        self.curve_corner(c.s2, save2[0] - sgap / 2 - shape['w_lo'], save2[1], round_c, 3)
        self.curve_corner(c.s2, save2[0] + sgap / 2 + shape['w_lo'], save2[1], round_c, 2)

    def round_corners_high_to_low(self, last_point_coord, radius):
        """
        Rounds the corners between high and low impedance sections of the multimode coupler.

        Args:
            last_point_coord (tuple): The (x, y) coordinates at the end of the high impedance section.
            radius (float): The radius for the rounded corners.

        This method draws the rounded corners connecting the high and low impedance sections
        of the multimode coupler, ensuring a smooth transition between the two sections.
        It uses the `curve_corner` function to create the appropriate corner geometry.
        """
        c = self.structure
        shape = self.shape
        sgap = shape['s_lo']
        save1 = last_point_coord
        round_c = radius

        # Connection between high and low impedance sections
        self.curve_corner(c.s2, save1[0] - sgap / 2 - shape['w_hi'], save1[1], round_c, 2)
        self.curve_corner(c.s2, save1[0] + sgap / 2 + shape['w_hi'], save1[1], round_c, 3)

        # Outer corners of low impedance section
        self.curve_corner(c.s2, save1[0] - sgap / 2 - shape['w_lo'], save1[1], round_c, 0)
        self.curve_corner(c.s2, save1[0] + sgap / 2 + shape['w_lo'], save1[1], round_c, 1)

    def draw(self):
        """
        Constructs a multimode coupler structure with alternating high and low impedance sections,
        rounded corners, and a flux loop, starting from the specified position.
            
        Returns:
            None
        Notes:
            - The method modifies the internal structure (`self.structure`) by adding CPW straight sections,
              rounded corners, and a flux loop.
            - The coupler consists of a sequence of alternating low and high impedance sections, each
              connected by rounded corners, followed by a flux loop at the end.
            - The method restores the original reference point and direction at the end of execution.
        """

        c = self.structure
        shape = self.shape
        round_c = self.round_c


        # section 0 to the low impedence section
        target = 1090+120+850+50
        c.s2.last = (self.xpos, self.ypos)
        c.s2.last_direction = 90
        CPWStraight(c.s2, target, pinw=shape['s_lo'], gapw=shape['w_hi'])
        last_point_coord = c.s2.last
        

        #### ---------------------------     Filter
        self.round_corners_high_to_low(last_point_coord, round_c)
        # entering the filter
        for i in range(4):
            # section 1 low impedance
            last_point_coord = self.draw_low_impedance_pads(last_point_coord)

            # Round the corners between low and high impedance sections
            self.round_corners_low_to_high(last_point_coord, round_c)

            # section 2 high impedence
            last_point_coord = self.draw_high_impedance_section(last_point_coord)

            # Round the corners between high and low impedance sections
            self.round_corners_high_to_low(last_point_coord, round_c)
       
        # Final low impedance section 
        last_point_coord = self.draw_low_impedance_pads(last_point_coord)

        # filer to SQUID
        width = shape['w_hi']
        sgap = shape['s_hi']
        free_y = 100+890+813-1100+661-1901+3080-0.107-0.5*0+2.5-250+920+30.5-3125.46-0.5-120+153-640.5+150-1400+423.3-50
        loop_curve = 15
        CPWStraight(c.s2, shape['length_hi'] + free_y, pinw=sgap, gapw=width)
        central_point_coord = c.s2.last

        # Draw the flux loop
        self.draw_flux_loop(central_point_coord, loop_curve, square=square)

        # set the last point to initial reference
        c.s2.last = self.pin_ref
        c.s2.last_direction = self.pin_direction
        return None
# # -------------
#     def draw_launchers(structure, d, exclude=[]):
#         c = structure

def draw_bridge_junction(structure, d, junc_correction, flip, type=0, draw_layer1=True):
    c = structure
    ## alignment boxes for e-beam
    draw_launchers(c, d, exclude=[2, 3, 8])
    draw_launchers(c, d, exclude=[2, 3, 8])
    draw_square_alignment_marks(c, flip)
    ## Junction parameters

    xpos = 5650/2.-100
    ypos = -100.

    shape = {}

    if draw_layer1 == True:

        # 1 draw a full CPW line

        shape['launcher_width'] = 500
        shape['launcher_length'] = 500
        shape['launcher_taper'] = 500
        shape['ground'] = 2825
        shape['gap2ground'] = 200*1.0
        shape['cpw_ground_width'] = 200*1.25

        shape['cpw_pin'] = 12*1.25
        shape['cpw_gap'] = 6*1.25
        shape['cpw_ground4'] = shape['cpw_ground_width']
        shape['cpw_ground4_length'] = 200*1.25

        shape['cps_gap2'] = 12*1.25
        shape['cps_w'] = (shape['cpw_ground4_length']*2+shape['cpw_gap']*2+shape['cpw_pin']-shape['cps_gap2'])/2
        shape['cps_port2_length'] = 400-100

        shape['cps_pin'] = 200
        shape['cps_gap'] = 3
        shape['cps_gap_target'] = 3
        shape['cps_pin_target'] = 3
        shape['cps_pin_taper'] = 2000
        shape['scpw_length'] = 1000



        # full_cpw(c, xpos, ypos, shape)

        # 2 aburptly cpw to cps

        # cpw2cps_v0(c, xpos, ypos, shape)

        # 3 tapered cpw to cps

        balun = MultimodeBalun(structure = c, xpos=xpos, ypos=ypos, shape=shape)
        balun_end_coord = balun.draw()
        # reff = cpw2cps_v1(c, xpos, ypos, shape)

        shape1 = {}
        shape1['s_lo'] = 3
        shape1['s_hi'] = 3
        shape1['w_lo'] = 2422.5
        shape1['w_hi'] = 3

        shape1['length_lo'] = 5500 + 600
        shape1['length_hi'] = 5500 + 600

        shape1['left_padx'] = 1900
        shape1['left_pady'] = 700
        shape1['right_padx'] = 1900
        shape1['right_pady'] = 1600 + 400
        shape1['bar_thick'] = 50
        shape1['bar_length'] = 500
        shape1['squid_space'] = 30
        shape1['wirebond_x'] = 750 - 200
        shape1['wirebond_y'] = 400
        shape1['wirebond_dis'] = 600
        shape1['flux_loop_x'] = 190 - 60
        shape1['flux_loop_y'] = 75 - 30
        shape1['squid_width'] = 4
        shape1['squid_x'] = 39
        shape1['squid_y'] = 14 + shape1['squid_width'] * 2

        # shape['jj_gap'] = 13
        shape1['tripplr_jj'] = False

        xpos = balun_end_coord[0]
        ypos = balun_end_coord[1] 
        flux_line = MultimodeFluxLine(structure=c, xpos=xpos, ypos=ypos, shape=shape1)
        flux_line.draw()
        
        xpos = 5450 / 2.
        ypos = 46000 - 17 + 3 + 890
        coupler_pad = MultimodeCouplerPad(structure=c, xpos=xpos, ypos=ypos, shape=shape1)
        coupler_pad.draw()

class random():
    def __init__(self):
        pass

def output_wafer(MaskName = 'TestMask', author="Eesh Gupta", draw_layer1=1, draw_layer2=0, 
                 two_layer = False):
    draw_layer1 = 1
    # Meant for writing junctions on a wafer with metallic base layer

    # m = WaferMask(MaskName, diameter=76200., flat_angle=270., flat_distance=37100., wafer_padding=2500,
    #               chip_size=(5450, 48300),
    #               dicing_border=200, etchtype=False, wafer_edge=True,
    #               dashed_dicing_border=80, ndashes=6, dice_corner=True, square_arr=False)
    m = WaferMask(MaskName, diameter=101600., flat_angle=270., flat_distance=48200., wafer_padding=2500,
                  chip_size=(5450, 48300),
                  dicing_border=200, etchtype=False, wafer_edge=True,
                  dashed_dicing_border=80, ndashes=6, dice_corner=True, square_arr=False)

    # 4in wafer -> diameter=101600, flat_distance=48200
    # 3in wafer -> diameter=76200., flat_distance=37100.
    # 2in wafer -> diameter=50800., flat_distance=24100.

    points = [(-24000., -4000.), (-24000., 4000.), (24000., -4000.), (24000., 4000.)]
    points_medium = [(-24000. + 1200, -4000.), (-24000. + 1200, 4000.), (24000. - 1200, -4000.), (24000. - 1200, 4000.)]
    points_small = [(-24000. + 2000, -4000.), (-24000. + 2000, 4000.), (24000. - 2000, -4000.), (24000. - 2000, 4000.)]

    # Create the alignment crosses on the wafer. NOTE: currently only works for solid = True setting.
    # Something goes wrong with the layers cause that's built in otherwise. Hacky solution:
    solid = True  # Eesh 
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

    mm_chips_base = MMChipsBase()
    d = mm_chips_base.set_mask_init()

    # --------------------------#

    # chip type 0,1,2,3,4

    chip_configs = [
        # (junc_correction, flip, CHIPNAME, type)
        ({'width': 0.4, 'gap': 0.2}, False, 'C1A', 1),
        ({'width': 0.4, 'gap': 0.2}, False, 'C1B', 1),
        ({'width': 0.45, 'gap': 0.2}, False, 'C1C', 1),
        ({'width': 0.5, 'gap': 0.2}, False, 'C1D', 1),
        ({'width': 0.55, 'gap': 0.2}, False, 'C1E', 1),
        ({'width': 0.6, 'gap': 0.2}, False, 'C1F', 1),
        ({'width': 0.13, 'gap': 0.2}, False, 'C1G', 0),
        ({'width': 0.14, 'gap': 0.2}, False, 'C1H', 0),
        ({'width': 0.12, 'gap': 0.2}, False, 'C1I', 0),
        ({'width': 0.13, 'gap': 0.2}, False, 'C1J', 1),
        ({'width': 0.14, 'gap': 0.2}, False, 'C1K', 1),
        ({'width': 0.35, 'gap': 0.2}, False, 'C1L', 1),
        ({'width': 0.25, 'gap': 0.2}, False, 'C1M', 1),
        ({'width': 0.2, 'gap': 0.2}, False, 'C1N', 0),
    ]

    chips = []
    for idx, (junc_correction, flip, CHIPNAME, type_) in enumerate(chip_configs):
        chip_id_loc = (100, 100) if idx < 2 else (200, 100)
        c = Chip(CHIPNAME, author=author, size=m.chip_size, mask_id_loc=(5800, 6430),
                 chip_id_loc=chip_id_loc, textsize=(70, 70), two_layer=two_layer, solid=solid)
        draw_bridge_junction(c, d, junc_correction, flip, type=type_)
        chips.append(c)
    

    if draw_layer1 == 0 and draw_layer2 == 1:
        label = False
    else:
        label = True

    # m.add_chip(ct, 5, label=label)
    for idx, c in enumerate(chips[:13]):
        m.add_chip(c, 1, label=label)

    return m


if __name__ == "__main__":

    m = output_wafer()
    m.save()

    print("\n\n Chip names are:")
    print("_____________________")
    # for name in chip_names:
    #     print(name)
    print("_____________________\n\n")

    sleep(.1)
    show_structure = True
    open_klayout = True
    show_wafer = True
    open_dwgviewer = False

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