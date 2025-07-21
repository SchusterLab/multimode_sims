# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 09:54:15 2021

@author: ziqian

Edited by: Eesh Gupta on 2025-07-16
"""

from mask_maker import * 
# from slab.circuits.mp_components import *

import numpy as np
square = 0

from mm_chips_base import MMChipsBase

class MultimodeCoupler(MMChipsBase): 

    def __init__(self, structure, xpos, ypos, loop_configs = None, junction_configs = None, optical_structure=0, jj = False):
        super().__init__()
        self.structure = structure
        # self.coupler = coupler
        self.optical_structure = optical_structure
        self.coords = {'origin': (xpos, ypos)}
        self.loop_configs = loop_configs
        self.junction_configs = junction_configs
        self.jj = jj



    def draw_loop_frame(self, optical_structure=0):
        c = self.structure
        xpos, ypos = self.coords['origin']
        loop = self.loop_configs

        pin_ref = c.s2.pin_layer.last
        pin_direction = c.s2.pin_layer.last_direction
        gap_ref = c.s2.gap_layer.last
        gap_direction = c.s2.gap_layer.last_direction

        # if optical_structure==1:
        #     c.s2.pin_layer.last = (xpos, ypos-26)
        #     c.s2.pin_layer.last_direction = 270
        #     CPWStraight(c.s2.pin_layer, 10, pinw=0, gapw=130 / 2)
        #     c.s2.gap_layer.last = (xpos, ypos - 26)
        #     c.s2.gap_layer.last_direction = 270
        #     CPWStraight(c.s2.gap_layer, 10, pinw=130, gapw=loop['gap_wide'])


        c.s2.pin_layer.last = (xpos, ypos)
        c.s2.pin_layer.last_direction = 270 # everything is drawn facing up  (not sure why 270 (though it would be 90))
        CPWStraight(c.s2.pin_layer, loop['width'], pinw=0, gapw=loop['x'] / 2) # top part of bar 
        c.s2.pin_layer.move(loop['y'])
        CPWStraight(c.s2.pin_layer, loop['width'], pinw=0, gapw=loop['x'] / 2) # bottom part of bar 

        c.s2.pin_layer.last = (xpos, ypos)
        c.s2.pin_layer.last_direction = 270
        CPWStraight(c.s2.pin_layer, 2*loop['width']+loop['y'], pinw=loop['x'], gapw=loop['overlap']) # left and right sides of bar 

        # draw gap layer
        c.s2.gap_layer.last = (xpos, ypos+loop['gap'])
        c.s2.gap_layer.last_direction = 270
        CPWStraight(c.s2.gap_layer, loop['gap'], pinw=0, gapw=loop['x'] / 2+loop['overlap']) #bottom bar outer: wider
        c.s2.gap_layer.move(loop['width'])
        CPWStraight(c.s2.gap_layer, loop['gap'], pinw=0, gapw=loop['x'] / 2) # bottom bar inner: thinner
        c.s2.gap_layer.move(loop['y']-loop['gap']*2)
        CPWStraight(c.s2.gap_layer, loop['gap'], pinw=0, gapw=loop['x'] / 2) # top bar inner: thinner
        c.s2.gap_layer.move(loop['width'])
        CPWStraight(c.s2.gap_layer, loop['gap'], pinw=0, gapw=loop['x'] / 2 + loop['overlap']) # top bar outer: wider

        c.s2.gap_layer.last = (xpos, ypos)
        c.s2.gap_layer.last_direction = 270
        CPWStraight(c.s2.gap_layer, 2 * loop['width'] + loop['y'], pinw=loop['x']+2*loop['overlap'], gapw=loop['gap_wide']) # side wider 
        c.s2.gap_layer.last = (xpos, ypos-loop['width']-loop['gap'])
        c.s2.gap_layer.last_direction = 270
        CPWStraight(c.s2.gap_layer, loop['y']-loop['gap']*2, pinw=loop['x'] - 2 * loop['gap_wide'],
                    gapw=loop['gap_wide'])

        ################################################
        c.s2.pin_layer.last = pin_ref
        c.s2.gap_layer.last = gap_ref
        c.s2.pin_layer.last_direction = pin_direction
        c.s2.gap_layer.last_direction = gap_direction


    def draw_djj_frame(self):
        '''
        Origin is at the center of top of bar 
        '''
        c = self.structure
        c = self.structure
        xpos, ypos = self.coords['origin']
        djj = self.junction_configs
        jj = self.jj


        pin_ref = c.s2.pin_layer.last
        pin_direction = c.s2.pin_layer.last_direction
        gap_ref = c.s2.pin_layer.last
        gap_direction = c.s2.pin_layer.last_direction
        ################################################

        # put junction gap
        gap_thin = 0.2
        
        # first cleanup the space

        c.s2.pin_layer.last_direction = 270
        c.s2.pin_layer.last = (xpos, ypos)
        total_junction_length = (djj['bar_length1'] + djj['bar_length2'] + djj['taper_length'] + djj['jj_length'] + djj['jj_gap'])
        ## draw the box that contains the junction
        CPWStraight(c.s2.pin_layer, djj['bar_width'], pinw=0,
                    gapw=total_junction_length/2)

        ## surrounding undercut (top and bottom) drawing down  (unnecessary?)
        c.s2.gap_layer.last_direction = 270
        c.s2.gap_layer.last = (xpos, ypos+gap_thin)

        CPWStraight(c.s2.gap_layer, gap_thin, pinw=0,
                    gapw=total_junction_length / 2) 
        c.s2.gap_layer.move(djj['bar_width'])
        CPWStraight(c.s2.gap_layer, gap_thin, pinw=0,
                    gapw=total_junction_length / 2)

        if jj:

            # Draw from left to right
            c.s2.pin_layer.last_direction = 0
            c.s2.pin_layer.last = (xpos-total_junction_length / 2, ypos-djj['bar_width']/2) # left edge 
            ## left bar 
            CPWStraight(c.s2.pin_layer, djj['bar_length1'], pinw=0,
                        gapw=djj['bar_width']/2)
            ## left taper
            CPWLinearTaper(c.s2.pin_layer, djj['taper_length'] , 0, 0, djj['bar_width'] / 2, djj['jj_width'] / 2)
            ## thin junction pin 
            CPWStraight(c.s2.pin_layer, djj['jj_length'], pinw=0,
                        gapw=djj['jj_width'] / 2)
            ## junction gap
            c.s2.pin_layer.move(djj['jj_gap'])
            ## right bar taper 
            CPWLinearTaper(c.s2.pin_layer, djj['bar_length2'] , 0, 0, djj['right_half_junction_width']/2, djj['bar_width'] / 2)

            # Draw undercut layer
            c.s2.gap_layer.last_direction = 0
            c.s2.gap_layer.last = (
            xpos - total_junction_length / 2, ypos - djj['bar_width'] / 2) # left edge
            ## left bar
            CPWStraight(c.s2.gap_layer, djj['bar_length1'], pinw=djj['bar_width'],
                        gapw=gap_thin)
            ## left taper
            CPWLinearTaper(c.s2.gap_layer, djj['taper_length'], djj['bar_width'], djj['jj_width'], gap_thin, gap_thin)
            ## thin junction pin
            CPWStraight(c.s2.gap_layer, djj['jj_length'], pinw=djj['jj_width'],
                        gapw=gap_thin)
            ## junction gap
            CPWStraight(c.s2.gap_layer, djj['jj_gap'], pinw=0,
                        gapw=djj['right_half_junction_width']/2+gap_thin)
            CPWLinearTaper(c.s2.gap_layer, djj['bar_length2'], djj['right_half_junction_width'], djj['bar_width'], gap_thin, gap_thin)

        ################################################
        c.s2.pin_layer.last = pin_ref
        c.s2.pin_layer.last_direction = pin_direction
        c.s2.pin_layer.last = gap_ref
        c.s2.pin_layer.last_direction = gap_direction

class MultimodeBalun(MMChipsBase):
    def __init__(self, structure, xpos, ypos, shape):
        super().__init__()
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

        # make the T's 
        c.s2.last_direction = 180
        c.s2.last = self.coords['left_bar_start']
        self._draw_t_block_pair(flip_dir=False)

        c.s2.last_direction = 0
        c.s2.last = self.coords['right_bar_start']
        self._draw_t_block_pair(flip_dir=False)
 
    def _draw_t_block_pair(self, flip_dir = False):
        self.c = self.structure
        x_last, y_last = self.c.s2.last
        if flip_dir:
            self.c.s2.last_direction = 180 + self.c.s2.last_direction
        CPWStraight(self.c.s2, 5, pinw=0, gapw=2.5)
        CPWStraight(self.c.s2, 5, pinw=0, gapw=7.5)
        if flip_dir:
            self.c.s2.last_direction = 180 + self.c.s2.last_direction
        self.c.s2.last = (x_last, y_last)

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
        self.CPWBendNew(c.s2, angle=-90, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi/2*radius # subtract the length of the bend

        # 3. First Straight section in the high impedance section
        sec2 = section_length
        CPWStraight(c.s2, sec2, pinw=sgap, gapw=width)
        remain = remain - sec2

        # 4. Second bend in the high impedance section
        self.CPWBendNew(c.s2, angle=180, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi * radius

        # 5. Second Straight section in the high impedance section (the largest)
        sec3 = section_length*2+radius*2 # larg length (2x low impedance pad length)
        CPWStraight(c.s2, sec3, pinw=sgap, gapw=width)
        remain = remain - sec3

        # 6. Third bend in the high impedance section
        self.CPWBendNew(c.s2, angle=-180, pinw=sgap, gapw=width, radius=radius,
                polyarc=1, segments=4,
                square=square)
        remain = remain - np.pi * radius

        # 7. Third Straight section in the high impedance section
        CPWStraight(c.s2, sec2, pinw=sgap, gapw=width)
        remain = remain - sec2
        self.CPWBendNew(c.s2, angle=90, pinw=sgap, gapw=width, radius=radius,
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
        self.CPWBendNew(c.s2, angle=90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_x'] / 2 - loop_curve * 2 - sgap / 2, pinw=0, gapw=width / 2)
        self.CPWBendNew(c.s2, angle=-90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_y'] - loop_curve * 2, pinw=0, gapw=width / 2)
        self.CPWBendNew(c.s2, angle=-90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_x'] - loop_curve * 2 + width, pinw=0, gapw=width / 2)

        # Save the starting point for the left half of the loop
        left_start = (centeral_coords[0] + sgap / 2 + width / 2, centeral_coords[1])
        c.s2.last_direction = 90
        c.s2.last = left_start

        # Left half of the flux loop
        self.CPWBendNew(c.s2, angle=-90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_x'] / 2 - loop_curve * 2 - sgap / 2, pinw=0, gapw=width / 2)
        self.CPWBendNew(c.s2, angle=90, pinw=0, gapw=width / 2, radius=loop_curve,
                   polyarc=1, segments=4, square=square)
        CPWStraight(c.s2, shape['flux_loop_y'] - loop_curve * 2, pinw=0, gapw=width / 2)
        self.CPWBendNew(c.s2, angle=90, pinw=0, gapw=width / 2, radius=loop_curve,
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

# Turn the function below into a method of a new class called MultimodeCouplerChip  
# The class should take parameters such as structure, d = defaults, flip
# There should be 3 methods: one that draws the alignment marks for ebeam, one that draws the coupler optical (layer 1), 
# and one that draws the coupler ebeam (layer 2)
class MultimodeCouplerChip(MMChipsBase):
    def __init__(self, structure,chip_defaults, flip=False):
        super().__init__()
        self.structure = structure
        self.d = chip_defaults
        self.flip = flip

    def draw_alignment_marks(self):
        c = self.structure
        ## alignment boxes for e-beam
        self.draw_launchers(c, self.d, exclude=[2, 3, 8])
        self.draw_launchers(c, self.d, exclude=[2, 3, 8])
        self.draw_square_alignment_marks(c, self.flip)

    def draw_coupler_optical(self, balun_configs, coupler_configs):
        """
        Draws the coupler optical (layer 1) for the multimode coupler chip.

        Args:
            balun_configs (dict): Dictionary containing all geometric parameters for the balun.
                Required keys include:
                    - 'launcher_width', 'launcher_length', 'launcher_taper', 'ground', 'gap2ground',
                      'cpw_ground_width', 'cpw_pin', 'cpw_gap', 'cpw_ground4', 'cpw_ground4_length',
                      'cps_gap2', 'cps_w', 'cps_port2_length', 'cps_pin', 'cps_gap', 'cps_gap_target',
                      'cps_pin_target', 'cps_pin_taper', 'scpw_length'
            coupler_configs (dict): Dictionary containing all geometric parameters for the coupler/flux line/pad.
                Required keys include:
                    - 's_lo', 's_hi', 'w_lo', 'w_hi', 'length_lo', 'length_hi', 'left_padx', 'left_pady',
                      'right_padx', 'right_pady', 'bar_thick', 'bar_length', 'squid_space', 'wirebond_x',
                      'wirebond_y', 'wirebond_dis', 'flux_loop_x', 'flux_loop_y', 'squid_width', 'squid_x',
                      'squid_y', 'tripplr_jj'

        Returns:
            None
        """
        c = self.structure

        # Unpack starting coordinates
        xpos = balun_configs.get('xpos', 5650/2.-100)
        ypos = balun_configs.get('ypos', -100.)

        # Draw balun
        balun = MultimodeBalun(structure=c, xpos=xpos, ypos=ypos, shape=balun_configs)
        balun_end_coord = balun.draw()

        # Draw flux line
        xpos = balun_end_coord[0]
        ypos = balun_end_coord[1]
        flux_line = MultimodeFluxLine(structure=c, xpos=xpos, ypos=ypos, shape=coupler_configs)
        flux_line.draw()

        # Draw coupler pad
        xpos = 5450 / 2.
        ypos = 46000 - 17 + 3 + 890
        coupler_pad = MultimodeCouplerPad(structure=c, xpos=xpos, ypos=ypos, shape=coupler_configs)
        coupler_pad.draw()


    def draw_coupler_ebeam(self, loop_configs, junction_configs, junc_correction):
        """
        Draws the coupler e-beam (layer 2) for the multimode coupler chip.

        Args:
            loop_configs (dict): Dictionary containing geometric parameters for the SQUID loop.
                Required keys include:
                    - 'x': float, horizontal dimension of the loop area
                    - 'y': float, vertical dimension of the loop area
                    - 'gap': float, gap width for the loop
                    - 'gap_wide': float, wide gap width for the loop
                    - 'width': float, width of the loop bar
                    - 'overlap': float, width of the side bars
            junction_configs (dict): Dictionary containing geometric parameters for the Dolan bridge junction.
                Required keys include:
                    - 'bar_width': float, width of the bar
                    - 'bar_length1': float, length of the left bar
                    - 'bar_length2': float, length of the right bar
                    - 'taper_length': float, length of the taper
                    - 'jj_length': float, length of the junction
                    - 'jj_width': float, width of the junction
                    - 'jj_gap': float, gap of the junction
                    - 'right_half_junction_width': float, width of the right half of the junction
            junc_correction (dict): Dictionary containing correction parameters for the junction.
                Required keys include:
                    - 'width': float, width correction for the junction
                    - 'gap': float, gap correction for the junction

        Returns:
            None
        """
        c = self.structure

        # Draw SQUID loop
        xpos = 5450 / 2.
        ypos = 46000 - 17 + 3 +64 #+ 890 + 824
        coupler = MultimodeCoupler(structure=c, xpos=xpos, ypos=ypos, loop_configs=loop_configs, optical_structure=0)
        coupler.draw_loop_frame(optical_structure=0)

        # Draw Dolan bridge junction (top)
        djj = dict(junction_configs)  # make a copy to avoid mutation
        djj['jj_width'] = junc_correction['width']
        djj['jj_gap'] = junc_correction['gap']
        coupler = MultimodeCoupler(structure=c, xpos=xpos, ypos=ypos, junction_configs=djj, jj=True)
        coupler.draw_djj_frame()
        shifted_ypos = ypos - loop_configs['y'] - loop_configs['width']
        coupler.coords['origin'] = (xpos, shifted_ypos)
        coupler.draw_djj_frame()

       
