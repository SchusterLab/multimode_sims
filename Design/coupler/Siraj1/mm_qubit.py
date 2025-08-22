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

class MultimodeQubit(MMChipsBase):
    def __init__(self, structure, xpos, ypos, optical_configs= None, junction_configs=None, flag = 1):
        super().__init__()
        self.c = structure
        self.xpos = xpos
        self.ypos = ypos
        self.optical_configs = optical_configs
        self.junction_configs = junction_configs
        self.flag = flag

        self.pin_ref = self.c.s2.last
        self.pin_direction = self.c.s2.last_direction

        # set starting position
        self.c.s2.last = (xpos, ypos)
        self.c.s2.last_direction = 90 + 180 * flag
        self.coords={}

    # Optical Functions ---------------------------------------
    def draw_optical(self):
        self.shape = self.optical_configs

        self._draw_left_leg()
        self._draw_t_block_pair(flip_dir=True)
        self.c.s2.move(self.shape['jjgap'])
        self._draw_t_block_pair()
        self._draw_right_leg()
        self._draw_corner_rounds()
        # self._draw_top_marker()

        # restore original pin reference
        self.c.s2.last = self.pin_ref
        self.c.s2.last_direction = self.pin_direction

    def _draw_left_leg(self):
        CPWStraight(self.c.s2, self.shape['y_left_thick'], pinw=0, gapw=self.shape['x_left_thick'] / 2)
        CPWStraight(self.c.s2, self.shape['y_left_thin'], pinw=0, gapw=self.shape['x_thin'] / 2)

    def _draw_right_leg(self):
        CPWStraight(self.c.s2, self.shape['y_right_thin'], pinw=0, gapw=self.shape['x_thin'] / 2)
        CPWStraight(self.c.s2, self.shape['y_right_thick'], pinw=0, gapw=self.shape['x_right_thick'] / 2)

    def _draw_t_block_pair(self, flip_dir = False):
        x_last, y_last = self.c.s2.last
        if flip_dir:
            self.c.s2.last_direction = 180 + self.c.s2.last_direction
        CPWStraight(self.c.s2, 5, pinw=0, gapw=2.5)
        CPWStraight(self.c.s2, 5, pinw=0, gapw=7.5)
        if flip_dir:
            self.c.s2.last_direction = 180 + self.c.s2.last_direction
        self.c.s2.last = (x_last, y_last)

    def _draw_corner_rounds(self):
        rL = 40
        rS = 10
        x0 = self.xpos
        y0 = self.ypos
        x1 = self.c.s2.last[0]
        y1 = self.c.s2.last[1]
        s = self.shape

        cc = self.curve_corner  # shorthand

        if self.flag == 0:
            # left thick corners
            cc(self.c.s2, x0 + s['x_left_thick'] / 2, y0, rL, 1)
            cc(self.c.s2, x0 - s['x_left_thick'] / 2, y0, rL, 0)
            cc(self.c.s2, x0 - s['x_left_thick'] / 2, y0 + s['y_left_thick'], rL, 3)
            cc(self.c.s2, x0 + s['x_left_thick'] / 2, y0 + s['y_left_thick'], rL, 2)

            # left thin
            cc(self.c.s2, x0 - s['x_thin'] / 2, y0 + s['y_left_thick'], rS, 1)
            cc(self.c.s2, x0 + s['x_thin'] / 2, y0 + s['y_left_thick'], rS, 0)

            # right thick
            cc(self.c.s2, x1 + s['x_left_thick'] / 2, y1, rL, 2)
            cc(self.c.s2, x1 - s['x_left_thick'] / 2, y1, rL, 3)
            cc(self.c.s2, x1 - s['x_left_thick'] / 2, y1 - s['y_right_thick'], rL, 0)
            cc(self.c.s2, x1 + s['x_left_thick'] / 2, y1 - s['y_right_thick'], rL, 1)

            # right thin
            cc(self.c.s2, x1 - s['x_thin'] / 2, y1 - s['y_right_thick'], rS, 2)
            cc(self.c.s2, x1 + s['x_thin'] / 2, y1 - s['y_right_thick'], rS, 3)
        else:
            # mirrored
            cc(self.c.s2, x0 + s['x_left_thick'] / 2, y0, rL, 2)
            cc(self.c.s2, x0 - s['x_left_thick'] / 2, y0, rL, 3)
            cc(self.c.s2, x0 - s['x_left_thick'] / 2, y0 - s['y_left_thick'], rL, 0)
            cc(self.c.s2, x0 + s['x_left_thick'] / 2, y0 - s['y_left_thick'], rL, 1)

            cc(self.c.s2, x0 - s['x_thin'] / 2, y0 - s['y_left_thick'], rS, 2)
            cc(self.c.s2, x0 + s['x_thin'] / 2, y0 - s['y_left_thick'], rS, 3)

            cc(self.c.s2, x1 + s['x_left_thick'] / 2, y1, rL, 1)
            cc(self.c.s2, x1 - s['x_left_thick'] / 2, y1, rL, 0)
            cc(self.c.s2, x1 - s['x_left_thick'] / 2, y1 + s['y_right_thick'], rL, 3)
            cc(self.c.s2, x1 + s['x_left_thick'] / 2, y1 + s['y_right_thick'], rL, 2)

            cc(self.c.s2, x1 - s['x_thin'] / 2, y1 + s['y_right_thick'], rS, 1)
            cc(self.c.s2, x1 + s['x_thin'] / 2, y1 + s['y_right_thick'], rS, 0)

    def _draw_top_marker(self):
        self.c.s2.move(9900)
        x0, y0 = self.c.s2.last

        def draw_marker_set(x, y):
            self.c.s2.last = (x, y)
            CPWStraight(self.c.s2, 100, pinw=0, gapw=200)
            CPWStraight(self.c.s2, 500, pinw=200, gapw=100)
            CPWStraight(self.c.s2, 100, pinw=0, gapw=200)

        draw_marker_set(x0, y0)

        if self.flag == 0:
            draw_marker_set(x0 + 600, y0)
            draw_marker_set(x0 - 600, y0)
            self.c.s2.last = (x0 - 450, y0 + 350)
            CPWStraight(self.c.s2, 350, pinw=0, gapw=50)
        else:
            draw_marker_set(x0 - 600, y0)
            draw_marker_set(x0 + 600, y0)
            self.c.s2.last = (x0 + 450, y0 - 350)
            CPWStraight(self.c.s2, 350, pinw=0, gapw=50)
    # -------------------------------------------------------------
    def draw_ebeam(self):
        structure = self.c
        xpos = self.xpos
        ypos = self.ypos + 5
        shape = self.junction_configs
        flag = self.flag

        c = structure
        pin_ref = c.s2.pin_layer.last
        pin_direction = c.s2.pin_layer.last_direction
        gap_ref = c.s2.gap_layer.last
        gap_direction = c.s2.gap_layer.last_direction

        gap_new = 0.7
        gap_thin = 0.1

        # Pin Layer 
        c.s2.pin_layer.last_direction = 90
        c.s2.pin_layer.last = (xpos, ypos)
        CPWStraight(c.s2.pin_layer, shape['bar_length'], pinw=0, gapw=shape['bar_width'] / 2)

        # Dolan bridge start 
        ## taper 1
        CPWLinearTaper(c.s2.pin_layer, shape['taper_length1'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['bar_width'] / 2,
                    stop_gapw=shape['thin_bar_width1'] / 2)
        
        ## thin bar 1
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length1'], pinw=0, gapw=shape['thin_bar_width1'] / 2)

        ## taper 
        CPWLinearTaper(c.s2.pin_layer, shape['taper_length2'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['thin_bar_width1'] / 2,
                    stop_gapw=shape['jj_width'] / 2)
        
        ## thin bar 
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length2'], pinw=0, gapw=shape['jj_width'] / 2)
        
        ##gap
        c.s2.pin_layer.move(shape['jj_gap'])

        ## thin bar 
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['top_bar_jj_width'] / 2)
        
        ## upper taper 
        CPWLinearTaper(c.s2.pin_layer, shape['top_bar_taper_length'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['top_bar_jj_width'] / 2,
                    stop_gapw=shape['bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, shape['bar_length'], pinw=0, gapw=shape['bar_width'] / 2)

        # ----------------Under cut 
        uc_width = shape['undercut_width']
        uc_width_pads = shape['undercut_width_pads'] 
        c.s2.gap_layer.last = (xpos, ypos - uc_width_pads)
        c.s2.gap_layer.last_direction = 90

        bar_taper_offset = 0.3 # offset to ensure the diagonal region has enough undercut
        

        ## bar
        CPWStraight(c.s2.gap_layer, uc_width_pads, pinw=0, gapw=uc_width_pads + shape['bar_width'] / 2)
        CPWStraight(c.s2.gap_layer, shape['bar_length'], pinw=shape['bar_width'], gapw=uc_width_pads)

        ## taper 
        CPWLinearTaper(c.s2.gap_layer, shape['taper_length1'] , start_pinw=shape['bar_width'], stop_pinw=shape['thin_bar_width1'],
                    start_gapw= uc_width_pads,
                    stop_gapw= uc_width_pads)
        CPWStraight(c.s2.gap_layer, 1, pinw=shape['thin_bar_width1'], gapw=uc_width_pads) # so that the undercut is wide enough at taper to straight section transition

        ## thin bar 
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length1']-1, pinw=shape['thin_bar_width1'], gapw=uc_width)

        ## taper 2
        CPWLinearTaper(c.s2.gap_layer, shape['taper_length2'] , start_pinw=shape['thin_bar_width1'], stop_pinw=shape['jj_width'],
                    start_gapw= uc_width,
                    stop_gapw= uc_width)
        # CPWStraight(c.s2.gap_layer, 1, pinw=shape['jj_width'], gapw=uc_width_pads) # so that the undercut is wide enough at taper to straight section transition
        
        ## thin bar2 
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length2'], pinw=shape['jj_width'], gapw=uc_width)
        

        ##gap
        CPWStraight(c.s2.gap_layer, shape['jj_gap'], pinw=0, gapw=uc_width + shape['top_bar_jj_width']/2)
        
        ## thin bar 
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'] -1, pinw=shape['top_bar_jj_width'] , gapw=uc_width)
        
        ## upper taper 
        CPWStraight(c.s2.gap_layer, 1, pinw=shape['top_bar_jj_width'], gapw=uc_width_pads) # so that the undercut is wide enough at taper to straight section transition
        CPWLinearTaper(c.s2.gap_layer, shape['top_bar_taper_length'] , 
                    start_pinw=shape['top_bar_jj_width'] ,
                    stop_pinw=shape['bar_width'] ,
                    start_gapw=uc_width_pads,
                    stop_gapw=uc_width_pads)
        CPWStraight(c.s2.gap_layer, shape['bar_length'] , pinw=shape['bar_width'], gapw=uc_width_pads)
        CPWStraight(c.s2.gap_layer, uc_width_pads, pinw=0, gapw=uc_width_pads + shape['bar_width'] / 2)



        c.s2.pin_layer.last = pin_ref
        c.s2.gap_layer.last = gap_ref
        c.s2.pin_layer.last_direction = pin_direction
        c.s2.gap_layer.last_direction = gap_direction
    
#"""For Boomrah wafer had to reorient the junctions so align with coupler junctions"""

    def draw_ebeam_horizontal(self):
        structure = self.c
        xpos = self.xpos
        ypos = self.ypos + 12
        self.coords['origin/left_pad_start'] = (xpos -12, ypos)
        shape = self.junction_configs
        flag = self.flag

        c = structure
        pin_ref = c.s2.pin_layer.last
        pin_direction = c.s2.pin_layer.last_direction
        gap_ref = c.s2.gap_layer.last
        gap_direction = c.s2.gap_layer.last_direction

        gap_new = 0.7
        gap_thin = 0.1

        # Pin Layer 
        c.s2.pin_layer.last_direction = 0
        c.s2.pin_layer.last = self.coords['origin/left_pad_start']
        CPWStraight(c.s2.pin_layer, shape['bar_length'], pinw=0, gapw=shape['bar_width'] / 2)

        # Dolan bridge start 
        ## left taper 
        CPWLinearTaper(c.s2.pin_layer, shape['taper_length'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['bar_width'] / 2,
                    stop_gapw=shape['jj_width'] / 2)
        
        ## thin bar 
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['jj_width'] / 2)
        
        ##gap
        c.s2.pin_layer.move(shape['jj_gap'])
        
        ## right
        CPWLinearTaper(c.s2.pin_layer, shape['taper_length'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['top_bar_jj_width'] / 2,
                    stop_gapw=shape['bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, shape['bar_length']/2, pinw=0, gapw=shape['bar_width'] / 2)
        self.coords['end_of_right_pad'] = c.s2.pin_layer.last

        # connect to top bar using cpw bend 
        c.s2.pin_layer.last_direction = 0
        self.CPWBendNew(c.s2.pin_layer, angle=180, pinw=0, gapw=shape['bar_width'] / 2, radius=shape['bend_radius'],
                polyarc=1, segments=4,
                square=square)
        c.s2.pin_layer.last_direction = 180
        CPWStraight(c.s2.pin_layer, shape['horizontal_l_arm'], pinw=0, gapw=shape['bar_width'] / 2)


        

        # ----------------Under cut 
        uc_width = shape['undercut_width']
        xpos, ypos = self.coords['origin/left_pad_start']
        c.s2.gap_layer.last = (xpos - uc_width, ypos )
        c.s2.gap_layer.last_direction = 0
        

        ## bar
        CPWStraight(c.s2.gap_layer, uc_width, pinw=0, gapw=uc_width + shape['bar_width'] / 2)
        CPWStraight(c.s2.gap_layer, shape['bar_length'], pinw=shape['bar_width'], gapw=uc_width)

        ## taper 
        CPWLinearTaper(c.s2.gap_layer, shape['taper_length'], start_pinw=shape['bar_width'], stop_pinw=shape['jj_width'],
                    start_gapw= uc_width,
                    stop_gapw= uc_width)
        
        ## thin bar 
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['jj_width'], gapw=uc_width)
        
        ##gap
        CPWStraight(c.s2.gap_layer, shape['jj_gap'], pinw=0, gapw=uc_width + shape['top_bar_jj_width']/2)
        
        ## right taper 
        CPWLinearTaper(c.s2.gap_layer, shape['taper_length'], 
                    start_pinw=shape['top_bar_jj_width'] ,
                    stop_pinw=shape['bar_width'] ,
                    start_gapw=uc_width,
                    stop_gapw=uc_width)
        CPWStraight(c.s2.gap_layer, shape['bar_length']/2, pinw=shape['bar_width'], gapw=uc_width)

        # bend to top bar 
        c.s2.gap_layer.last_direction = 0
        self.CPWBendNew(c.s2.gap_layer, angle=180, pinw=shape['bar_width'], gapw = uc_width, radius=shape['bend_radius'],
                polyarc=1, segments=4,
                square=square)
        c.s2.pin_layer.last_direction = 180
        CPWStraight(c.s2.gap_layer, shape['horizontal_l_arm'], pinw=shape['bar_width'], gapw=uc_width)
        CPWStraight(c.s2.gap_layer, uc_width, pinw=0, gapw=shape['bar_width']/2+uc_width)

        # CPWStraight(c.s2.gap_layer, uc_width, pinw=0, gapw=uc_width + shape['bar_width'] / 2)

        ## connect to l arm 
        # xpos,ypos = self.coords['end_of_right_pad']
        # c.s2.gap_layer.last = (xpos, ypos - shape['bar_width']/2)



        c.s2.pin_layer.last = pin_ref
        c.s2.gap_layer.last = gap_ref
        c.s2.pin_layer.last_direction = pin_direction
        c.s2.gap_layer.last_direction = gap_direction

    def draw_ebeam_horizontal_small_loop(self):
        structure = self.c
        xpos = self.xpos - 10
        ypos = self.ypos + 6
        self.coords['origin/left_pad_start'] = (xpos, ypos)
        shape = self.junction_configs
        flag = self.flag

        c = structure
        pin_ref = c.s2.pin_layer.last
        pin_direction = c.s2.pin_layer.last_direction
        gap_ref = c.s2.gap_layer.last
        gap_direction = c.s2.gap_layer.last_direction

        gap_new = 0.7
        gap_thin = 0.1

        # Pin Layer 
        c.s2.pin_layer.last_direction = 90
        c.s2.pin_layer.last = self.coords['origin/left_pad_start']
        CPWStraight(c.s2.pin_layer, shape['bar_length'], pinw=0, gapw=shape['bar_width'] / 2)

        # Dolan bridge start 
        ## left taper 
        CPWLinearTaper(c.s2.pin_layer, shape['taper_length'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['bar_width'] / 2,
                    stop_gapw=shape['top_bar_jj_width'] / 2)

        ## slightly thinner bar  
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['top_bar_jj_width'] / 2)

        ## thiner bar to gap connection
        c.s2.pin_layer.last_direction = 90
        
        self.CPWBendNew(c.s2.pin_layer, angle=-90, pinw=0, gapw=shape['top_bar_jj_width'] / 2, radius=shape['bend_radius'],
                polyarc=1, segments=4,
                square=square)
        c.s2.pin_layer.last_direction = 0
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['top_bar_jj_width'] / 2)
        CPWLinearTaper(c.s2.pin_layer, shape['thin_bar_length'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['top_bar_jj_width'] / 2,
                    stop_gapw=shape['jj_width'] / 2)
        
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['jj_width'] / 2)
        
        ##gap
        c.s2.pin_layer.move(shape['jj_gap'])

        ## gap to taper connection
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['top_bar_jj_width'] / 2)
        self.CPWBendNew(c.s2.pin_layer, angle=90, pinw=0, gapw=shape['top_bar_jj_width'] / 2, radius=shape['bend_radius'],
                polyarc=1, segments=4,
                square=square)
        c.s2.pin_layer.last_direction = 90
        CPWStraight(c.s2.pin_layer, shape['thin_bar_length'], pinw=0, gapw=shape['top_bar_jj_width'] / 2)
        
        ## right
        CPWLinearTaper(c.s2.pin_layer, shape['taper_length'], start_pinw=0, stop_pinw=0,
                    start_gapw=shape['top_bar_jj_width'] / 2,
                    stop_gapw=shape['bar_width'] / 2)
        CPWStraight(c.s2.pin_layer, shape['bar_length'], pinw=0, gapw=shape['bar_width'] / 2)
        self.coords['end_of_right_pad'] = c.s2.pin_layer.last

        # connect to top bar using cpw bend 
        # c.s2.pin_layer.last_direction = 0
        # self.CPWBendNew(c.s2.pin_layer, angle=180, pinw=0, gapw=shape['bar_width'] / 2, radius=shape['bend_radius'],
        #         polyarc=1, segments=4,
        #         square=square)
        # c.s2.pin_layer.last_direction = 180
        # CPWStraight(c.s2.pin_layer, shape['horizontal_l_arm'], pinw=0, gapw=shape['bar_width'] / 2)


        

        # ----------------Under cut 
        uc_width = shape['undercut_width']
        xpos, ypos = self.coords['origin/left_pad_start']
        c.s2.gap_layer.last = (xpos, ypos -uc_width)
        c.s2.gap_layer.last_direction = 90
        

        ## bar
        CPWStraight(c.s2.gap_layer, uc_width, pinw=0, gapw=uc_width + shape['bar_width'] / 2)
        CPWStraight(c.s2.gap_layer, shape['bar_length'], pinw=shape['bar_width'], gapw=uc_width)

        CPWLinearTaper(c.s2.gap_layer, shape['taper_length'],
                    start_pinw=shape['bar_width'],
                    stop_pinw=shape['top_bar_jj_width'], 
                    start_gapw=uc_width,
                    stop_gapw=uc_width)

        ## slightly thinner bar  
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['top_bar_jj_width'], gapw=uc_width)

        ## thiner bar to gap connection
        c.s2.gap_layer.last_direction = 90
        
        self.CPWBendNew(c.s2.gap_layer, angle=-90, pinw=shape['top_bar_jj_width'] , gapw=uc_width, radius=shape['bend_radius'],
                polyarc=1, segments=4,
                square=square)
        
        c.s2.gap_layer.last_direction = 0
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['top_bar_jj_width'], gapw=uc_width)
        CPWLinearTaper(c.s2.gap_layer, shape['thin_bar_length'],
                    start_pinw=shape['top_bar_jj_width'] ,
                    stop_pinw=shape['jj_width'] , 
                    start_gapw=uc_width,
                    stop_gapw=uc_width)
        



        # ## taper 
        # CPWLinearTaper(c.s2.gap_layer, shape['taper_length'], start_pinw=shape['bar_width'], stop_pinw=shape['jj_width'],
        #             start_gapw= uc_width,
        #             stop_gapw= uc_width)
        
        # ## thin bar 
        # CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['jj_width'], gapw=uc_width)
        
        # ## thin bar to gap connection
        # c.s2.gap_layer.last_direction = 90
        # self.CPWBendNew(c.s2.gap_layer, angle=-90,  pinw=shape['jj_width'], gapw=uc_width, radius=shape['bend_radius'],
        #         polyarc=1, segments=4,
        #         square=square)
        c.s2.gap_layer.last_direction = 0
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['jj_width'], gapw=uc_width)


        ##gap
        CPWStraight(c.s2.gap_layer, shape['jj_gap'], pinw=0, gapw=uc_width + shape['top_bar_jj_width']/2)
        
        ## gap to taper connection
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['top_bar_jj_width'], gapw=uc_width)
        self.CPWBendNew(c.s2.gap_layer, angle=90, pinw=shape['top_bar_jj_width'], gapw=uc_width, radius=shape['bend_radius'],
                polyarc=1, segments=4,
                square=square)
        c.s2.gap_layer.last_direction = 90
        CPWStraight(c.s2.gap_layer, shape['thin_bar_length'], pinw=shape['top_bar_jj_width'], gapw=uc_width)


        ## right taper 
        CPWLinearTaper(c.s2.gap_layer, shape['taper_length'], 
                    start_pinw=shape['top_bar_jj_width'] ,
                    stop_pinw=shape['bar_width'] ,
                    start_gapw=uc_width,
                    stop_gapw=uc_width)
        CPWStraight(c.s2.gap_layer, shape['bar_length'], pinw=shape['bar_width'], gapw=uc_width)
        CPWStraight(c.s2.gap_layer, uc_width, pinw=0, gapw=uc_width + shape['bar_width'] / 2)



        c.s2.pin_layer.last = pin_ref
        c.s2.gap_layer.last = gap_ref
        c.s2.pin_layer.last_direction = pin_direction
        c.s2.gap_layer.last_direction = gap_direction

class MultimodeQubitChip(MMChipsBase):
    def __init__(self, structure,chip_defaults, flip=False):
        super().__init__()
        self.structure = structure
        self.d = chip_defaults
        self.flip = flip
        self.chipInit(self.structure, defaults=self.d)

    # def draw_alignment_marks(self):
    #     c = self.structure
    #     ## alignment boxes for e-beam
    #     self.draw_launchers(c, self.d, exclude=[2, 3, 8])
    #     self.draw_launchers(c, self.d, exclude=[2, 3, 8])
    #     self.draw_square_alignment_marks(c, self.flip)
    def draw_alignment_marks(self, draw_marks = False):
        c = self.structure
        ## alignment boxes for e-beam
        self.draw_launchers(c, self.d, exclude=[2, 3, 8])
        self.draw_launchers(c, self.d, exclude=[2, 3, 8])

        # draw the chip to cavity connection marker 
        c.s2.last = (1000, 9800 - 50)
        CPWStraight(c.s2, 100, pinw=0, gapw=1200 / 2)
        # c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)

        if draw_marks:
            self.draw_square_alignment_marks(c, self.flip)

       
    def draw_square_alignment_marks(self,structure, flip):
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
            # c.s2.last = (1000, 9800 - 50)
            # CPWStraight(c.s2, 100, pinw=0, gapw=1500 / 2)
            # c.s2.last = (c.size[0] / 2 - 1900, c.size[1] - 600)


# break down this function into 2 subfunctions: one for optical and one ebeam . 
# shape = qubit_optical_configs  should be renamed. So too for shape1 = test_optical_configs, junction_configs
# make sure in doc strings to include all parameters. 
#integrate these functions in multimode qubit chip class
    def draw_multimode_qubit_optical(self, qubit_optical_configs, flag):
        """
        Draws the optical (layer 1) features for the multimode qubit.

        Args:
            structure: The chip structure object.
            qubit_optical_configs (dict): Dictionary containing geometric parameters for the multimode qubit optical layer.
                Required keys include:
                    - 'x_thin', 'y_right_thin', 'x_right_thick', 'y_right_thick',
                    'y_left_thin', 'x_left_thick', 'y_left_thick', 'jjgap'
            xpos (float): X position for the qubit.
            ypos (float): Y position for the qubit.
            flag (int): Orientation flag (0 or 1).

        Returns:
            None
        """
        
        flag = 1
        if self.flip:
            flag = 0
        xpos = 1000.-200.
        ypos = 17700.-200.+7780
        if self.flip:
            ypos = 18100 - ypos
        
        qubit = MultimodeQubit(self.structure, xpos, ypos, qubit_optical_configs, flag= flag)
        qubit.draw_optical()
    
    def draw_multimode_qubit_optical(self, qubit_optical_configs, flag, xpos = 1000 - 200, ypos = 17700 - 200 + 7780):
        """
        Draws the optical (layer 1) features for the multimode qubit.

        Args:
            structure: The chip structure object.
            qubit_optical_configs (dict): Dictionary containing geometric parameters for the multimode qubit optical layer.
                Required keys include:
                    - 'x_thin', 'y_right_thin', 'x_right_thick', 'y_right_thick',
                    'y_left_thin', 'x_left_thick', 'y_left_thick', 'jjgap'
            xpos (float): X position for the qubit.
            ypos (float): Y position for the qubit.
            flag (int): Orientation flag (0 or 1).

        Returns:
            None
        """
        
        flag = 1
        if self.flip:
            flag = 0
        xpos = xpos
        ypos = ypos
        if self.flip:
            ypos = 18100 - ypos
        
        qubit = MultimodeQubit(self.structure, xpos, ypos, qubit_optical_configs, flag= flag)
        qubit.draw_optical()

    def draw_multimode_qubit_ebeam(self, junction_configs, junc_correction, flag = 1, 
                                   xpos = 1000 - 200, ypos = 17700 - 200 + 7780):
        """
        Draws the e-beam (layer 2) features for the multimode qubit, such as the junction frame.

        Args:
            junction_configs (dict): Dictionary containing geometric parameters for the Dolan bridge junction.
                Required keys include:
                    - 'shrink_length', 'overlap', 'overlap_width', 'offset',
                    'mjj_x_width', 'mjj_y_width', 'mjj_x_length', 'mjj_y_length',
                    'mjj_x_bar_length', 'mjj_y_bar_length', 'mjj_bar_width', 'gap'
            xpos (float): X position for the junction.
            ypos (float): Y position for the junction.
            flag (int): Orientation flag (0 or 1).

        Returns:
            None
        """
        # flag: direction of the junction
        xpos = xpos
        ypos = ypos -21 # 14200.-720+7700-21
        flag1 = 0
        if self.flip:
            flag1 = 1
            xpos = 1000. - 200.
            ypos = 18030-(14200. - 720)
        print('flag1', flag1    )
        junction_configs['jj_width'] = junc_correction['width']
        junction_configs['jj_gap'] = junc_correction['gap']
        qubit = MultimodeQubit(self.structure, xpos, ypos, junction_configs= junction_configs ,flag= flag)
        # qubit.draw_ebeam_horizontal_small_loop()
        qubit.draw_ebeam()

        
    # draw_junction_frame(self.structure, xpos, ypos, junction_configs, flag)
