from phidl import Device, geometry as pg, quickplot as qp, Layer, LayerSet

class MattResonator(Device):
    """ Draws one of Matt's resonators """

    def __init__(self, substrate_layer=0, device_layer=1, chip_x = 7, chip_y = 49, params = None, type_num = 0, dicing = True, undercut = False):

        super().__init__(name = 'tripole_stripline')

        self.substrate_layer = substrate_layer
        self.device_layer = device_layer

        if params == None:
            self.params =  {'cpl': 0.5e3,
                            'L': 14e3,
                            'w1': 10,
                            'w2': 400,
                            'd1': 10,
                            'd2': 1200} # in um
        else:
            self.params = params
        dicing_border = self._draw_chip_dicing_border()
        self._draw_chips()
        self._draw_metal(dicing_border)
        
        #self._add_dicing_marks()
        if undercut == True:
            self._draw_device_undercut()
        self._add_chip_info(dicing_border, type_num)
        #if dicing == True:
        #    self._add_dicing_lines()
        self.move((-self.x, -self.y))

  
    def _add_dicing_marks(self):
        dicing_lane_width = 450

        xmax = self.xmax
        ymax = self.ymax
        xmin = self.xmin
        ymin = self.ymin

        self << pg.cross(length=250, width=100, layer=3).move((xmax + dicing_lane_width/2, ymax + dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmax + dicing_lane_width/2, ymin - dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmin - dicing_lane_width/2, ymax + dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmin - dicing_lane_width/2, ymin - dicing_lane_width/2)) # used to be 100 20
    
    def _draw_chips(self):
        TSL = Device('tripole_stripline')
        ls = LayerSet()
        ls.add_layer(name = 'substrate', gds_layer = 8, gds_datatype = 8,
                    description = 'substrate', color = 'green')
        # create substrate
        substrate = Device('substrate')
        substrate_x = (0, 4e3, 4e3, 0)
        substrate_y = (0, 0, 45e3, 45e3)
        substrate_obj = substrate.add_polygon([substrate_x, substrate_y], layer = ls['substrate'])
        #substrate_obj.move([0, 0])
        #pattern_inv = pg.boolean(substrate, TSL, operation='A-B', layer = ls['metal'])

        self << substrate

    def _draw_chip_dicing_border(self):
        TSL = Device('tripole_stripline')
        ls = LayerSet()
        ls.add_layer(name = 'substrate', gds_layer = 0, gds_datatype = 0,
                    description = 'substrate', color = 'green')
        # create substrate
        substrate = Device('substrate')
        substrate_x = (0, 4e3+100, 4e3+100, 0)
        substrate_y = (0, 0, 45e3+100, 45e3+100)
        substrate_obj = substrate.add_polygon([substrate_x, substrate_y], layer = ls['substrate'])
        substrate_obj.move([-50, -50])
        #pattern_inv = pg.boolean(substrate, TSL, operation='A-B', layer = ls['metal'])

        self << substrate
        return substrate_obj

    def _draw_metal(self, dicing_border):  
        TSL = Device('tripole_stripline')

        ls = LayerSet() # Create a blank LayerSet
        
        ls.add_layer(name = 'metal', gds_layer = 1, gds_datatype = 0,
                    description = 'metal', color = 'red')
        # create substrate


        strip1_x = (0, self.params['w1'], self.params['w1'], 0)


        strip1_y = (0, 0, self.params['cpl']+self.params['L'], self.params['cpl']+self.params['L'])
        strip1 = TSL.add_polygon([strip1_x, strip1_y], layer = ls['metal'])

        # creating strip2 (middle)
        strip2_x = (0, self.params['w2'], self.params['w2'], 0)
        strip2_y = (0, 0, self.params['L'], self.params['L'])
        strip2 = TSL.add_polygon([strip2_x, strip2_y], layer = ls['metal'])
        strip2.move([self.params['w1']+self.params['d1'], self.params['cpl']])

        # creating strip3 (outer)
        strip2_x = (0, self.params['w2'], self.params['w2'], 0)
        strip2_y = (0, 0, self.params['L'], self.params['L'])
        strip2 = TSL.add_polygon([strip2_x, strip2_y], layer = ls['metal'])
        strip2.move([self.params['w1']+self.params['d1']+self.params['w2']+self.params['d2'], self.params['cpl']])
        
        # invert the pattern
        TSL.move([(-.72+2)*1e3, (15.5+6)*1e3]) # add 1.5 as their is a 1.5mm groove in the plate on the back of cavity, add .5 as the antenna on simulation is .5 off
        pattern_inv = pg.boolean(dicing_border, TSL, operation='A-B', layer = ls['metal'])
        
        
        
        self << pattern_inv
        

    def _draw_device_undercut(self):
        TSL = Device('tripole_stripline')

        ls = LayerSet() # Create a blank LayerSet
        ls.add_layer(name = 'metal', gds_layer = 1, gds_datatype = 0,
                    description = 'metal', color = 'red')
        # create substrate


        strip1_x = (0, self.params['w1'], self.params['w1'], 0)


        strip1_y = (0, 0, self.params['cpl']+self.params['L'], self.params['cpl']+self.params['L'])
        strip1 = TSL.add_polygon([strip1_x, strip1_y], layer = ls['metal'])

        # creating strip2 (middle)
        strip2_x = (0, self.params['w2'], self.params['w2'], 0)
        strip2_y = (0, 0, self.params['L'], self.params['L'])
        strip2 = TSL.add_polygon([strip2_x, strip2_y], layer = ls['metal'])
        strip2.move([self.params['w1']+self.params['d1'], self.params['cpl']])

        # creating strip3 (outer)
        strip2_x = (0, self.params['w2'], self.params['w2'], 0)
        strip2_y = (0, 0, self.params['L'], self.params['L'])
        strip2 = TSL.add_polygon([strip2_x, strip2_y], layer = ls['metal'])
        strip2.move([self.params['w1']+self.params['d1']+self.params['w2']+self.params['d2'], self.params['cpl']])
        TSL.move([(-.72+2)*1e3, (15.5+6)*1e3]) # add 1.5 as their is a 1.5mm groove in the plate on the back of cavity, add .5 as the antenna on simulation is .5 off

        TSL_undercut = Device('tripole_stripline_undercut')

        ls.add_layer(name = 'substrate', gds_layer = 0, gds_datatype = 0,
                    description = 'substrate', color = 'green')
        ls.add_layer(name = 'undercut', gds_layer = 2, gds_datatype = 0,
                    description = 'undercut', color = 'red')
        # create substrate

        undercut = 1
        # strip1
        strip1_x = (-undercut, self.params['w1']+undercut, self.params['w1']+undercut, -undercut)
        strip1_y = (-undercut, -undercut, self.params['cpl']+self.params['L']+undercut, self.params['cpl']+self.params['L']+undercut)
        strip1_undercut = TSL_undercut.add_polygon([strip1_x, strip1_y], layer = ls['undercut'])

        # creating strip2 (middle)
        strip2_x = (-undercut, self.params['w2']+undercut, self.params['w2']+undercut, -undercut)
        strip2_y = (-undercut, -undercut, self.params['L']+undercut, self.params['L']+undercut)
        strip2_undercut = TSL_undercut.add_polygon([strip2_x, strip2_y], layer = ls['undercut'])
        strip2_undercut.move([self.params['w1']+self.params['d1'], self.params['cpl']])

        # creating strip3 (outer)
        strip2_x = (-undercut, self.params['w2']+undercut, self.params['w2']+undercut, -undercut)
        strip2_y = (-undercut, -undercut, self.params['L']+undercut, self.params['L']+undercut)
        strip2_undercut = TSL_undercut.add_polygon([strip2_x, strip2_y], layer = ls['undercut'])
        strip2_undercut.move([self.params['w1']+self.params['d1']+self.params['w2']+self.params['d2'], self.params['cpl']])
        TSL_undercut.move([(-.72+2)*1e3, (15.5+6)*1e3]) # add 1.5 as their is a 1.5mm groove in the plate on the back of cavity, add .5 as the antenna on simulation is .5 off

        print('hi')
        pattern_inv = pg.boolean(TSL_undercut, TSL, operation='A-B', layer = ls['undercut'])

        self << pattern_inv

    def _add_chip_info(self, type_num):
        # Create a new device for the text
        TSL = Device('tripole_stripline')
        ls = LayerSet() # Create a blank LayerSet
        ls.add_layer(name = 'text', gds_layer = 3, gds_datatype = 0,
                    description = 'text', color = 'green')
        # Create a layer for the text
        text_layer = 3  # Use a simple tuple for layer, datatype
        
        # Generate the text
        text1 = pg.text(text=f'v2 070825 MC', size=300, layer=ls['text'])
        text1.move([700,1000])
        
        # Add the text to the text_device
        text_ref1 = TSL.add_ref(text1)
        
        remainder = type_num%9
        dash_val = type_num//9
        text2 = pg.text(text=f'T{remainder+1}_{dash_val+1}', size=1000, layer=ls['text'])
        text2.move([300,2000])
        

        # Add the text to the text_device
        text_ref2 = TSL.add_ref(text2)

        # Add the text_device to the main device (self)
        self.add_ref(TSL)
        
        self << TSL
        # If you need to position the text at a specific location:
        # text_ref.move([x, y])  # Uncomment and replace x,y with coordinates
    

    def _add_dicing_lines(self):
        TSL = Device('tripole_stripline')
        ls = LayerSet()
        ls.add_layer(name = 'dicing', gds_layer = 4, gds_datatype = 0,
                    description = 'dicing', color = 'purple')
        # create substrate
        dicing = Device('dicing')
        dicing_x = (-0.05e3, 4.05e3, 4.05e3, -.05e3)
        dicing_y = (-0.05e3, -0.05e3, 45.05e3, 45.05e3)
        dicing_obj = dicing.add_polygon([dicing_x, dicing_y], layer = ls['dicing'])
        #substrate_obj.move([0, 0])
        #pattern_inv = pg.boolean(substrate, TSL, operation='A-B', layer = ls['metal'])
        #dicing.move([1, 0])
        self << dicing


class BlankSquares(Device):
    """ Draws a blank square for the resonator """

    def __init__(self, substrate_layer=0, device_layer=1, params = None, metal = False, size = 7):

        super().__init__(name = 'blank_square')

        self.substrate_layer = substrate_layer
        self.device_layer = device_layer

        if params == None:
            self.params =  {'size': 10e3} # in um
        else:
            self.params = params

        self._draw_chips(metal, size)
        self._add_dicing_marks()

    def _draw_chips(self, metal, size):
        TSL = Device('tripole_stripline')
        ls = LayerSet()
        ls.add_layer(name = 'blank substrate', gds_layer = 5, gds_datatype = 0,
                    description = 'blank substrate', color = 'green')
        
        ls.add_layer(name = 'blank metal', gds_layer = 6, gds_datatype = 0,
                    description = 'blank metal', color = 'red')
        size = size * 1e3
        # create substrate]
        if metal == True:
            layer_type = 'blank substrate'
        else:
            layer_type  = 'blank metal'   
        layer = Device(layer_type)
        layer_x = (0, size, size, 0)
        layer_y = (0, 0, size, size)
        print(layer_x)
        layer_obj = layer.add_polygon([layer_x, layer_y], layer = ls[layer_type])
        self << layer
        #substrate_obj.move([0, 0])
        #pattern_inv = pg.boolean(substrate, TSL, operation='A-B', layer = ls['metal'])
        
    def _add_dicing_marks(self):
        dicing_lane_width = 0

        xmax = self.xmax
        ymax = self.ymax
        xmin = self.xmin
        ymin = self.ymin

        self << pg.cross(length=250, width=100, layer=3).move((xmax + dicing_lane_width/2, ymax + dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmax + dicing_lane_width/2, ymin - dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmin - dicing_lane_width/2, ymax + dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmin - dicing_lane_width/2, ymin - dicing_lane_width/2)) # used to be 100 20
    
    

class JustDicingMarks(Device):
    """ Draws one of Matt's resonators """

    def __init__(self, substrate_layer=0, device_layer=1, chip_x = 7, chip_y = 49, params = None, type_num = 0, dicing = True):

        super().__init__(name = 'tripole_stripline')

        self.substrate_layer = substrate_layer
        self.device_layer = device_layer

        if params == None:
            self.params =  {'cpl': 0.5e3,
                            'L': 14e3,
                            'w1': 10,
                            'w2': 400,
                            'd1': 10,
                            'd2': 1200} # in um
        else:
            self.params = params
        self._draw_chips()
        self._draw_metal()
        #self._add_dicing_marks()
        self._draw_device_undercut()
        self._add_chip_info(type_num)
        #if dicing == True:
        #    self._add_dicing_lines()
        self.move((-self.x, -self.y))
        
    def _add_dicing_marks(self):
        dicing_lane_width = 450

        xmax = self.xmax
        ymax = self.ymax
        xmin = self.xmin
        ymin = self.ymin

        self << pg.cross(length=250, width=100, layer=3).move((xmax + dicing_lane_width/2, ymax + dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmax + dicing_lane_width/2, ymin - dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmin - dicing_lane_width/2, ymax + dicing_lane_width/2))
        self << pg.cross(length=250, width=100, layer=3).move((xmin - dicing_lane_width/2, ymin - dicing_lane_width/2)) # used to be 100 20