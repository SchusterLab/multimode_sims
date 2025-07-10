# MaskMaker Module Summary

This document provides a comprehensive summary of all classes and functions in the MaskMaker module, which is used for designing and generating mask layouts for superconducting circuits.

## Table of Contents
- [Error Handling](#error-handling)
- [Point-wise Operations](#point-wise-operations)
- [Mask and Chip Generation](#mask-and-chip-generation)
- [CPW Components](#cpw-components)
- [Capacitors and Couplers](#capacitors-and-couplers)
- [Channel Components](#channel-components)
- [Alignment and Text Components](#alignment-and-text-components)
- [Utility Functions](#utility-functions)

## Error Handling

### MaskError
An exception class raised when invalid parameters are used in MaskMaker functions.
- **Parameters**: `value` - Error message string

## Point-wise Operations

### distance(tuple1, tuple2)
Calculates the Euclidean distance between two points.

### ang2pt(direction, distance)
Converts an angle and distance to a point offset (dx, dy).
- **Parameters**:
  - `direction` - Angle in degrees
  - `distance` - Distance to travel

### rotate_pt(p, angle, center=(0, 0))
Rotates a point around a center point by a specified angle.
- **Parameters**:
  - `p` - Point to rotate (x, y)
  - `angle` - Rotation angle in degrees (CCW)
  - `center` - Center of rotation (default: origin)

### rotate_pts(points, angle, center=(0, 0))
Rotates an array of points around a center point.

### translate_pt(p, offset)
Translates a point by an offset.
- **Parameters**:
  - `p` - Point to translate (x, y)
  - `offset` - Translation offset (x, y)

### translate_pts(points, offset)
Translates an array of points by an offset.

### orient_pt(p, angle, offset)
Rotates a point by an angle and then translates it.
- **Parameters**:
  - `p` - Point to orient (x, y)
  - `angle` - Rotation angle in degrees
  - `offset` - Translation offset (x, y)

### orient_pts(points, angle, offset)
Orients an array of points by rotating and translating them.

### scale_pt(p, scale)
Scales a point by a scale factor.
- **Parameters**:
  - `p` - Point to scale (x, y)
  - `scale` - Scale factors (sx, sy)

### scale_pts(points, scale)
Scales an array of points by a scale factor.

### mirror_pt(p, axis_angle, axis_pt)
Mirrors a point about a line at a specified angle passing through a specified point.
- **Parameters**:
  - `p` - Point to mirror (x, y)
  - `axis_angle` - Angle of the mirror axis in degrees
  - `axis_pt` - Point on the mirror axis

### mirror_pts(points, axis_angle, axis_pt)
Mirrors an array of points about a line.

## Mask and Chip Generation

### WaferMask
Class for placing chips on a wafer with a flat. Contains functions to layout chips, add chips to the mask, and create a manifest.
- **Key Parameters**:
  - `name` - Name of the mask
  - `diameter` - Wafer diameter
  - `flat_angle` - Angle of the wafer flat
  - `flat_distance` - Distance from center to flat
  - `chip_size` - Size of individual chips
  - `dicing_border` - Width of dicing border
  - `etchtype` - Type of etching (True for standard with dicing borders)

### Chip
Class which contains structures for a chip design.
- **Key Parameters**:
  - `name` - Name of the chip
  - `author` - Author of the chip design
  - `size` - Size of the chip (x, y)
  - `two_layer` - Layer index of second layer 
  - `solid`- Layer Index (0)
  - `layer` - Layer for the chip

### Structure
Keeps track of current location and direction for drawing structures.
- **Key Parameters**:
  - `chip` - Chip to add the structure to
  - `start` - Starting point
  - `direction` - Starting direction in degrees
  - `layer` - Layer for the structure
  - `defaults` - Default values for the structure

### ChipBorder
Creates a border around a chip for dicing.
- **Parameters**:
  - `chip` - Chip to add the border to
  - `border_thickness` - Thickness of the border

### DashedChipBorder
Creates a dashed border around a chip for e-beam drawing and dicing.
- **Key Parameters**:
  - `chip` - Chip to add the border to
  - `border_thickness` - Thickness of the border
  - `dash_width` - Width of dashes
  - `dash_length` - Length of dashes
  - `ndashes` - Number of dashes
  - `dice_corner` - Whether to add corner marks

## CPW Components

### Launcher
Creates a CPW launcher pad with a taper to the main transmission line.
- **Key Parameters**:
  - `structure` - Structure to add the launcher to
  - `flipped` - Whether the launcher is flipped
  - `pad_length` - Length of the pad
  - `taper_length` - Length of the taper
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap

### Box
Creates a one-layer box that can be launched asymmetrically.
- **Key Parameters**:
  - `structure` - Structure to add the box to
  - `length` - Length of the box
  - `width` - Width of the box
  - `offset` - Offset from the current position

### CoupledStraight
Creates a straight section of coupled CPW transmission line.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `length` - Length of the section
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap
  - `center_gapw` - Width of the center gap

### LabLogo
Creates the Schuster Lab logo.
- **Key Parameters**:
  - `structure` - Structure to add the logo to
  - `logo_width` - Width of the logo
  - `logo_starting_point` - Starting point for the logo

### ExtendedFluxLine
Creates a split flux bias design with the option to extend the flux line from the ground plane.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `flux_extention` - Extension of the flux line
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap
  - `flux_length` - Length of the flux line
  - `flux_seperation` - Separation of the flux line

### CPWStraight
Creates a straight section of CPW transmission line.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `length` - Length of the section
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap

### CPWQubitBox
Creates a straight section of CPW transmission line with fingers in the ground plane to add a capacitor.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `fingerlen` - Length of the fingers
  - `fingerw` - Width of the fingers
  - `finger_gapw` - Gap width between fingers
  - `finger_no` - Number of fingers
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap

### CPWLinearTaper
Creates a section of CPW which linearly tapers from one set of dimensions to another.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `length` - Length of the taper
  - `start_pinw` - Starting pin width
  - `stop_pinw` - Ending pin width
  - `start_gapw` - Starting gap width
  - `stop_gapw` - Ending gap width

### CPWTaper
Simplified version of CPWLinearTaper that scales by a ratio.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `length` - Length of the taper
  - `ratio` - Scaling ratio
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap

### CPWBend
Creates a CPW bend.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `turn_angle` - Angle of the turn in degrees
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap
  - `radius` - Radius of the bend

### CPWL
Combination of CPWStraight, CPWBend, and another CPWStraight.
- **Key Parameters**:
  - `s` - Structure to add the component to
  - `d1` - Length of the first straight segment
  - `t1` - Turn angle of the bend
  - `r1` - Bend radius
  - `d2` - Length of the second straight segment

### CPWSturn
Combination of CPWStraight, CPWBend, CPWStraight, CPWBend, and another CPWStraight.
- **Key Parameters**:
  - `s` - Structure to add the component to
  - `d1`, `d2`, `d3` - Lengths of straight segments
  - `t1`, `t2` - Turn angles
  - `r1`, `r2` - Bend radii

### CPWWiggles
Creates CPW wiggles (meanders).
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `num_wiggles` - Number of wiggles
  - `total_length` - Total length of the meander
  - `start_up` - Whether to start with a CCW turn
  - `radius` - Radius of the bends

### CoupledBend
Creates a coupled CPW bend.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `turn_angle` - Angle of the turn in degrees
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap
  - `center_gapw` - Width of the center gap
  - `radius` - Radius of the bend

### CoupledWiggles
Creates coupled CPW wiggles (meanders).
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `num_wiggles` - Number of wiggles
  - `total_length` - Total length of the meander
  - `start_up` - Whether to start with a CCW turn
  - `radius` - Radius of the bends
  - `center_gapw` - Width of the center gap

### CPWWigglesByLength
Updated version of CPWWiggles which is more general, allowing for different starting angles and symmetry options.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `num_wiggles` - Number of wiggles
  - `total_length` - Total length of the meander
  - `start_bend_angle` - Starting bend angle
  - `symmetric` - Whether the meander is symmetric

### CPWRightJoint
Creates a sharp right angle in a CPW.
- **Key Parameters**:
  - `s` - Structure to add the component to
  - `CCW` - Whether the turn is counterclockwise
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap

### RightJointWiggles
Creates square wiggles for faster simulations.
- **Key Parameters**:
  - `s` - Structure to add the component to
  - `total_length` - Total length of the wiggles
  - `num_wiggles` - Number of wiggles
  - `radius` - Radius parameter

### ChannelWigglesByLength
Channel version of CPWWigglesByLength.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `num_wiggles` - Number of wiggles
  - `total_length` - Total length of the meander
  - `start_bend_angle` - Starting bend angle
  - `symmetric` - Whether the meander is symmetric
  - `channelw` - Width of the channel

### CPWWigglesByArea
Creates CPW wiggles which fill a specified area.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `length` - Length of the area
  - `width` - Width of the area
  - `start_up` - Whether to start with a CCW turn
  - `radius` - Radius of the bends

### CPWTee
Creates a Tee structure with padding.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `stub_length` - Length of the stub
  - `feed_length` - Length of the feed line
  - `flipped` - Whether the Tee is flipped
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap

## Capacitors and Couplers

### CPWGapCap
Creates a CPW gap capacitor (a gap in the center pin).
- **Key Parameters**:
  - `GAP` - Gap distance
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap
  - `capacitance` - Assumed capacitance value

### CPWInductiveShunt
Creates an inductive shunt with multiple fingers.
- **Key Parameters**:
  - `num_segments` - Number of finger segments
  - `segment_length` - Length of each finger
  - `segment_width` - Width of each finger
  - `segment_gap` - Gap between fingers
  - `taper_length` - Length of the taper
  - `pinw` - Width of the center pin
  - `inductance` - Assumed inductance value

### CPWFingerCap
Creates a CPW finger capacitor.
- **Key Parameters**:
  - `num_fingers` - Number of interdigitated fingers
  - `finger_length` - Length of each finger
  - `finger_width` - Width of each finger
  - `finger_gap` - Gap between fingers
  - `taper_length` - Length of the taper
  - `gapw` - Width of the gap
  - `capacitance` - Assumed capacitance value

### CPWFingerCapInside
Inside version of CPWFingerCap for two-layer capacitors.
- **Parameters**: Same as CPWFingerCap

### CPWLCoupler
Creates a structure coupled to a CPW via an L coupler, used for medium to high Q hangers.
- **Key Parameters**:
  - `coupler_length` - Length of the coupler
  - `separation` - Separation distance
  - `flipped` - Whether the coupler is flipped
  - `padding_type` - Type of padding
  - `pad_to_length` - Length to pad to
  - `pinw` - Width of the center pin
  - `gapw` - Width of the gap
  - `radius` - Radius of the bend
  - `capacitance` - Assumed capacitance value

### ChannelCouplerLayer
Channel version of CPWLCoupler for use in two-layer masks.
- **Parameters**: Similar to CPWLCoupler

### FingerCoupler
Creates a CPWTee plus finger capacitor.
- **Key Parameters**:
  - `structure` - Structure to add the component to
  - `cap_desc` - Capacitor description
  - `stub_length` - Length of the stub
  - `padding_length` - Length of padding
  - `flipped` - Whether the coupler is fl
