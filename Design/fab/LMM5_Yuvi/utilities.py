'''
Goal: Provide utility functions which help make wafer in phidl 
'''
import phidl.geometry as pg
import os 
import numpy as np

def add_dxf_chip_to_wafer(wafer, dxf_filename, path, chip_width = 0, idx= 0, layer_in=0, 
                          layer_out=5, x_offset=0, y_offset=0, add_bool = True, 
                          cell_name = None, spacing = 100):
    """
    Imports a? GDS file converted from a DXF filename, copies a layer, and places it on the wafer.
    Args:
        wafer: phidl Device to add chip to
        dxf_filename: str, filename ending in .dxf
        path: str, directory containing the GDS files
        chip_width: int, spacing between chips
        idx: int, index for chip placement
        layer_in: int, layer to copy from
        layer_out: int, layer to copy to
        x_offset: int, additional x offset for chip placement
        y_offset: int, additional y offset for chip placement
        add_bool: bool, whether to add the chip to the wafer or not
    """
    gds_filename = dxf_filename.replace('.dxf', '.gds')
    file_path = os.path.join(path, gds_filename)
    print(f"Importing {file_path}")
    imported = pg.import_gds(file_path, cellname=cell_name)
    copied = pg.copy_layer(imported, layer=layer_in, new_layer=layer_out)
    x_pos = idx * (chip_width) + (x_offset)
    y_pos = y_offset#+200

    if add_bool:
        wafer << copied.move((x_pos, y_pos))
    return wafer, copied


def add_dxf_chip_to_wafer_junc(wafer, dxf_filename, path, chip_width, idx, gap_layer_in=np.int64(2), pin_layer_in=np.int64(3),
                               dicing_layer_in=np.int64(4),
                          gap_layer_out=2, pin_layer_out = 3, dicing_layer_out=4,
                            x_offset=0, y_offset=0, add_bool = True):
    """
    Imports a? GDS file converted from a DXF filename, copies a layer, and places it on the wafer.
    Args:
        wafer: phidl Device to add chip to
        dxf_filename: str, filename ending in .dxf
        path: str, directory containing the GDS files
        chip_width: int, spacing between chips
        idx: int, index for chip placement
        layer_in: int, layer to copy from
        layer_out: int, layer to copy to
        x_offset: int, additional x offset for chip placement
        y_offset: int, additional y offset for chip placement
        add_bool: bool, whether to add the chip to the wafer or not
    """
    gds_filename = dxf_filename.replace('.dxf', '.gds')
    file_path = os.path.join(path, gds_filename)
    print(f"Importing {file_path}")
    imported = pg.import_gds(file_path)
    # copied = pg.copy_layer(imported, layer=layer_in, new_layer=layer_out)
    print("Layers in imported DXF:", imported.layers)
    copied_gap = pg.copy_layer(imported, layer=gap_layer_in, new_layer=gap_layer_out)
    # copied = pg.copy_layer(copied, layer=pin_layer_in, new_layer=pin_layer_out)
    copied_pin = pg.copy_layer(imported, layer=pin_layer_in, new_layer=pin_layer_out)
    copied_dicing = pg.copy_layer(imported, layer=dicing_layer_in, new_layer=dicing_layer_out)
    
    x_pos = idx * chip_width + x_offset
    y_pos = y_offset

    # make a bounding box for the chip
    if add_bool:
        # wafer << copied.move((x_pos, y_pos))
        wafer << copied_gap.move((x_pos, y_pos))
        wafer << copied_pin.move((x_pos, y_pos))
        wafer << copied_dicing.move((x_pos, y_pos))
    print('returning dcing layer too')
    return wafer, copied_gap, copied_pin, copied_dicing

# for multimode qubit 
# inverts qubit and add bounding box for laser dicing
def add_qubit_bbox_and_negative_mask(wafer, qubit_chip_obj, qubit_chip_width, idx, xpos_offset, ypos_offset):
    qubit_length = 16550  # Example length for qubit chip, adjust as needed
    test_structure_length_down = 9850
    test_structure_length_up = 2100

    xpos = xpos_offset + (idx * qubit_chip_width) + 200
    ypos = ypos_offset + 200

    # qubit_bbox = pg.rectangle(size=(qubit_chip_width, qubit_length), layer=13)
    # wafer << qubit_bbox.move((xpos -200, ypos-200  + test_structure_length_down))
    qubit_bbox = pg.rectangle(size=(qubit_chip_width, qubit_length + test_structure_length_down - 250), layer=13)
    wafer << qubit_bbox.move((xpos -200, ypos-200+100))
    # test_structure_bbox = pg.rectangle(size=(qubit_chip_width, test_structure_length_down), layer=13)
    # wafer << test_structure_bbox.move((xpos - 200, ypos - 200))

    qubit_chip_bbox = pg.rectangle(size=(qubit_chip_width, qubit_length + test_structure_length_down + test_structure_length_up), layer=14)
    qubit_chip_bbox.move((-200, -200))

    qubit_negative = pg.boolean(A=qubit_chip_bbox, B=qubit_chip_obj, operation='xor', layer=5)
    wafer << qubit_negative.move((xpos, ypos))
    # wafer<< qubit_chip_obj.move((xpos, ypos))

    return wafer, qubit_negative, qubit_bbox


# import phidl.device as pd
# import os

# Assume these helper functions are defined elsewhere in your code
# from your_helpers import add_dxf_chip_to_wafer, add_qubit_bbox_and_negative_mask
# from your_helpers import add_dicing_lanes, add_dxf_chip_to_wafer_junc

# def add_mm_qubit_array(
#     wafer,
#     dxf_filenames,
#     dxf_path,
#     layout_mode,
#     top_left_pos=(1000, 25870),
#     spacing=100,
#     chip_width=2000,
#     chip_length=28497.05,
#     chip_top_dead_space=2250,
#     dicing_lane_margin=100,
# ):
#     """
#     Imports an array of qubit chips from DXF files, places them on a wafer,
#     rotates them, and adds associated geometry like bounding boxes and dicing lanes.

#     This function handles the specific coordinate transformation where chips are
#     designed vertically and then rotated -90 degrees into a horizontal array.

#     Args:
#         wafer (Device): The main PHIDL Device to which the chips will be added.
#         dxf_filenames (list of str): A list of the DXF file names to import.
#         dxf_path (str): The path to the directory containing the DXF files.
#         layout_mode (str): Determines the geometry to add. Can be 'optical' or 'ebeam'.
#         top_left_pos (tuple): The (x, y) coordinate for the top-left corner of the
#                               *final, rotated* array of chips.
#         spacing (float): Spacing between adjacent chips in the array.
#         chip_width (float): The width of a single qubit chip.
#         chip_length (float): The total length of a single qubit chip DXF.
#         chip_top_dead_space (float): Unused space at the top of the chip to be excluded
#                                      from some calculations (e.g., dicing).
#         dicing_lane_margin (float): Additional margin for dicing lane calculations.

#     Returns:
#         Device: The modified wafer Device with the qubit array added.

#     Requires the following helper functions to be defined:
#         - add_dxf_chip_to_wafer()
#         - add_qubit_bbox_and_negative_mask()
#         - add_dicing_lanes()
#         - add_dxf_chip_to_wafer_junc()
#     """
#     # The original code defines the final desired coordinates and then calculates
#     # the necessary pre-rotation GDS coordinates. We do the same here.
#     # Final desired position (x_final, y_final) is achieved by rotating the
#     # GDS position (x_gds, y_gds) by -90 degrees.
#     # Rotation formula for -90 deg: (x_final, y_final) = (y_gds, -x_gds)
#     # Therefore, y_gds = x_final and x_gds = -y_final.
#     x_gds_offset = -top_left_pos[1]
#     y_gds_offset = top_left_pos[0]

#     print(f"Calculated pre-rotation GDS offsets: x={x_gds_offset}, y={y_gds_offset}")

#     if layout_mode == 'optical':
#         # Create and place each chip and its bounding box
#         for idx, file in enumerate(dxf_filenames):
#             # This helper should add the DXF to the wafer and return a reference to it
#             wafer, qubit_obj = add_dxf_chip_to_wafer(
#                 wafer, file, dxf_path, chip_width, idx,
#                 x_offset=x_gds_offset, y_offset=y_gds_offset, add_bool=False
#             )
#             # This helper creates the negative mask and bounding box
#             wafer, qubit_negative, qubit_bbox = add_qubit_bbox_and_negative_mask(
#                 wafer, qubit_obj, chip_width, idx,
#                 xpos_offset=x_gds_offset + idx * spacing, ypos_offset=y_gds_offset
#             )
            
#             # Rotate the added components into their final position
#             # Note: Depending on your helpers, you might need to rotate qubit_obj too.
#             # Assuming the helpers add references to the wafer that can be rotated.
#             qubit_negative.rotate(-90)
#             qubit_bbox.rotate(-90)

#         # Add dicing lanes for the entire array
#         num_chips = len(dxf_filenames)
#         box_size = (chip_width, chip_length - chip_top_dead_space - dicing_lane_margin)
        
#         # Calculate position for the dicing lane geometry before rotation
#         x_pos_dicing_gds = x_gds_offset - spacing / 2
#         y_pos_dicing_gds = y_gds_offset
        
#         print(f"Dicing lane pre-rotation GDS offsets: x={x_pos_dicing_gds}, y={y_pos_dicing_gds}")

#         dicing_lane_geometry = add_dicing_lanes(
#             x_start=x_pos_dicing_gds,
#             y_start=y_pos_dicing_gds,
#             num_cols=num_chips,
#             col_spacing=spacing,
#             box_size=box_size,
#             num_rows=1 # Based on the original code
#         )
#         dicing_lane_geometry.rotate(-90)
#         wafer << dicing_lane_geometry

#     elif layout_mode == 'ebeam':
#         for idx, file in enumerate(dxf_filenames):
#             # This helper adds the e-beam layers (gap, pin, dicing)
#             wafer, qubit_gap, qubit_pin, qubit_dicing = add_dxf_chip_to_wafer_junc(
#                 wafer, file, dxf_path, chip_width, idx,
#                 x_offset=x_gds_offset, y_offset=y_gds_offset, add_bool=True
#             )
            
#             # Rotate the components into their final position
#             qubit_gap.rotate(-90)
#             qubit_pin.rotate(-90)
#             qubit_dicing.rotate(-90)
#     else:
#         raise ValueError(f"Invalid layout_mode: '{layout_mode}'. Must be 'optical' or 'ebeam'.")

#     return wafer
# ----------------- coherence qubit utilities -----------------
def add_coherence_qubits(
    wafer_device,
    x_offset,
    y_offset,
    qubit_gds=None,
    bbox_size=(7000, 7000),
    spacing=None,
    num_qubits=1,
    layer_in = 5, 
    layer_out = 5,
    bbox = False,
    origin_at_chip_botleft_corner = False,
):
    if spacing is None:
        spacing = 7000

    # Initialize arrays for X and Y coordinates
    blank_x_arr = np.array([0] * num_qubits)
    blank_y_arr = np.array([0] * num_qubits)

    blank_x_arr = blank_x_arr  + x_offset
    blank_y_arr = blank_y_arr + y_offset
    # blank_x_arr = #[blank_x_arr[i] + (i) * spacing for i in range(num_qubits)]
    blank_y_arr = [blank_y_arr[i] + (i) * spacing for i in range(num_qubits)]

    for i in range(num_qubits):
        # qubit_raw = pg.import_gds(qubit_gds)
        qubit_loc = (blank_x_arr[i] + bbox_size[0] / 2, blank_y_arr[i] + bbox_size[1] / 2)
        print("Qubit location before adjustment:", qubit_loc)
        print("BBox size:", bbox_size)
        bbox_loc = (blank_x_arr[i] + bbox_size[0] / 2 , blank_y_arr[i] + bbox_size[1] / 2)
        if origin_at_chip_botleft_corner:
            qubit_loc = (qubit_loc[0] - 1 * bbox_size[0] / 2, qubit_loc[1] - 1 * bbox_size[1] / 2)
            print("Adjusted qubit location for bottom-left origin:", qubit_loc)

        qubit_raw = import_gds_robust(qubit_gds, spacing=0)
        qubit_pos = pg.copy_layer(qubit_raw, layer=layer_in, new_layer=layer_out)
        qubit_bbox = pg.rectangle(size=(bbox_size[0], bbox_size[1]), layer=13).move((-bbox_size[0] / 2, -bbox_size[1] / 2))

        wafer_device << qubit_pos.move(qubit_loc)
        # print(f"Placed qubit at x: {blank_x_arr[i] + bbox_size[0] / 2}, y: {blank_y_arr[i] + bbox_size[1] / 2}")
        if bbox:
            wafer_device << qubit_bbox.move(bbox_loc)

def place_chip_row(
    wafer: pg.Device,
    gds_files: list,
    start_position: tuple,
    chip_size: tuple = (6000, 9000),
    spacing: float = 100,
    layers: tuple = (200, 5),
    origin_at_chip_botleft_corner: bool = False,
    bbox: bool = True
):
    """
    Places a horizontal row of chips on a wafer and adds dicing lanes.

    Args:
        wafer (pd.Device): The main PHIDL Device to which geometry will be added.
        gds_files (list): A list of strings, where each string is the
                          identifier for a chip GDS file.
        start_position (tuple): The (x, y) coordinate for the bottom-left
                                corner of the first chip in the row.
        chip_size (tuple): The (width, height) of a single chip.
        spacing (float): The spacing between adjacent chips.
        layers (tuple): A tuple containing the (layer_in, layer_out) to be
                        passed to the chip-adding function.
    """
    num_chips = len(gds_files)
    if num_chips == 0:
        print("Warning: No GDS files provided to place_chip_row.")
        return

    x_start, y_start = start_position
    chip_width, _ = chip_size
    layer_in, layer_out = layers
    # if layers is 2D list, handle that
    separate_layers_bool = False
    print(layers)
    if isinstance(layer_in, list) and isinstance(layer_out, list):
        print("Using separate layers for each chip.")
        separate_layers_bool = True
        # don't use layer in , layer out, just use layers[i] in the loop below
    
    # Calculate the horizontal pitch (center-to-center or edge-to-edge distance)
    pitch_x = chip_width + spacing

    # --- 1. Place all the chips in a row ---
    print(f"Placing a row of {num_chips} chips...")
    for i, file_name in enumerate(gds_files):
        if file_name is None or file_name == "":
            print(f"Skipping chip {i+1}/{num_chips} due to empty filename.")
            continue
        # Calculate the position for the current chip
        current_x = x_start + i * pitch_x
        print(f"Placing chip {i+1}/{num_chips} at x={current_x}, y={y_start} using file: {file_name}")
        
        # Call your function to add the actual chip geometry
        add_coherence_qubits(
            wafer_device=wafer, 
            x_offset=current_x,
            y_offset=y_start,
            qubit_gds=file_name,
            bbox_size=chip_size,
            # Pass along other necessary parameters from your original code
            num_qubits=1, # This function places one chip at a time
            bbox=bbox, 
            layer_in=layers[i][0] if separate_layers_bool else layer_in, 
            layer_out=layers[i][1] if separate_layers_bool else layer_out,
            origin_at_chip_botleft_corner=origin_at_chip_botleft_corner
        )

    # --- 2. Add dicing lanes around the entire row ---
    # Dicing lanes should start 'spacing' distance away from the first chip's edge
    dicing_grid_x_start = x_start - spacing/2
    dicing_grid_y_start = y_start - spacing/2

    if bbox:
    
        dicing_lane_geometry = add_dicing_lanes(
            box_x_offset=dicing_grid_x_start,
            box_y_offset=dicing_grid_y_start,
            num_rows=1,
            num_cols=num_chips,
            box_size=chip_size,
            spacing=spacing
        )
        wafer << dicing_lane_geometry

def place_chip_column(
    wafer: pg.Device,
    gds_files: list,
    start_position: tuple,
    chip_size: tuple = (6000, 9000),
    spacing: float = 100,
    layers: tuple = (200, 5),
    origin_at_chip_botleft_corner: bool = False,
    bbox: bool = True
):
    """
    Places a vertical column of chips on a wafer and adds dicing lanes.

    Args:
        wafer (pd.Device): The main PHIDL Device to which geometry will be added.
        gds_files (list): A list of strings, where each string is the
                          identifier for a chip GDS file.
        start_position (tuple): The (x, y) coordinate for the bottom-left
                                corner of the first (lowest) chip in the column.
        chip_size (tuple): The (width, height) of a single chip.
        spacing (float): The spacing between adjacent chips.
        layers (tuple): A tuple containing the (layer_in, layer_out) to be
                        passed to the chip-adding function.
    """
    num_chips = len(gds_files)
    if num_chips == 0:
        print("Warning: No GDS files provided to place_chip_row.")
        return

    x_start, y_start = start_position
    _, chip_length = chip_size
    layer_in, layer_out = layers
    
    # Calculate the horizontal pitch (center-to-center or edge-to-edge distance)
    pitch_y = chip_length + spacing

    # --- 1. Place all the chips in a row ---
    print(f"Placing a row of {num_chips} chips...")
    for i, file_name in enumerate(gds_files):
        if file_name is None or file_name == "":
            print(f"Skipping chip {i+1}/{num_chips} due to empty filename.")
            continue
        # Calculate the position for the current chip
        current_y = y_start + i * pitch_y
        print(f"Placing chip {i+1}/{num_chips} at x={x_start}, y={current_y} using file: {file_name}")

        # Call your function to add the actual chip geometry
        add_coherence_qubits(
            wafer_device=wafer, 
            x_offset=x_start,
            y_offset=current_y,
            qubit_gds=file_name,
            bbox_size=chip_size,
            # Pass along other necessary parameters from your original code
            num_qubits=1, # This function places one chip at a time
            bbox=bbox, 
            layer_in=layer_in, 
            layer_out=layer_out,
            origin_at_chip_botleft_corner=origin_at_chip_botleft_corner
        )

    # --- 2. Add dicing lanes around the entire row ---
    # Dicing lanes should start 'spacing' distance away from the first chip's edge
    dicing_grid_x_start = x_start - spacing/2
    dicing_grid_y_start = y_start - spacing/2
    if bbox:
        dicing_lane_geometry = add_dicing_lanes(
            box_x_offset=dicing_grid_x_start,
            box_y_offset=dicing_grid_y_start,
            num_rows=num_chips,
            num_cols=1,
            box_size=chip_size,
            spacing=spacing, 
        )
        wafer << dicing_lane_geometry

# ----------------- dicing lane utilities -----------------
from phidl import Device
import phidl.geometry as pg

# Dummy create_boxes function for completeness
def create_boxes(num_rows, num_cols, box_size, spacing, wafer, x_offset=0, y_offset=0):
    array_width = num_cols * (box_size[0] + spacing) - spacing
    array_height = num_rows * (box_size[1] + spacing) - spacing
    start_x = x_offset - array_width / 2
    start_y = y_offset - array_height / 2
    for row in range(num_rows):
        for col in range(num_cols):
            box_geo = pg.rectangle(size=box_size, layer=99)
            box_ref = wafer << box_geo
            box_ref.xmin = start_x + col * (box_size[0] + spacing)
            box_ref.ymin = start_y + row * (box_size[1] + spacing)
    return wafer

# This function is correct and creates a centered cross.
def create_cross_device(length, linewidth):
    """
    Creates a PHIDL Device containing a single, centered cross shape.
    """
    D = Device('cross')
    vertical_bar = pg.rectangle(size=(linewidth, length))
    vertical_bar.center = (0, 0)
    D << vertical_bar
    horizontal_bar = pg.rectangle(size=(length, linewidth))
    horizontal_bar.center = (0, 0)
    D << horizontal_bar
    return D


def add_dicing_lanes(box_x_offset, box_y_offset, num_cols, spacing, box_size, num_rows, positive_cross_no_lanes = False, 
                     dicing_layer = 5, crosses_dist = 1000):
    """
    Creates dicing lanes with two sets of crosses:
    1. One cross precisely at each lane intersection (chip corner).
    2. Crosses at 1mm intervals along the lanes, between the corners.

    box_x_offset, box_y_offset now correspond to the lower-left corner of the first box
    (i.e. the corner of the first box), not the center of the whole array.
    """
    lanes_to_expose = Device('LANES')
    crosses_to_subtract = Device('CROSSES')

    # --- 1. Define lane geometry ---
    total_width = num_cols * (box_size[0] + spacing) - spacing
    total_height = num_rows * (box_size[1] + spacing) - spacing
    dicing_lane_width = 100

    full_lane_width = total_width + spacing
    full_lane_height = total_height + spacing

    # Interpret box_x_offset/box_y_offset as the corner (lower-left) of the first box
    left_edge_x = box_x_offset
    bottom_edge_y = box_y_offset

    # vertical lanes at every box boundary (num_cols + 1 boundaries)
    vertical_lane_x_coords = [left_edge_x + i * (box_size[0] + spacing) for i in range(num_cols + 1)]
    # horizontal lanes at every box boundary (num_rows + 1 boundaries)
    horizontal_lane_y_coords = [bottom_edge_y + i * (box_size[1] + spacing) for i in range(num_rows + 1)]

    # centers for lane rectangles
    center_x = left_edge_x + total_width / 2
    center_y = bottom_edge_y + total_height / 2

    # Create visual lane rectangles
    for x in vertical_lane_x_coords:
        lane_geo = pg.rectangle(size=(dicing_lane_width, full_lane_height), layer=5)
        lane_ref = lanes_to_expose << lane_geo
        lane_ref.center = (x, center_y)

    for y in horizontal_lane_y_coords:
        lane_geo = pg.rectangle(size=(full_lane_width, dicing_lane_width), layer=5)
        lane_ref = lanes_to_expose << lane_geo
        lane_ref.center = (center_x, y)

    # --- 2. Place crosses based on the two rules ---

    # RULE 1: Place one cross at every corner intersection.
    for x in vertical_lane_x_coords:
        for y in horizontal_lane_y_coords:
            cross_ref = crosses_to_subtract << create_cross_device(100, 20)
            cross_ref.center = (x, y)

    # RULE 2: Place crosses at 1mm intervals BETWEEN corners.

    # Along Vertical Lanes (between bottom_edge_y and bottom_edge_y + total_height)
    y_start = bottom_edge_y
    num_y_steps = int(total_height / crosses_dist)
    for x in vertical_lane_x_coords:
        for j in range(1, num_y_steps):  # skips the start and end points
            y = y_start + j * crosses_dist
            cross_ref = crosses_to_subtract << create_cross_device(100, 20)
            cross_ref.center = (x, y)

    # Along Horizontal Lanes (between left_edge_x and left_edge_x + total_width)
    x_start = left_edge_x
    num_x_steps = int(total_width / 
                      crosses_dist)
    for y in horizontal_lane_y_coords:
        for j in range(1, num_x_steps):  # skips the start and end points
            x = x_start + j * crosses_dist
            cross_ref = crosses_to_subtract << create_cross_device(100, 20)
            cross_ref.center = (x, y)

    # --- 3. Perform the final boolean subtraction ---
    final_lanes = pg.boolean(
        A=lanes_to_expose,
        B=crosses_to_subtract,
        operation='not',
        layer=dicing_layer
    )
    if positive_cross_no_lanes:
        final_lanes = pg.boolean(
            A=crosses_to_subtract,B = crosses_to_subtract, operation='or', layer=dicing_layer
        )
    return final_lanes

# ----------------- import gds utilities -----------------
# from phidl import Device as pd
import gdstk  # PHIDL's backend, used for reading GDS library info
import phidl as p 


def import_gds_robust(
    filename: str, 
    cellname: str = None,
    spacing: float = 20, 
    **kwargs
) -> Device:
    """
    A robust wrapper for phidl.geometry.import_gds that handles the case
    of multiple top-level cells by packing them into a single Device.

    Args:
        filename (str): Path to the GDS file.
        cellname (str, optional): If specified, imports only this cell.
                                  Defaults to None.
        spacing (float): Spacing used by p.pack() if multiple cells are packed.
        **kwargs: Additional keyword arguments passed to pg.import_gds().

    Returns:
        Device: The imported geometry as a PHIDL Device.
    """
    # --- Case 1: The user specified a cell. This is the most direct path. ---
    if cellname is not None:
        print(f"[import_gds_robust] Importing specific cell: '{cellname}'")
        return pg.import_gds(filename=filename, cellname=cellname, **kwargs)

    # --- Case 2: Auto-detection. Try the standard import first. ---
    try:
        # This is the "happy path" for clean files with 0 or 1 top-level cell.
        print("[import_gds_robust] Attempting standard import...")
        return pg.import_gds(filename=filename, **kwargs)
    except ValueError as e:
        # If the happy path fails, check if it's the specific error we can handle.
        if "multiple top-level cells" in str(e):
            print("[import_gds_robust] Multiple cells detected. Switching to pack mode.")
            
            # Use gdstk *only* to get the list of names. No fixing.
            lib = gdstk.read_gds(filename)
            top_cell_names = [cell.name for cell in lib.top_level()]
            print(f"[import_gds_robust] Top-level cells found: {top_cell_names}")
            return pg.import_gds(filename=filename, cellname="chip", **kwargs)
            
            # if not top_cell_names:
            #     print("[import_gds_robust] Warning: File has no top-level cells.")
            #     return Device(f"EMPTY_{os.path.basename(filename)}")

            # print(f"[import_gds_robust] Found cells to pack: {top_cell_names}")
            
            # devices_to_pack = []
            # for name in top_cell_names:
            #     # Import each cell one by one using the standard function
            #     d = pg.import_gds(filename=filename, cellname=name, **kwargs)
            #     devices_to_pack.append(d)
            
            # # Create a container and pack the imported devices into it
            # container = Device(f"PACKED_{os.path.basename(filename).split('.')[0]}")
            # references = p.pack(device_list=devices_to_pack, spacing=spacing)
            # container.add(references)
            
            # return container
        else:
            # If it's a different ValueError (or any other error), we can't handle it.
            # Re-raise it so the user knows about the different problem.
            raise e