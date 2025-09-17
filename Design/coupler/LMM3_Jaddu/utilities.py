'''
Goal: Provide utility functions which help make wafer in phidl 
'''
import phidl.geometry as pg
import os 
import numpy as np

def add_dxf_chip_to_wafer(wafer, dxf_filename, path, chip_width = 0, idx= 0, layer_in=0, 
                          layer_out=5, x_offset=0, y_offset=0, add_bool = True):
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
    copied = pg.copy_layer(imported, layer=layer_in, new_layer=layer_out)
    x_pos = idx * chip_width + x_offset
    y_pos = y_offset

    # make a bounding box for the chip
    if add_bool:
        wafer << copied.move((x_pos, y_pos))
    return wafer, copied


def add_dxf_chip_to_wafer_junc(wafer, dxf_filename, path, chip_width, idx, gap_layer_in=np.int64(2), pin_layer_in=np.int64(3),
                          gap_layer_out=2, pin_layer_out = 3,  x_offset=0, y_offset=0, add_bool = True):
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
    
    x_pos = idx * chip_width + x_offset
    y_pos = y_offset

    # make a bounding box for the chip
    if add_bool:
        # wafer << copied.move((x_pos, y_pos))
        wafer << copied_gap.move((x_pos, y_pos))
        wafer << copied_pin.move((x_pos, y_pos))
    return wafer, copied_gap, copied_pin