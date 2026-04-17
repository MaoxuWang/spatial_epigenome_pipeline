import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import threshold_local
from PIL import Image
import cv2
import pandas as pd
import imutils
import argparse
import math
import os
import xml.etree.ElementTree as ET
import json
from argparse import ArgumentParser


FLANK_X = 100  #  consistent with perl script
FLANK_Y = 100 
GRID_SIZE = 20 
SVG_SIZE = 1200


def parse_args(parser):
    parser.add_argument("-i", "--svg_path", required=True, help="input svg file")
    parser.add_argument("-r", "--img_res", type=int, 
        default = 1300,
        help="image resolution (the smaller one)")
    parser.add_argument('-o', '--outdir', type=str, default="output directory")
    parser.add_argument('--decode_file', required = True)
    parser.add_argument('--in_tissue', action='store_true',
        help="Enable tissue filtering (default: False)")

        
    return parser.parse_args()



def parse_svg_rects(svg_path):
    """Explain rect form SVG"""
    tree = ET.parse(svg_path)
    root = tree.getroot()
    ns = {'ns': 'http://www.w3.org/2000/svg'}
    
    rects = []
    for rect in root.findall('.//ns:rect', ns):
        # Exclude black background
        if rect.get('fill') == 'black' and \
           float(rect.get('width', 0)) == 990 and \
           float(rect.get('height', 0)) == 990:
            continue
            
        rect_data = {
            'x': float(rect.get('x', 0)),
            'y': float(rect.get('y', 0)),
            'opacity': float(rect.get('opacity', 1.0))
        }
        rects.append(rect_data)
    return rects

def create_opacity_matrix(rects):
    """Opacity matrix creation"""
    matrix = np.zeros((1200, 1200), dtype=np.uint8)
    for rect in rects:
        x = int(rect['x'])
        y = int(rect['y'])
        # fill opacity 10x10 square（0-255）
        matrix[y:y+10, x:x+10] = int(rect['opacity'] * 255)
    return matrix

def detect_tissue_region(opacity_matrix, in_tissue):
    """OpenCV detect ROI"""
    # Otsu automatic search for cutoff
    if not in_tissue:
        # all region valid except for non-captured area
        return opacity_matrix.astype(np.uint8)
    else:
        _, thresh = cv2.threshold(opacity_matrix, 0, 255, 
                                cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        
        # morphology continuous area
        kernel = np.ones((5,5), np.uint8)
        thresh = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
        thresh = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel)
    
    return thresh

def get_tissue_barcodes(rects, mask, in_tissue):
    """extract barcodes from mask"""
    barcodes = []
    for rect in rects:
        # calculate original coordinates (consistent with Perl script generation logic)
        x_orig = int((rect['x'] - FLANK_X) / GRID_SIZE) + 1
        y_orig = int((rect['y'] - FLANK_Y) / GRID_SIZE) + 1
        
        # calculate center coordinates
        cx = int(rect['x'] + 5)  
        cy = int(rect['y'] + 5)
        
        if (100 <= cx < 1090) and (100 <= cy < 1090):
            # convert to cropped coordinates
            cx_cropped = cx - 100
            cy_cropped = cy - 100
            # check whether coordinates are valid
            if (0 <= cx_cropped < 990) and (0 <= cy_cropped < 990):
                if mask[cy_cropped, cx_cropped] != 0:
                    barcodes.append(f"{x_orig}x{y_orig}")
    return barcodes


def main():
    """
        generate in-tissue or whole spatial metadata for spatial genomics
    """
    # Load parameters
    parser = ArgumentParser(description=main.__doc__)
    args = parse_args(parser)

    outdir = args.outdir
    decode_file = args.decode_file
    barcode_list = os.path.join(outdir, "barcode_ROI.txt")
    out_png = os.path.join(outdir, "tissue_lowres_image.png")
    in_tissue = args.in_tissue

    # Processing
    # Setp1: process svg file
    rects = parse_svg_rects(args.svg_path)

    # Step2: create opacity matrix
    opacity_matrix = create_opacity_matrix(rects)

    # Step3: Detect tissue region
    tissue_mask = detect_tissue_region(opacity_matrix, in_tissue)

    tissue_mask_cropped = tissue_mask[100:1090, 100:1090]
    # resize
    scale_factor = args.img_res / 990
    tissue_mask_resized = cv2.resize(tissue_mask_cropped, 
                                    (args.img_res, args.img_res),
                                    interpolation=cv2.INTER_NEAREST)
    cv2.imwrite(out_png, tissue_mask_resized)



    # Step4: get valid barcodes
    valid_barcodes = get_tissue_barcodes(rects, tissue_mask_cropped, in_tissue)

    print(f"Found {len(valid_barcodes)} tissue barcodes")
    print(f"Mask saved to {out_png}")

    # Convert the list into a DataFrame
    position = pd.DataFrame(valid_barcodes, columns=['Values'])

    # Add a new column 'in_tissue'
    position['in_tissue'] = 1

    # Set the row index to the values of the first column
    position.index = position['Values']

    # Read the spatial_barcodes.txt file
    BC = pd.read_csv(decode_file, sep='\t', header=None)

    # Create a new column 'position' by concatenating the values of the second and third columns
    BC['position'] = BC[1].astype(str) + 'x' + BC[2].astype(str)
    BC.loc[:, [0, "position"]].to_csv(os.path.join(outdir, "spatial_barcodes.txt"), index=False, header=True)

    # Set the row index to the values of the 'position' column
    BC.index = BC['position']

    # Merge the two dataframes
    df_in_tissue = pd.merge(position, BC, left_index=True, right_index=True, how='outer')

    # Select only the required columns
    df_in_tissue = df_in_tissue[[0, 'in_tissue', 2, 1]]

    # Replace NaN values with 0
    # df_in_tissue = df_in_tissue.fillna(0)
    # df_in_tissue['in_tissue'] = 1 if not in_tissue else df_in_tissue['in_tissue'].fillna(0)
    df_in_tissue.loc[df_in_tissue.index.isin(valid_barcodes), 'in_tissue'] = 1

    # Subtract 1 from the values of the third and fourth columns
    df_in_tissue[1] = df_in_tissue[1] - 1
    df_in_tissue[2] = df_in_tissue[2] - 1


    # 计算低分辨率图像坐标
    pixel_ratio = args.img_res / 990
    df_in_tissue[3] = (df_in_tissue[2] * GRID_SIZE + 1) * pixel_ratio  # x坐标
    df_in_tissue[4] = (df_in_tissue[1] * GRID_SIZE + 1) * pixel_ratio  # y坐标

    # 坐标取整
    df_in_tissue[3] = np.ceil(df_in_tissue[3])
    df_in_tissue[4] = np.ceil(df_in_tissue[4])

    # Create the directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    df_in_tissue = df_in_tissue.set_axis(['barcode', 'in_tissue', 'array_row', 'array_column', 'pxl_col_in_fullres', 'pxl_row_in_fullres'], axis=1)

    df_in_tissue.to_csv(os.path.join(outdir, "tissue_positions_list.csv"), index=False, header=False)
    df_in_tissue[df_in_tissue['in_tissue'] == 1].loc[:, ['barcode']].to_csv(os.path.join(outdir, "ROI_barcodes.txt"), index=False, header=False)
    # ## 3. Check_spots_alignment
    im = cv2.imread(out_png)
    # im_res = im
    im_res = cv2.resize(im, (args.img_res, args.img_res))

    spots = pd.read_csv(os.path.join(outdir, "tissue_positions_list.csv"), header=None)
    spots = spots.set_axis(['barcode', 'in_tissue', 'array_row', 'array_column', 'pxl_col_in_fullres', 'pxl_row_in_fullres'], axis=1)

    spot_size = 6 * pixel_ratio  # 原始10像素按比例缩放
    spot_width = int(spot_size)
    spot_height = int(spot_size)

    for i, spot in spots.iterrows():
        startX = int(round(spot['pxl_row_in_fullres']))
        startY = int(round(spot['pxl_col_in_fullres']))
        width = int(round(spot_width)*2)
        height = int(round(spot_height)*2)
        
        if(spot['in_tissue'] == 1):
            cv2.rectangle(im_res, (startX, startY), (startX+width, startY+height), (0, 255, 0), 1)
        
    cv2.imwrite(os.path.join(outdir, 'tissue_lowres_image_align_check.png'), im_res) 

    # ## 4. generate json
    # generate scalefactors_json.json

    scalefactors = {"spot_diameter_fullres": spot_size, 
                    "tissue_hires_scalef": 1.0, 
                    "fiducial_diameter_fullres": spot_size, 
                    "tissue_lowres_scalef": 1.0}
    with open(os.path.join(outdir, 'scalefactors_json.json'), 'w') as outfile:
        json.dump(scalefactors, outfile)


if __name__ == '__main__':
    main()