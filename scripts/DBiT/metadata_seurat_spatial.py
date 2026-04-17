from argparse import ArgumentParser
import pandas as pd
import xml.etree.ElementTree as ET
import os
import json
import sys 


def parse_args(parser):
    parser.add_argument('--decode_file', required = True)
    parser.add_argument('--svg', required = True)
    parser.add_argument('-o', '--outdir', type=str, default="output directory")

    return parser.parse_args()


def read_svg(svg_file):

    svg = ET.parse(svg_file)
    root = svg.getroot()
    return root 


def main():
    parser = ArgumentParser(description="Generate metadata information for spatial mapping")
    args = parse_args(parser)

    barcodes = pd.read_csv(args.decode_file, sep='\t',  lineterminator='\n', header=None)
    barcodes.columns = ['barcode', 'x', 'y']
    barcodes.index = barcodes['x'].astype(str) + 'x' + barcodes['y'].astype(str)
    barcodes = barcodes.drop(columns=['x', 'y'])
    barcodes['index'] = barcodes.index

    # Prepare tissue_positions_list.csv file that compatible with Seurat
    spots = pd.DataFrame(columns=['in_tissue', 
                                'array_row',
                                'array_column', 
                                'pxl_row_in_fullres',
                                'pxl_col_in_fullres'])

    col_reversed = False

    spot_width = 0
    spot_height = 0

    root = read_svg(args.svg)
    for i, g in enumerate(root.findall("{http://www.w3.org/2000/svg}rect")):
        if i == 0:
            continue
        x = float(g.attrib['x'])
        y = float(g.attrib['y'])
        spot_width = float(g.attrib['width'])
        spot_height = float(g.attrib['height'])
        in_tissue = 1
        # if(g.attrib['class'] == 'st1'):
        #     in_tissue = 1
        # else:
        #     in_tissue = 0
        
        row = int((x - 100) / 20) + 1
        col = int((y - 100) / 20) + 1
        if(col_reversed == True):
            col = 50-col+1
        
        x_c = x + spot_width/2
        y_c = y + spot_height/2
        
        spots = pd.concat([spots, pd.DataFrame({'in_tissue': in_tissue, 'array_row': row, 'array_column': col, 
                            'pxl_row_in_fullres': int(round(x_c)),                            
                            'pxl_col_in_fullres': int(round(y_c))}, index = [0])], ignore_index=True)
    spots.index = spots['array_column'].astype(str) + 'x' + spots['array_row'].astype(str)
    spots = pd.concat([barcodes.drop('index', axis=1), spots], axis=1)
    spots.to_csv(os.path.join(args.outdir, 'tissue_positions_list.csv'), index=False, header=False)

    # Generate scalefactors_json.json file
    scalefactors = {"spot_diameter_fullres": spot_width, 
                    "tissue_hires_scalef": 1.0, 
                    "fiducial_diameter_fullres": spot_width, 
                    "tissue_lowres_scalef": 1.0}

    with open(os.path.join(args.outdir, 'scalefactors_json.json'), 'w') as outfile:
        json.dump(scalefactors, outfile)
    print(f'metadata for seurat spatial object is ready in {args.outdir}!', file=sys.stderr)

if __name__ == '__main__':
    main()