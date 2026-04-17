import argparse
import os
import sys
import cv2
import tifffile


if __name__ == '__main__':
    print(sys.argv)
    parser = argparse.ArgumentParser(description="Generate fluoresence")

    parser.add_argument('--seg_channel', type=int, help="seg_channel", required=True)
    parser.add_argument('--ori_img_path', type=str, help="FL ori image ", required=True)
    parser.add_argument('--out_dir', type=str, help="out image dir", required=True)
    args = parser.parse_args()

    img_path = args.ori_img_path
    seg_channel = args.seg_channel
    out_dir = args.out_dir
    print("img_path:{}\nseg_channel:{}\nout_dir:{}".format(img_path, seg_channel, out_dir))

    img = tifffile.imread(img_path)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if len(img.shape) == 2:
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    tifffile.imwrite(os.path.join(out_dir, "fluorescence.tif"), cv2.cvtColor(img[..., seg_channel], cv2.COLOR_GRAY2BGR)
                     , compression="jpeg")
