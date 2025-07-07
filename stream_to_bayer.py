#!/usr/bin/env python3

# --------------------------------------------------------------------------------
#
# stream_to_bayer.py
#
# Reference software to demultiplex ODR lossless 
# streams from the output of msss_odr.py into Bayer CCD form
#
# Release 2025-04-01
#
# MIT License
#
# Copyright (c) 2025 Malin Space Science Systems
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
import numpy as np
import argparse


#-------------------------------------------------------------
#
# This function interleaves stream order from decode_dat.py
# into natural CCD Bayer order
#
# The CCD CFA Bayer pattern
#
#   Line 
#   0    R0 G0 R0 G0 R0 G0 ...
#   1    G1 B1 G1 B1 G1 B1 ...
#   2    R2 G2 R2 G2 R2 G2 ...
#   3    G3 B3 G3 B3 G3 B3 ...
#   4    R4 G4 R4 G4 R4 G4 ...
#   5    G5 B5 G5 B5 G5 B5 ...
#   6    R6 G6 R6 G6 R6 G6 ...
#   7    G7 B7 G7 B7 G7 B7 ...
#
#  The output of decode_dat.py is de-interleaved
#  such that like-kind Bayer bands are consecutive
#
#   R0 R2 R4 R6
#   G0 G2 G4 G6
#   G1 G3 G5 G7
#   B1 B3 B5 B7
#   R8 R8 ...
#
#



def raw_to_bayer( raw_lines:np.array,  width:int, height:int ):

  # groups of 4

  LINES_PER_GROUP = 4         

  # each band will be a quarter image

  BAND_WIDTH  = width//2     
  BAND_HEIGHT = height//2   

  # reshape 1D stream to
  #   RRRR
  #   GGGG
  #   GGGG
  #   BBBB

  raw_lines = raw_lines.reshape( ( BAND_HEIGHT, BAND_WIDTH*LINES_PER_GROUP) )

  # slice out each band

  R   = raw_lines[0::LINES_PER_GROUP].reshape(BAND_HEIGHT,BAND_WIDTH)
  G1  = raw_lines[1::LINES_PER_GROUP].reshape(BAND_HEIGHT,BAND_WIDTH)
  G2  = raw_lines[2::LINES_PER_GROUP].reshape(BAND_HEIGHT,BAND_WIDTH)
  B   = raw_lines[3::LINES_PER_GROUP].reshape(BAND_HEIGHT,BAND_WIDTH)

  # Interleave to CFA pattern
  #   R G
  #   G B

  image = np.zeros( (height,width), np.uint8 )
  image[0::2,0::2] = R
  image[1::2,0::2] = G1
  image[0::2,1::2] = G2
  image[1::2,1::2] = B

  return image


def main():

  parser = argparse.ArgumentParser(description='Interleave raw data into a Bayer pattern image')
  parser.add_argument('infile', help='Input file path (raw data from decode_dat.py)')
  parser.add_argument('outfile', help='Output file path (Bayer image)')
  parser.add_argument('width', type=int, help='Width of the image')
  parser.add_argument('height', type=int, help='Height of the image')

  args = parser.parse_args()

  try:
    with open(args.infile, 'rb') as f:
      raw_data = np.fromfile(f, dtype=np.uint8)  # Assuming uint8 data, adjust if needed

  except FileNotFoundError:
    print(f'Error: Input file "{args.infile}" not found')
    return
  except Exception as e:
    print(f'An error occurred while reading the input file: {e}')
    return

  try:
    bayer_image = raw_to_bayer(raw_data, height=args.height, width=args.width)
    bayer_image.tofile( args.outfile )  
    print(f'Image processed and saved to "{args.outfile}"')

  except ValueError as e:
    print( e )
    print('Error: Input data size does not match specified width and height')
  except Exception as e:
    print(f'An error occurred during processing: {e}')

if __name__ == '__main__':
    main()


