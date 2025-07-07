#!/usr/bin/env python3

# --------------------------------------------------------------------------------
#
# msss_odr.py
#
# Reference software to decode an MSSS ODR MiniHeader and decompress 
# an MSSS Original Data Record Payload
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
import argparse
from bitstring import ConstBitStream as bitstream
import huffman_table as huffman

#-------------------------------------------------------------
#
# Original Data Record MSSS Mini-Header 
# Bit structure breakdown referencing SIS section 4.4
# 

# Table 4.4.3

start_fields = [
  ( 'camera_product_id',      'uint:32' ),
  ( 'magic_initial',          'hex:32'  ),     # always 0xFF00F0CA
  ( 'spacecraft_clock_start', 'uint:32' ),    # DEA sclk
  ( 'detector_erase_count',   'uint:16' ),    # vertical register flush
  ( None,                     'pad:4'   ),      # padding
  ( 'CCD_state',              'uint:4'  ),     # table 4.4-4
  ( 'led1_control',           'uint:1'  ),     # 0 off, 1 on
  ( 'led2_control',           'uint:1'  ),     # 0 off, 1 on
  ( 'led3_control',           'uint:1'  ),     # 0 off, 1 on
  ( 'video_readout_mode',     'uint:1'  ),     # 0 off, 1 on
  ( 'clkdiv2',                'uint:1'  ),     # Table 4.4-4
  ( 'long_integration_mode',  'uint:1'  ),     # 0 off, 1 on
  ( 'test_mode',              'uint:1'  ),     # 0 off, 1 on
  ( 'clkdiv1',                'uint:1'  ),     # Table 4.4-4
]

image_fields = [
  ( 'filter_number',          'uint:8'  ),     # 0-7 filter index
  ( 'exposure_duration',      'uint:24' ),     # exposure in ms*10
  ( 'first_sample',           'uint:8'  ),     # start x/8
  ( 'first_line',             'uint:8'  ),     # start y/8
  ( 'samples',                'uint:8'  ),     # width/8
  ( 'lines',                  'uint:8'  ),     # height/8
]


# auto-focus

af_fields = [
  ( 'autofocus_init_position', 'uint:15' ),
  ( 'autofocus_step_size',     'uint:10' ),
  ( 'autofocus_step_count',    'uint:6'  ),
  ( 'autofocus_zstack',        'uint:1'  ),
]


# auto-exposure

ae_fields = [
  ( 'auto_exposure_data_cut',            'uint:8'  ),
  ( 'auto_exposure_pixel_fraction',      'uint:8'  ),
  ( 'auto_exposure_percent',             'uint:8'  ),
  ( 'auto_exposure_max_iteration_count', 'uint:8'  )
]


# zstack command fields
# overlap af_fields and ae_fields of image

zstack_fields = [
  ( 'focus_merge_start_camera_product_id',  'uint:32' ),
  ( 'focus_merge_depth',                    'uint:8'  ),
  ( None,                                   'pad:22'  ), # parms padding
  ( 'focus_merge_blending_flag',            'uint:1'  ), # 0 off, 1 on
  ( 'focus_merge_registration_flag',        'uint:1'  )  # 0 off, 1 on
]


# image compression, see SIS 4.4.4

compression_fields = [
  ( None,                     'pad:8'   ),   # compression A
  ( None,                     'pad:8'   ),   # compression B
  ( 'image_compression_mode', 'uint:8'  ),   # compression C
  ( 'image_quality',          'uint:8'  ),   # compression D
  ( None,                     'pad:24'  ),   # compression EFG
  ( 'companding_table',       'uint:8'  )    # compression H
]

flags_status_fields = [
  ( None,                       'pad:1'   ),   # padding
  ( 'led_uv',                   'uint:1'  ),   # led_bit 2
  ( 'led_vis2',                 'uint:1'  ),   # led_bit 1
  ( 'led_vis1',                 'uint:1'  ),   # led_bit 0
  ( None,                       'pad:1'   ),   # undefined
  ( 'hall_filter',              'uint:1'  ),   # hall sensor state for Mastcam filter wheel
  ( 'hall_cover',               'uint:1'  ),   # hall sensor state for MAHLI cover
  ( 'hall_focus',               'uint:1'  ),   # hall sensor state for focus mechanism
  ( 'instrument_serial_number', 'uint:24' ),   # DEA serial number
]


# mechanism motor position

mech_fields = [
  ( 'focus_position_count',  'int:32'  ),   # focus motor count
  ( 'zoom_position_count',   'int:16'  ),   # zoom motor count
  ( 'filter_position_count', 'int:16'  ),   # filter wheel motor count
]

miniheader_part2_fields = [
  ( 'dc_offset',          'uint:32' ), 
  ( 'intial_size',        'uint:32' ),
  ( 'magic_final',        'hex:32'  )    # always 0x1010CC28
]


# fields combined to decode  an image product
image_product_fields = start_fields + image_fields + af_fields + ae_fields + compression_fields +  \
  flags_status_fields + mech_fields + miniheader_part2_fields

# fields combined to decode a zstack product
zstack_product_fields = start_fields + image_fields + zstack_fields + compression_fields +  \
  flags_status_fields + mech_fields + miniheader_part2_fields

MAGIC_INITIAL = 'ff00f0ca'

# ------------------------------------------------------------------------

# decode bit stream using key and bit format array
# return as dictionary

def ingest( stream:bitstream, decoder:tuple ) ->dict:

  header = {}
  for name, fmt in decoder:
    value = stream.read( fmt )
    if name: header[name] = value
  return header


# Parse MSSS Mini-Header

def read_header( stream:bitstream ) ->dict:

  # parse bit stream as image
  image_header = ingest( stream, image_product_fields )

  # re-parse bit stream as zstack
  stream.pos = 0 # start again to re-parse
  zstack_header = ingest( stream, zstack_product_fields )

  # simple, mostly-reliable, heuristic to determine product type (image or zstack)
  # use of Telemetry application_id in xml label is recommended
  if zstack_header['focus_merge_start_camera_product_id'] == 0: return image_header
  if zstack_header['focus_merge_depth'] in range(1,9): return zstack_header

  return image_header


def decode_huffman( width:int, height:int, stream:bitstream ) ->dict:

  # sync codes precede each 2 line group (a frame)
  SYNC_CODE          = 0xFFFF0000
  LINES_PER_FRAME    = 2
  FRAMES_PER_IMAGE   = height//LINES_PER_FRAME
  BYTES_PER_FRAME    = LINES_PER_FRAME*width


  def get_symbol():
    TRAVERSE_LEFT_BRANCH  = 0b10
    TRAVERSE_RIGHT_BRANCH = 0b01
    node = 0

    for n in range(len(stream)-stream.pos):
      bit = stream.read('uint:1')
      if bit == 0:
        if huffman.traverse[node] & TRAVERSE_RIGHT_BRANCH:
          node = huffman.left[node]   # branch
        else: 
          return huffman.left[node]   # leaf
      else:
        if huffman.traverse[node] & TRAVERSE_LEFT_BRANCH:
          node = huffman.right[node]  # branch
        else: 
          return huffman.right[node]  # leaf

    raise ValueError( 'stream consumed prematurely' )

  frame_buf = bytearray(BYTES_PER_FRAME)

  for frame in range(FRAMES_PER_IMAGE):

      # one sync per frame
      value = stream.read('uint:32')
      if value != SYNC_CODE: 
        raise ValueError( f'Unexpected value 0x{value:x} found instead of sync code 0x{SYNC_CODE:x}' )

      pixel_value = 0

      for n in range(BYTES_PER_FRAME):
         
        value = get_symbol()
        pixel_value = (pixel_value + value) & 0xff   # difference coding
        frame_buf[n] = pixel_value

      yield frame_buf

      # advance to next whole byte boundary
      stream.pos = (stream.pos+7) & ~7

      # advance to next word boundary
      stream.pos = (stream.pos+31) & ~31



# Read image payload after MSSS Mini-Header

def read_image( stream:bitstream, header:dict ):

  cmode   = header['image_compression_mode'] 
  quality = header['image_quality']
  width   = header['samples'] * 8
  height  = header['lines'] * 8
  lut     = header['companding_table']

  print( f'Dimensions {width}x{height}' )

  if lut == 0xff:

    # ----------------------------------------------
    # RAW, 16 bit uncompressed image SIS 4.4.5.1

    print( 'RAW 16' )
    print( 'Writing RAW bytes to img_16bit.raw' )
    open( 'img_16bit.raw','wb').write( stream.tobytes() )

  elif cmode == 0 and quality == 0:

    # ----------------------------------------------
    # RAW, 8 bit uncompressed image SIS 4.4.5.2
    # ----------------------------------------------

    print( 'Writing RAW8 bytes to img_8bit.raw' )
    open( 'img_8bit.raw','wb').write( stream.tobytes() )

  elif cmode == 0xff:

    # ----------------------------------------------
    # Losslessly compressed image, SIS 4.4.5.3 
    # ----------------------------------------------

    print( f'Writing lossless to img_8bit.raw' )

    # output file will contain non-interleaved line group ordering
    with open('img_8bit.raw','wb') as outfile:
      for frame in decode_huffman( width, height, stream ):
        outfile.write( frame )


  elif cmode in ( 0, 1, 2 ):

    # ----------------------------------------------
    # JPEG compression, SIS 4.4.5.4
    # ----------------------------------------------

    SOI = '0xFFD8' # JPEG SOI marker
    offsets = stream.findall( SOI, bytealigned=True) # start of each JPEG in bits
    offsets = list(offsets) + [len(stream)]          # include final bit at EOI
    segments = zip(offsets,offsets[1:])              # as segment pairs

    for frame_num, seg in enumerate(segments):

      out = f'img_{frame_num}.jpg'
      print( f'Writing JPEG frame {frame_num} to {out}' )
      open( out, 'wb').write( stream[seg[0]:seg[1]].tobytes() )

  else:
    print( f'unknown encoding, lut={lut} cmode={cmode}', file=sys.stderr )


def read_file( testfile:str ):

  with open( testfile, 'rb' ) as infile:

    # read file stream as bits
    stream = bitstream( infile )

    # header first
    header = read_header( stream )
    for k, v in header.items(): print( f'{k:35} {v}' )

    # valid?
    if header['magic_initial'] != MAGIC_INITIAL:
      raise ValueError( f'magic_initial is {header["magic_initial"]} not {MAGIC_INITIAL}' )

    # then subsequent image payload, always at byte 64
    read_image( stream[64*8:], header ) # shallow copy



def main():

  parser = argparse.ArgumentParser(description='Decode MSSS ODR DAT files')
  parser.add_argument('infile', help='Input file path (DAT)' )

  args = parser.parse_args()

  try:
    read_file( args.infile )

  except FileNotFoundError:
    print(f'Error: Input file "{args.infile}" not found')
    return
  except Exception as e:
    print(f'An error occurred while reading the input file: {e}')
    return

if __name__ == '__main__':
    main()
