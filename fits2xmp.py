#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:39:39 2024

@author: ms33sm
"""

"""
TODO
    done make switch for WxH of output thumbnail
    done finish resize_filter - set it from command line
    done make switch for STRETCH functions
    !!! Add to list keys from ALL headers in get_FITS_data() function
    !!! If set -FITcol or -XMPcol in commandpromt, make sure that another is set automatically as well
    !!! Make Class ?
    !!! debug mode ?
    !!! add test files
    
"""
#### imports 
import os
import sys
import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path

from PIL import Image, ImageOps

from astropy.io import fits
from astropy.visualization import (AsinhStretch, ContrastBiasStretch,
                                   HistEqStretch,
                                   LinearStretch,
                                   LuptonAsinhStretch,
                                   PowerDistStretch,
                                   PowerStretch,
                                   SinhStretch,
                                   ManualInterval,
                                   SqrtStretch,
                                   MinMaxInterval,
                                   # SymmetricInterval,
                                   # LogStretch
                 )

# try:
#     ## another stretch from : https://github.com/LCOGT/auto_stretch
#     from auto_stretch import apply_stretch
#     auto_stretch = True
# except:
#     auto_stretch = False

try:
    import exiftool
    pyexiftool_present = True
except:
    pyexiftool_present = False
    print("Using exiftool executable directly. Slow!")
    print("Consider to install pyexiftool from: https://github.com/smarnach/pyexiftool/tree/master")
    print("     pip install pyexiftool --user")

global exit_program
exit_program = False

resize_image = False
resize_filter = None
separation_string = "<"
comment_string = "#"

def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"{string} not found.")

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"Path: {string} is not valid folder path.")

def set_filter(use_filter):
    match use_filter:
        case 'NEAREST':
            resize_filter = Image.NEAREST
        case 'BILINEAR':
            resize_filter = Image.BILINEAR
        case 'HAMMING':
            resize_filter = Image.HAMMING
        case 'BICUBIC':
            resize_filter = Image.BICUBIC
        case 'LANCZOS':
            resize_filter = Image.LANCZOS
        case 'None':
            resize_filter = None
    return resize_filter

def conversion_strings_fromFits2XMP(configfile, comment_char='#', split_char=separation_string):
    """
    Read configfile for future conversion between FITS keywors in header
    and xmp keywords.
        
    Parameters
    ----------
    configfile : textfile
        Example of config file:
            # FITS header = XMP equivalent
            # keywords should be separated by character "split_char='<'" // without quotes
            #   FITS keywords are not case sensitive: BItpix is the same as bitpix
            PIXBIT < -XMP-tiff:BitsPerSample
            History < -XMP-photoshop:History
            
    comment_char : char, optional
        DESCRIPTION. The default is '#'.

    Returns
    -------
    fromFITStoXMP : list
        Pair of equivalent FITS and XMP keywords.
    """
    fromFITStoXMP = []
    with open(configfile) as file:
        for line in file:
            
            line = line.strip()
            # print(line)
            line = line.split(comment_char)[0] # remove comment
            if not (line.startswith((comment_char,' ')) or line == ''):
                line = line.strip().strip()
                line = line.split(split_char)
                fromFITStoXMP.append([line[0].strip(),line[1].strip()]) #storing everything in memory!
    return fromFITStoXMP

def show_keys_pairs(filelist):
    global exit_program
    if filelist is None:
        filelist="fits2xmp.args"
    
    print("Listing pairs for configuration file: {}".format(filelist))
    pairs = conversion_strings_fromFits2XMP(filelist)
    # print("pairs:{}",pairs)
    for pair in pairs:
        print(pair)
    exit_program = True

def manual():
    # print("Maual how to use this python script.")
    print()
    print("Make XMP files from FITS files in current working directory, run:")
    print("     python fits2xmp.py")
    print()
    print("Script looks for configuration file to convert FITS tags to XMP tags.")
    print("     For structure of this file, see provided file: fits2xmp.args")
    print("     Add desired pairs of tags in form: FITS_tag < XMP_tag")
    print()

### initialize parser and parse input parameters
parser = argparse.ArgumentParser(
        prog='fits2xmp',
        description='Converts FITS file\'s header to XMP sidecar. Optionally, generates thumbnail-like file of FITS\'s file.',
        usage=manual(),
     )

parser.add_argument("-i", "--input_file",   nargs='+', help="Input FITS file(s). One or separated by space.", default=[])
parser.add_argument("-e", "--exiftool_path", type=str, help="Path to exiftool executable.", default="exiftool")
parser.add_argument("-c", "--config_file",   type=file_path, help="Path to fits2xmp.args configuration file.", default="./fits2xmp.args")
parser.add_argument("-FITcol", "--config_fit_col", type=int, choices=[0, 1], help="Column of FITS tag in config file.", default=1)
parser.add_argument("-XMPcol", "--config_xmp_col", type=int, choices=[0, 1], help="Column of XMP tag in config file.", default=0)

parser.add_argument("-wh", "--write_header", help="Write whole FITS\'s header to -XMP-dc:Description", action="store_true")

parser.add_argument("-n",  "--no_xmp",       help="Disable xmp file output. XMP files will not be written.", action="store_false")
parser.add_argument("-fits", "--fits_extensions", action='append', default=['fit', 'fits', 'FIT', 'FITS'])

parser.add_argument("-f", "--root_folder",   type=dir_path, help="Root folder of FITS files.", default=os.getcwd())
parser.add_argument("-p", "--files_pattern", type=str, help="Pattern for file filtering.\n Example: * - select all files; *stacked* - select only files with name containing word stacked", default='*')

parser.add_argument("-r", "--recursive",  help="Performs recursive search in all folders in cwd or folder defined by --root_folder.", action="store_true")

parser.add_argument("-t", "--thumbnails", help="Make thumbnails-like file for each found FITS file.", action="store_true")
parser.add_argument("-te","--thumbnails_extension", type=str, help="Make thumbnails file of chosen type. Default: xyz.jpg", default="jpg")
parser.add_argument("-w", "--overwrite",  help="Overwrite thumbnail file.", action="store_true")

parser.add_argument("-u", "--flip_image",   help="Flip thumbnail file.",    action="store_false")
parser.add_argument("-m", "--mirror_image", help="Mirror thumbnail file.",  action="store_true")
parser.add_argument("-rgb", "--rgb_image",  help="Save color thumbnail even FITS file is only black and white.", action="store_true")

parser.add_argument("-lp", "--file_config_keys_pairs", nargs='+', type=show_keys_pairs, help="Show equivalent key pair for config file and exit. Default: fits2xmp.args")
parser.add_argument("-lf", "--list_files_only", type=int, choices=[0, 1], help="Only lists files that will be processed", default=0)       

parser.add_argument("-s", "--separation_string", type=str, help="String/character separating EXIF tags and FITS header keys in config file (fits2xmp.args).", default="<")
parser.add_argument("-k", "--comment_string", type=str, help="String/character used mark comments in EXIF tags and FITS header keys in config file (fits2xmp.args).", default="#")

parser.add_argument("-I", "--Information", help="Show configuration and exit (no action made). Use for command prompt arguments testing.", action="store_true")


parser.add_argument("-R", "--resize_picture", type=int, nargs=2, help="Resize thumbnail-like picture to passed int numbers: witdh height. By default thumbnail will have size of original FITS file.")
parser.add_argument("-F", "--resize_filter", type=str, choices=['NEAREST', 'BILINEAR', 'HAMMING', 'BICUBIC', 'LANCZOS', 'None'], default='LANCZOS', help="Filter for resized thumbnail-like picturet. Default: LANCZOS.")
parser.add_argument("-S", "--stretch_function", type=str, choices=['AsinhStretch',
                                                                   'ContrastBiasStretch', 
                                                                   'HistEqStretch',
                                                                   'LinearStretch', 
                                                                   # 'LogStretch', 
                                                                   'LuptonAsinhStretch',
                                                                   'PowerDistStretch',
                                                                   'PowerStretch',
                                                                   'SinhStretch',
                                                                   'SqrtStretch',
                                                                   'MinMaxInterval',
                                                                   'ManualInterval',
                                                                   'ManualSigma33',
                                                                   ], default='ManualSigma33', help="Stretch function from astropy.vizualization. Default: ManualSigma33.")
parser.add_argument("-P", "--stretch_parameters", type=float, nargs=2, help="Parameters for stretch functions that need them. Ignored otherwise. Default: [1, 2]")

args = parser.parse_args()
    
################### Set parsed values #################
files=args.input_file
exiftool_executable=args.exiftool_path #"exiftool" # If nothing specified it is assumed that exiftool executable is in $PATH
configfile=args.config_file #'fits2xmp.args' # Config file where FITS to XMP pairs are set
make_xmp = args.no_xmp #True # If True, the XMP files are made
config_fits_col=args.config_fit_col # column of fits tags in configfile
config_exif_col=args.config_xmp_col #0 # column of xmp  tags in configfile
add_header = args.write_header #True # if True, add whole fits header to: -XMP-dc:Description tag

### !!! potentially set following 2 variables from command line
compatible_with_commercial_programs=True
overwrite_original=True # for exiftool to not make xmp file with new extension: xyz.xmp_original

fits_extensions=args.fits_extensions #['fit', 'fits', 'FIT', 'FITS']
root_folder = args.root_folder #os.getcwd()
find_file_pattern=args.files_pattern #'*' # default: "*"; what files want to filter; "*" all files with fits_extensions e.g: *stacked* find files containting word: stacked
find_files_recursively=args.recursive #False

force_folder_pattern = '' # "all in one": folder pattern for glob.glob(): e.g: "/root/home/dir/astro/**/*.fits" 

make_thumbnails = args.thumbnails #False
thumbnails_extension = args.thumbnails_extension #".jpg"
thumbnail_overwrite = args.overwrite #True
flip_image_h = args.flip_image #True # default: True (keep same orientation as fits file); flip final image horizontally
mirror_image_v = args.mirror_image #False # default: False; mirror final image
save_image_as_RGB = args.rgb_image #False # default: False; if True convert final image to RGB even it is gray image
stretch_function = args.stretch_function #"ManualInterval"

separation_string=args.separation_string
comment_string=args.comment_string

# xxx = args.file_config_keys_pairs
list_files_only = args.list_files_only

if args.stretch_parameters is not None:
    stretch_parameters = args.stretch_parameters
else:
    stretch_parameters = [1, 2]

information=args.Information

if args.resize_picture is not None:
    resize_image = True
    thumbnails_width = args.resize_picture[0]
    thumbnails_height = args.resize_picture[1]
    
if args.resize_filter is not None:
    resize_filter = set_filter(args.resize_filter)

def print_parsed_arguments():
    print("Command prompt parsed arguments:")
    print(f"  path to exiftool executable        : {args.exiftool_path}")
    print( "                                         change with: -e /path/to/exiftool)")
    print(f"  path to configuration file         : {args.config_file}")
    print( "                                         (change with: -c /path/to/configfile)")
    print(f"  write xmp file                     : {args.no_xmp}\t(change with -n)")
    
    print(f"  extensions of FITS files           : {args.fits_extensions}")
    print( "                                         (add fItS with: -fits fItS)")
    
    print(f"  root folder for FITS files         : {args.root_folder}")
    print( "                                         (change with: -f /path/to/fits/folders)")
    print(f"  write whole FITS header to XMP file : {args.write_header} \t(change with: -wh)")
    print(f"  FITS tag column in config file      : {args.config_fit_col} \t(change to Y with: -FITcol Y)")
    print(f"  XMP  tag column in config file      : {args.config_xmp_col} \t(change to P with: -XMPcol P)")
    print(f"  pattern to select files with word   : {args.files_pattern}\t(change to *NEW* with: -p *NEW*)")
    print(f"  look for FITS files recursively     : {args.recursive}\t(change with: -r)")
    print(f"  make thumbnail-like file for FITS   : {args.thumbnails}\t(change with: -t)")
    print(f"  extension of thumbnails-like file   : {args.thumbnails_extension}\t(change to tiff with: -te tiff)")
    print(f"  overwrite tumbnails-like files      : {args.overwrite}\t(change with: -w)")
    print(f"  flip thumbnail-like file            : {args.flip_image}\t(change with: -u)")
    print(f"  mirror thumbnail-like file          : {args.mirror_image}\t(change with: -m)")
    print(f"  make RGB thumbnails-like file       : {args.rgb_image}\t(change with: -rgb)")
    print(f"  process individual files            : {args.input_file}\t(add files with: -i file1 file2)")
    print(f"  show parsed arguments info          : {args.Information}\t(show with: -I)")
    
    print(f"  list FITS / XMP pair keys           : {args.file_config_keys_pairs}\t(show with -lp xyzfile.args)")
    print(f"  only list all files to process      : {args.list_files_only}\t(show all files to process with -lf 1 )")
    print(f"  string separating tags in args file : {args.separation_string}\t(change to : with -s ':')")
    print(f"  string used as comments in args file: {args.comment_string}\t(change to : with -k ':')")
    
    print(f"  width and height of thumbnail file  : {args.resize_picture}\t(change to : with -R intWidth intHeight)")
    print(f"  filter for thumbnail file resizing  : {args.resize_filter}\t(change to : with -F Filter)")
    print(f"  stretch function for thumbnail file : {args.stretch_function}\t(change to : with -S StretchFunction)")
    print(f"  stretch function parameters         : {stretch_parameters}\t(change to : with -P real real)")

    exit()

### !!! TODO: Add to list keys from ALL headers 
def get_FITS_data(fits_file, _index=0):
    """
    Read fits_file and return header, keys and hdul in the first header.

    Parameters
    ----------
    fits_file : file
        FITS file.

    Returns
    -------
    fits_header, fits_keys_list, hudl
        DESCRIPTION.stretch

    """
    # check if file exists    
    # get fits keywords from first header
    try:
        hdul = fits.open(fits_file)
        hdr = hdul[_index].header
        fits_keys = list(hdr.keys())
        return hdr, fits_keys, hdul
    except:
        return None, None, None

### some stretches available in astropy.vizualization
def stretch(data, stretch_function="ManualInterval", sigma=[1,2]):       
    match stretch_function:

        case "AsinhStretch":
            # print("AsinhStretch")
            stretch = AsinhStretch()
            return stretch(data, sigma[0])
        
        case "ContrastBiasStretch":
            # print("AsinhStretch")
            stretch = ContrastBiasStretch(sigma[0], sigma[1])
            return stretch(data)

        case "HistEqStretch":
            # print("HistEqStretch")
            stretch = HistEqStretch(np.array(data))
            return stretch(data)
        
        case "LinearStretch":
            # print("LinearStretch")
            stretch = LinearStretch(sigma[0], sigma[1])
            return stretch(data)  
        
        # case "LogStretch":
        #     # print("LogStretch")
        #     stretch = LogStretch(sigma[0])
        #     return stretch(data)  
        
        case "LuptonAsinhStretch":
            # print("LogStretch")
            stretch = LuptonAsinhStretch(sigma[0], sigma[1])
            return stretch(data)
        
        case "PowerDistStretch":
            # print("LogStretch")
            stretch = PowerDistStretch(sigma[0])
            return stretch(data)
        
        case "PowerStretch":
            # print("LogStretch")
            stretch = PowerStretch(sigma[0])
            return stretch(data)
        
        case "SinhStretch":
            # print("SinhStretch")
            stretch = SinhStretch(sigma[0])
            return stretch(data)
        
        case "SqrtStretch":
            # print("sqrt")
            stretch = SqrtStretch()
            return stretch(data) 
        
        case "MinMaxInterval":
            # print("sqrt")
            stretch = MinMaxInterval()
            return stretch(data) 
        
        case "ManualInterval":
            # print("sqrt")
            stretch = ManualInterval(sigma[0], sigma[1])
            return stretch(data) 
        
        case "ManualSigma33":
            # print("ManualInterval")
            average=np.average(data)
            stdev=np.std(data)           
            stretch = ManualInterval(vmin=average-sigma[0]*stdev, vmax=average+sigma[1]*stdev)
            return stretch(data)          
        
        # case "SymmetricInterval":
        #     # print("sqrt")
        #     stretch = ManualInterval(sigma[0], sigma[1])
        #     return stretch(data) 
        
        # case "AutoStretch":
        #     if auto_stretch:
        #         return apply_stretch(data)
        #     else:
        #         print("No stretch applied. Missing auto_stretch package.")
        #         return data
            
        # If an exact match is not confirmed, this last case will be used if provided
        case _:
            # print("Default stretch")
            average=np.average(data)
            stdev=np.std(data)           
            stretch = ManualInterval(vmin=average-sigma[0]*stdev, vmax=average+sigma[1]*stdev)
            return stretch(data)
        
def list_files(root_folder, find_file_pattern, find_files_recursively):
    if (force_folder_pattern == ''):
        if root_folder == '':
            root_folder = "."
        tmp_find_path = root_folder
        if find_files_recursively:
            tmp_find_path +="/**/" # needs to be here for glob.glob() recursive search 
        else:
            tmp_find_path +="/"
        tmp_find_path += find_file_pattern
        tmp_find_path += "."
        for fit_ext in fits_extensions:
            tmp_path=tmp_find_path+fit_ext
            # print(tmp_path)
            files.extend( glob.glob(tmp_path, recursive=True) )
    else:
        files.extend( glob.glob(force_folder_pattern, recursive=True) )
        
    return files

###############################################################################

### exit program if required befor by seting "show_keys_pairs" etc.
if exit_program:
    sys.exit(0)
    
if information:
    print_parsed_arguments()
    
convpairs=conversion_strings_fromFits2XMP(configfile) # conversion pairs
  
print("Collecting files ...")

if files == []:
    files = list_files(root_folder, find_file_pattern, find_files_recursively)

if list_files_only:
    for f in files:
        print(f)
    print("Listed {} files to process.".format(len(files)))
    sys.exit(0)
   
current_file_count=1
total_files=len(files)
print("Found {} files to process.".format(total_files))

print("Processing ...")

### Start pyexiftool if present in system
et = None
if pyexiftool_present:
    et = exiftool.ExifToolHelper()
    et.run()

### Start FITS to XMP and thumbnails generation if desired
fits_keys = None
for fit_file in files:
    if not os.path.isfile(fit_file):
        print(f"File: {fit_file} was not found. Skipping it ...")
        break
    
    print('{} of {}'.format(current_file_count, total_files) )
    current_file_count += 1
        
    hdr,fits_keys, hudl = get_FITS_data(fit_file)

    if fits_keys is None:
        continue
    
    xmp_fields_and_values=[] 
    for fitsitem in fits_keys:
        if fitsitem != '':
            for convitem in convpairs:
                if convitem[config_fits_col].upper() == fitsitem.upper():
                    xmp_fields_and_values.append([fitsitem, convitem[config_exif_col]])
    
    filename, file_extension = os.path.splitext(fit_file)

    if make_thumbnails:
        try:
            if thumbnails_extension[0] != '.':
                # add '.' at beginning
                thumbnails_extension = thumbnails_extension[:0] + '.' + thumbnails_extension[0:]
            
            thumbnail_name=filename+thumbnails_extension
            
            if (not os.path.isfile(thumbnail_name) or (thumbnail_overwrite) ):

                image_data = hudl[0].data
                thumbnail= np.zeros(image_data.shape, 'uint8')
                
                print(stretch_function)

                if len(image_data.shape) == 2:
                    thumbnail = np.round(255*stretch(image_data/np.max(image_data), stretch_function=stretch_function,sigma=stretch_parameters) ).astype(np.uint8)
                else:
                    for i in range(len(image_data)):
                        thumbnail[i,:,:]=np.round(255*stretch(image_data[i]/np.max(image_data[i]), stretch_function=stretch_function)).astype(np.uint8)

                    thumbnail = np.moveaxis(thumbnail, 0, -1) # Image.fromarray needs array like: [x,y,RGB]

                im = Image.fromarray( thumbnail )
                
                if resize_image:
                    im = im.resize( (thumbnails_width, thumbnails_height), resample=resize_filter )
                    print(resize_filter)

                if save_image_as_RGB:
                    im = im.convert('RGB')

                if flip_image_h:
                    im = ImageOps.flip(im)

                if mirror_image_v:
                    im = ImageOps.mirror(im)
                    
                print(f'Saving file {thumbnail_name}')
                im.save(thumbnail_name)
        except:
            pass
                    
    if compatible_with_commercial_programs:
        fname = filename
    else:
        fname = filename + file_extension

    if make_xmp:
        fname = fname + '.xmp'
        Path(fname).touch()
        for tag_to_add in convpairs:
            if tag_to_add[config_fits_col].upper() in fits_keys:
                tag_from_fits=str(hdr[tag_to_add[config_fits_col]]).strip()
                tag_from_config=tag_to_add[config_exif_col].strip()
                
                # if GPS coordinates are present and need to be written, they are converted to "sky" coordinates
                if tag_to_add[config_exif_col].strip() == '-XMP-exif:GPSLongitude':
                    longitude = float(hdr[tag_to_add[config_fits_col]])
                    print("Long before: ", longitude)
                    if longitude > 180.0:
                            longitude = longitude - 360 
                    print("Long after: ", longitude)
                    tag_from_fits = str(-1*longitude)
                    
                tag=tag_from_config + '=' + tag_from_fits
                tag = tag.replace("'","").replace("\r\n","").replace("\n"," ")
                
                if overwrite_original:
                    if pyexiftool_present:
                        et.execute( "-q" ,"-overwrite_original", tag , fname)
                    else:
                        subprocess.call([exiftool_executable, "-q","-overwrite_original", tag, fname ])
                else:
                    if pyexiftool_present:
                        # cmd = tag + fname
                        et.execute(tag.encode("utf-8") , fname,  raw_bytes=True)
                    else:
                        subprocess.call([exiftool_executable, tag, fname ])
            
        if add_header:
            header_command="-XMP-dc:Description"+"="+'"'+str(hdr)+'"'
            if overwrite_original:
                try:
                    if pyexiftool_present:
                        et.execute( "-q" ,"-overwrite_original", header_command , fname)
                    else:
                        subprocess.call([exiftool_executable, "-q","-overwrite_original", header_command, fname ])
                except:
                    pass
            else:
                try:
                    if pyexiftool_present:
                        et.execute(header_command , fname)
                    else:
                        subprocess.call([exiftool_executable, header_command, fname ])
                except:
                    pass

if pyexiftool_present:
    et.stop()