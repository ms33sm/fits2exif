### **Table of Contents** *generated with [md-toc 9.0.0](https://pypi.org/project/md-toc/)*
<!-- md_toc -s 2  --in-place github README_doc.md -->

<!--TOC-->

- [Fits2Exif Documentation](#fits2exif-documentation)
  - [Arguments](#arguments)
  - [Custom tags pairs -- args file](#custom-tags-pairs----args-file)
    - [fits2exif.args file](#fits2exifargs-file)
  - [GPS coordinates](#gps-coordinates)

<!--TOC-->

# Fits2Exif Documentation

Fits2Exif script use python's astropy package to access Fits header(s) keys and exiftool executable to write XMP site cards. It is possible, to some extent, modify script's behavior by passing arguments to it on command line. 

## Arguments


| Switch  |   Option      | Example              | Comment |  
| ----- | --------------  | -------              |-------- |
| -h    | -\-help         | -h                   |Show help|
| -i    | -\- input_file  | -i file1 file2       | Operation will be performed only on passed file(s)
| -e    | -\-exiftool_path| -e /path/to/exiftool | Path to exiftool executable |
| -c    | -\-config_file  | -c /path/to/file.args| Path to args file|
| -FITcol|-\-config_fit_col| -FITcol 0           | Fits key column in args file|
| -XMPcol| -\-config_xmp_col| -XMPcol 1          | Exif tag column in args file|
| -n     | -\-config_xmp_col| -n                 | Disable writing XMP file|
| -wh   | -\-write_header   | -wh                | Write Fits header to -XMP-dc:Description tag in XMP file |
|-fits  | -\-fits_extensions| -fits fITS FiTs    | Change Fits file extension in case they have nonstandard extension like: *.fITS|
| -f | -\-root_folder | -f /path/to/files/folder | Path to folders where FITS files are located |
| -p | -\-files_pattern | -p \*stacked\* | Only files containing patter "stacked" will be selected
| -r | -\-recursive | -r | Will search for all sub-folders in root directory set by -f, or current working directory if no root directory is specified.
| -t | -\-thumbnails | -t | Make thumbnails files of extension set by -te, or by default jpg.
| -te | -\-thumbnail_extension | -te tiff | Make thumbnails files with specified extension. Default is jpg. 
| -w | -\-overwrite | -w | Overwrite thumbnail files. Useful to correct stretch of FITS files (all, or passed in -i switch)
| -u | -\-flip_image | -u | Flip thumbnail files (up side down). |
| -m | -\-mirror_image | -n | Mirror thumbnail files (left to right). |
| -rgb | -\-rgb_image | -rgb | Make thumbnail files in RGB format, despite of FITS images to be not colour.
| -lp | -\-file_config_keys_pairs | -lp tag_pairs_file.args | List pairs of FITS header key and corresponding Exif tags for information and exits. No files written. Default: fits2xmp.args|
| -lf | -\-list_files_only | -lf 1 | Only lists files that will be processed. Default: 0, files will be processed. |
| -s | -\-separation_string | -s "<" |String/character separating EXIF tags and FITS header keys in config file (fits2xmp.args). |
| -k | -\-comment_string    | -k "#" | String/character used mark comments in EXIF tags and FITS header keys in config file (fits2xmp.args). |
| -I | -\-Information | -I | Show configuration and exit (no action made). Use for command prompt arguments testing. |
| -R | -\-resize_picture | -R 480 560 | Set width and height of thumbnail. Default: same size as original FITS file. |
| -F | -\-resize_filter | -F BICUBIC | Filter used for resizing thumbnail. Default: LANCZOS |
| -S | -\-stretch_function | -S HistEqStretch | Stretch function used make stretched thumbnails. For functions required parameters, set them with -P. Default: ManualSigma33 |
| -P | -\-stretch_parameters | -P 1.2 1.2 | Stretch function parameters. 

**Table 1:** Table of arguments.

"Fits key" - "Exif tag" pairs are defined in args file (by default: fits2xmp.args). To add custom pairs, please edit this default file or make your own file and pass it to fits2exif.py by passing "-c" or "-\-config_file" option.
```
python3 fits2exif.py -c my_own_pair.args
```

or equivalently

```
python3 fits2exif.py --config_file my_own_pair.args
```

## Custom tags pairs -- args file

If you have own FITS/EXIF tag pair that you would like to incorporate, this is the place to start. 

It is possible define custom FITS/EXIF tag pairs in args file. Structure of the file is as follows (if defaults settings are used):
* the first column is the EXIF tag <!-- (change to second column by passing -FITcol 1 to script) -->
* next is the separation character(s) "<" <!-- (change to first by passing -XMPcol 0) -->
* third is FITS header key

Any line starting with "#" (or character/string set by -k "character" switch) is treated as comment.

Default separation character is "<". If you use own separation character in args file, inform script about it by passing -s "own_separatio_character" as script argument.

For completeness, it is possible change order of pair tags by setting switch -XMPcol 1 and -FITcol 1. Currently, both columns should be set.

First place to get motivation for custom FITS/EXIF tags pair is your FITS file header. There are [standard FITS keywords](https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html) but you can find many [others](https://fits.gsfc.nasa.gov/fits_dictionary.html), or even unique, in your FITS file. The biggest problem is to find corresponding alternative as [EXIFtool](https://exiftool.org/TagNames/index.html) tag that will be properly interpreted by photo management software. 

If you have your (standard or custom) FITS key and corresponding EXIFtool tag that are not listed in this repository's args file(s), please share them here to help other people make their life easier.

### fits2exif.args file


| FITS header key | EXIFtool tag |
| :--------------- | -------- |
| Focallen | -FocalLength |
| Gain     | -ISO |
| Instrume | -Model    |
| Isospeed | -ISO |
| **XMP-exif tags** | |
| Comment  | -XMP-exif:UserComment |
| Dec      | -XMP-exif:GPSLatitude |
| Ra       | -XMP-exif:GPSLongitude|
| Naxis1   | -XMP-exif:ExifImageWidth |
| Naxis2   | -XMP-exif:ExifImageHeight|
| Bitpix   | -XMP-tiff:BitsPerSample |
| Object   | -XMP-pmi:ObjectType |
| Observer | -XMP-xmp:Author |
| History  | -XMP-photoshop:History |

**Table 2:** Table of FITS/EXIF pairs used in default args file.

## GPS coordinates

All header keywords are copied as is to corresponding tag, except GPS coordinates. To keep view as we see it on the sky, the coordinates are shifted as follows:
$$
longitude = longitude - 180
$$