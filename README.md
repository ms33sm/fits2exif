# Fits2Exif / Fits2XMP

Python script  allowing photo management software to access  information in astronomy Fits files by converting selected Fits' header keys to "semi-equivalent" Exif's tags saved as XMP site cards.

## 1. Introduction

Fits files are typically used in astrophotography to store pictures and various information (like: used equipment; photographed object, its coordinates etc.) Access and effectively manage these information specialized software is necessary.

On other hand, there are many photography programs capable of tagging, selecting or editing photos. However, these programs are missing access to Fits information, or missing the information "understanding".

Purpose of this python script is to bridge gap between "Fits language" and photo management software. "Translating" Fits information to photography language allow amateur astrophotographers unleash powerful software features to astronomical photos.

## 2. Features

* Convert Fits' header keys to Exif tags
* Fits header key - Exif tag pair is user selectable
* Exif tags are written in XMP site card 
* Produce thumbnails of Fits file (if desired)
* Thumbnails can be stretch by selectable stretch functions
* Run on all major platforms (Linux, Mac, Windows), where exiftool and python can be installed.

## 3. Getting started

### Dependencies
	
This script requires two programs to be present (not necessary installed) on operating system.
* python 
* [Exiftool](https://exiftool.org/install.html) executable
* Optional: [pyexiftool](https://github.com/sylikc/pyexiftool) speeds writing XMP files significantly. 
 
It is possible to install most important python dependences to user's directory by executing following command:
```bash
pip3 install pyexiftool pillow astropy --user
```
Each distribution provides their own packages. Install them from your distro repositories, if you know how, or by executing above command as normal user. It is your choice. 

### Installation

No installation is required to produce XMP site cards. 

Python file from this repository is sufficient, however, if you want to practice with test Fits files and get default Fits - Exif file, download the ZIP file (clicking on right top icon: <>Code)
```bash
cd to\downloaded\zip_file\folder
unzip fits2exif-main.zip
```
or clone the repository
```
git clone git.XYZ
```
Go to newly created folder (unzip-ed or clone-ed)
```
cd fits2exif-main
```

### Usage
Produce XMP site cards for all Fits file in the same folder as python script, run following command

```
python3 fits2xmp.py
```
If you need to make also stretched jpg thumbnails files for all Fits files in directory, along the XMP files, run
```
python3 fits2xmp.py -t 
```

For additional usage like:
-  define custom pairs of Fits header key - Exif tag
-  pass individual Fits files to process
-  set pattern for Fits files selection
-  set path to exiftool executable
-  set "root" directory

refer to [documentation](./documentation/README_doc.md) file.

## License

## Contribution

<!---
## Items to introduce
* FITS2XMP / Fits2Exif
* Description (what it does)
* Motivation (?) (Maintain Fits <-> Exif equivalents)
* Features
	- Run on all major platforms (Linux, Mac, Windows), where python and exiftool can be installed.
	- 
* Relation to another project? (may be in Motivation)
* Getting started
	- Prerequisites
	- Dependencies
	- Installation
	- How to run (Documentation; link to different page?; Include screenshot of runig program/results. Include Examples?)
* License (?)
* Contribution (?)
* [Link](http://localhost.com "Title")
* Make photo / logo 
-->