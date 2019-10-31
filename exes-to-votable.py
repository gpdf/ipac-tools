# Simple example of a format-conversion program for SOFIA EXES 1D spectral data.
#
# Initially handles only "mrgordspec"-format data, file code MRD.
#
# Resulted from R&D for how to use Firefly to display some of IRSA's holdings
# that have inadequate internal metadata to permit useful "generic" display.
#
# Copyright (c) Gregory Dubois-Felsmann, Caltech/IPAC, 2019

import sys, os
import argparse

from astropy.io.votable.tree import VOTableFile, Resource, Table, Group, Param, Field
from astropy.io import fits

def make_fits_header_group( hdr, table=None ):
    """
    Take an Astropy FITS HDU object as input and returns an Astropy VOTable
    `Group` object that contains all the raw headers from the HDU as Params.

    Makes no attempt to interpret the headers; the point of this is just to
    get a complete copy.  Usually the resulting <GROUP> will be included by
    the caller in a <TABLE> derived from the HDU's data content, so a Table
    argument is provided for this to allow the parent Table to be available
    at construction time, as required by the Group and Table constructors.
    """
    g = Group( table, 
               ID='FITS-Headers', name='FITS-Headers', 
               ucd='meta.fits' )
    g.description = 'Verbatim copy of input FITS headers'

    for i in range(len(hdr.cards)):
        # * access is integer-indexed so that we capture multiple headers
        #   of the same key correctly, if they occur
        # * as obvious as it may seem, we cannot use the Astropy header
        #   dictionary's value as the Param value, because it is round-tripped
        #   through Astropy's idea of how to do type conversion.  We need to
        #   save the literal original format of the header card.
        card  = str(hdr.cards[i])
        hdrid = f"FITS-{i}"

        # make the <PARAM> for this header
        p = Param( table, ID=hdrid, name=hdrid, value=card, arraysize='*' )

        g.entries.append(p)

    return g


### Here beginneth the main program

# Define arguments: for now, a single input filename, and an "overwrite output" flag.
parser = argparse.ArgumentParser(description='Convert an EXES L3 spectral' + \
                                     ' FITS file to a VOTable for Firefly display')
parser.add_argument('input_fits', type=argparse.FileType('rb'),
                    help='path to file to be processed')
parser.add_argument('--force', '-f', action='store_true',
                    help='overwrite output file (suffix .votable)')

args = parser.parse_args()

# Echo the arguments (for now #DEBUG)
print('file to process is <' + args.input_fits.name + '>')

outputfile = os.path.splitext(os.path.basename(args.input_fits.name))[0] + '.votable'
print('output file is <' + outputfile + '>')

mode = 'wb' if args.force else 'xb'

print('... which can' + ( '' if args.force else ' not' ) + ' be overwritten; mode: ', mode)

# Open the input file.  Perhaps this can be made into a loop over multiple files later.
try:
    hdulist = fits.open( args.input_fits )
except Exception as exc:
    raise RuntimeError('Error reading input file ' + args.input_fits.name + ' as FITS.') from exc

# Right now we only know how to handle a very specific file format from the EXES instrument.
# Not entirely clear how to test for it up front.
if len(hdulist)!=1:
    print(sys.argv[0], 'can only process FITS files with a single, primary HDU.',file=sys.stderr)
    sys.exit(1)

# Process the FITS headers into a <GROUP> in the first <TABLE>.
hdr = hdulist[0].header

# Set up the main structure of the VOTable file:
vt = VOTableFile()

#   Create a single <RESOURCE> in the file
r = Resource()
vt.resources.append(r)

#   Place an empty (for now) <TABLE> in the file
t = Table(vt)
r.tables.append(t)

#   Create the special <GROUP> for the raw FITS headers
g = make_fits_header_group( hdr, t )
t.groups.append(g)


# Access the actual file data.  We can only have a primary HDU - this means it must be an IMAGE HDU.
pixel_data = hdulist[0].data

#   Here is the place for some introspection, if any is possible, into the data format.
#   For now this is hard-coded for one particular EXES file format.

#   Confirm that the data HDU shape is as expected.  This kind of EXES file has four stripes of data:
#     0: wavenumber
#     1: flux per wavenumber bin
#     2: flux error
#     3: reference spectrum
assert pixel_data.shape[0] == 4

#   Confirm that the data are 64-bit floats.
assert hdr['BITPIX'] == -64

fits_pixel_type = ''
if hdr['BITPIX'] == -64:
    fits_pixel_type = 'double'
elif hdr['BITPIX'] == -32:
    fits_pixel_type = 'float'
else:
    raise RuntimeError('Unrecognized BITPIX value for file converstion')

#   Create the FIELD elements
f1 = Field(t, name='wavenumber',
           datatype=fits_pixel_type, ucd='em.wavenumber;em.MIR', unit='cm-1')
f1.description = 'wavenumber merged over orders'

f2 = Field(t, name='intensity',
           datatype=fits_pixel_type, ucd='phot.flux.density;em.MIR', unit='erg.s-1.cm-1.sr-1')

f3 = Field(t, name='intensity_err',
           datatype=fits_pixel_type, ucd='stat.error;phot.flux.density;em.MIR', unit='erg.s-1.cm-1.sr-1')
f3.description = 'error (standard deviation)'

f4 = Field(t, name='transmission_reference',
           datatype=fits_pixel_type, ucd='phys.transmission;em.MIR', unit='')
f4.description = 'reference modeled atmospheric transmission spectrum'

#   Add the fields to the table
t.fields.extend( [f1, f2, f3, f4] )

#   Allocate the space for the four columns in the table
t.create_arrays( pixel_data.shape[1] )

#   Copy over each row of the table.  Is there a better way?
for i in range(pixel_data.shape[1]):
    values = tuple(pixel_data[:,i])
    t.array[i] = values

# Output the XML VOTable file
with open(outputfile,mode) as f:
    vt.to_xml(f)

sys.exit(0)
