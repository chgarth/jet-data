#! /usr/bin/env python3

import sys
import nrrd
import vtk
from vtk.util.numpy_support import *
import numpy as np

STEP = int(sys.argv[1])
TIME = float(sys.argv[2])
REFINE = int(sys.argv[3])
SPLIT = 1
VSIZE = 2**REFINE
GMINS = (-0.5, -0.5, -0.5)


# read all blocks into one large array
data = np.empty( (VSIZE, 2*VSIZE, VSIZE, 3), dtype=np.float32 )
data[:] = np.nan

for n in range(0,2*8**SPLIT):
    blk, header = nrrd.read( f'tmp_{STEP:04d}_{n:04d}.nrrd' )

    mins = header['axis mins']
    maxs = header['axis maxs']

    imins = [ int((m-o)*VSIZE) for m, o in zip(mins, GMINS) ]
    imaxs = [ int((m-o)*VSIZE) for m, o in zip(maxs, GMINS) ]
    isize = [a-b for a,b in zip(imaxs,imins)]

    data[ imins[0]:imaxs[0], imins[1]:imaxs[1], imins[2]:imaxs[2], : ] = np.reshape( blk, (isize[0], isize[1], isize[2], 3) )

# save array as VTI file

array = numpy_to_vtk(data.ravel(), array_type=vtk.VTK_FLOAT)
array.SetNumberOfComponents(3)
array.SetName("velocity")

time = vtk.vtkFloatArray()
time.InsertNextTuple( [TIME] )
time.SetName("TimeValue")

vti = vtk.vtkImageData()
vti.SetOrigin( -0.5, -0.5, -0.5 )
vti.SetDimensions( VSIZE, 2*VSIZE, VSIZE )
vti.SetSpacing( 1.0/(VSIZE-1), 2.0/(2*VSIZE-1), 1.0/(VSIZE-1) )

vti.GetPointData().AddArray(array)
vti.GetFieldData().AddArray(time)

dsw = vtk.vtkXMLImageDataWriter()
dsw.SetFileName( sys.argv[4] )
dsw.SetCompressorTypeToLZMA()
dsw.SetInputData(vti)
dsw.Write()

print( f"converted step {STEP}", file=sys.stderr )