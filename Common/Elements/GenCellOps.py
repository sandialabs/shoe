#!/usr/bin/python2.2
#
# Copyright 2012 Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the
# U.S. Government. Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that this Notice and any
# statement of authorship are reproduced on all copies.
#
#
#	Run this script with 4 filenames:
#	1. The list of vtkCellOps member variables
#	2. The file specifying the member values for each cell type
#	3. The name of the output file to hold the structure instantiation code
#	Here's an example:
# python2.2 GenCellOps.py CellOpsMembers.txt CellOpsSpecifications.txt vtkCellOpsImpl.cxx

import re
import sys

#import sqrt for functions in Pad1DArray and Pad2DArray
from math import sqrt

def Pad1DArray( s, nr, pad=-1 ):
  reOpenCurly = re.compile( '{' )
  reCloseCurly = re.compile( '}' )
  stmt = reOpenCurly.sub( '[', reCloseCurly.sub( ']', s ) )
  vararr = eval(stmt) # now we have an array
  vararr[len(vararr):] = [ pad for i in range(nr - len(vararr)) ]
  return re.sub( '\[', '{', re.sub( '\]', '}', str(vararr) ) )

def Pad2DArray( s, nc, nr, pad=-1 ):
  reOpenCurly = re.compile( '{' )
  reCloseCurly = re.compile( '}' )
  stmt = reOpenCurly.sub( '[', reCloseCurly.sub( ']', s ) )
  vararr = eval(stmt) # now we have an array
  vararr[len(vararr):] = [ [] for i in range(nr - len(vararr)) ]
  for row in vararr:
    row[len(row):] = [ pad for i in range(nc - len(row)) ]
  return re.sub( '\[', '{', re.sub( '\]', '}', str(vararr) ) )

def WriteOrderedMembers( s, i, p, d ):
	print >> cellout, "static vtkCellOps %s%s%sOps = {" % ( s, i, p )
	for m in members:
		padded = m.ljust( memberpad )
		if d.has_key(m):
			print >> cellout, '\t/* %s */ %s' % ( padded, d[m] ),
		else:
			print >> cellout, '\t/* %s */ 0' % padded,
		if m == lastmember:
			print >> cellout
		else:
			print >> cellout, ','
	print >> cellout, '};'
	print >> cellout

def WriteMemberInit( s, i, p, d ):
	print >> cellout, '  int %s%s%sClassId = RegisterCellType( (vtkCellType) { %s, %s, %s }, &%s%s%sOps );' % ( s, i, p, s, i, p, s, i, p )
	print >> cellout, '  (void)%s%s%sClassId;' % ( s, i, p );

def SaveMemberValues( t, d ):
	if t[0] == '*':
		sa = shapes
	else:
		sa = [ t[0] ]

	if t[1] == '*':
		ia = interps
	else:
		ia = [ t[1] ]

	if t[2] == '*':
		pa = prodspaces
	else:
		pa = [ t[2] ]

	#print " Over %s, %s, %s " % ( sa, ia,pa )

	for s in sa:
		for i in ia:
			for p in pa:
				for k in d.keys():
					AllCellOps[ (s, i, p) ][ k ] = d[ k ]

shapes = ( 'Point', 'Curve', 'Rod', 'Tube', 'Triangle', 'Quadrilateral', 'TriangleShell', 'QuadrilateralShell', 'Tetrahedron', 'Hexahedron', 'Wedge', 'Pyramid' )
interps = ( 'Legendre', 'Lagrange', 'BSpline' )
prodspaces = ( 'MaxTotalOrder', 'TruncatedTotalOrder', 'Tensor' )

header_file = 'vtkCellOps.h'
if len( sys.argv ) > 1:
	header_file = sys.argv[1]

opspec_file = 'CellOpsSpecifications.txt'
if len( sys.argv ) > 2:
	opspec_file = sys.argv[2]

impl_file = 'vtkCellOpsImpl.cxx'
if len( sys.argv ) > 3:
	impl_file = sys.argv[3]

hdrfile = file( header_file, 'r' )
defines = file( opspec_file, 'r' )
cellout = file( impl_file, 'w' )

header = ''
header_arr = hdrfile.readlines()
for line in header_arr:
	header += line

reComm = re.compile( "//.*\n" )
header = reComm.sub( '', header )

reClss = re.compile( "struct[ \t]+VTK_SNL_COMMON_EXPORT[ \t]+vtkCellOps[ \t\n]*\{(.*\n)*}[ \t]*;" )
match = reClss.search( header )
if match == None:
	print 'Error!!!'
header = header[match.regs[0][0]:match.regs[0][1]]

# Here we make an assumption about the input file...
# All member variables must be declared before any
#	members that are function pointers. Otherwise,
#	the two searches we perform below will end up reordering
# entries, which defeats the entire purpose of this script.
reMembV = re.compile( '\n[ \t]*((const|unsigned)[ \t]+)*(void|int|double|long|char|vtkIdType)([ \t]*\*)*[ \t]*([A-Za-z_][A-Za-z0-9_]*)[ \t]*(?:\[[ \t]*[0-9]*[ \t]*\][ \t]*)*;' );
membVarList = reMembV.findall( header )

reMembFP = re.compile( '\n[ \t]*([A-Za-z0-9_\*]+[ \t]+)*\([ \t]*\*[ \t]*([A-Za-z_][A-Za-z0-9_]*)[ \t]*\)' );
membFuncPtrList = reMembFP.findall( header )

members = []
for mv in membVarList:
	members.append( mv[-1] )
for fp in membFuncPtrList:
	members.append( fp[-1] )

memberpad = max( map( len, members ) )
lastmember = members[ len(members)-1 ]

#print "Members are:"
#for m in members:
#	print m

dm = {};
reDefn = re.compile( "\([ \t]*(\*|[A-Za-z0-9_]+)[ \t]*,[ \t]*(\*|[A-Za-z0-9_]+)[ \t]*,[ \t]*(\*|[A-Za-z0-9_]+)[ \t]*\)[ \t]*\n" )
reMemb = re.compile( "[ \t]*([A-Za-z0-9_]+)[ \t]*:(.*)\n" )

AllCellOps = {}
for S in shapes:
	for I in interps:
		for P in prodspaces:
			AllCellOps[ (S,I,P) ] = {}

eof = 0
dm = {}
while not eof:
	d = defines.readline()
	if d == '':
		eof = 1
	elif d[0] == '(':
		if len(dm) > 0:
			SaveMemberValues( (DomainShape, PolyInterpolant, PolyProdSpace), dm )
		m = reDefn.search( d )
		if m != None:
			DomainShape = m.groups()[0]
			PolyInterpolant = m.groups()[1]
			PolyProdSpace = m.groups()[2]
		dm = {}
	elif (d[0] == '#') or (d[0] == '\n'):
		continue
	else:
		m = reMemb.search( d )
		if m == None:
			print >> cellout, '#error "No match in %s for <%s>"' % (opspec_file,d)
		else:
			member = m.groups()[0]
			if member == 'EdgeCornerNodes':
			  value = Pad2DArray( m.groups()[1], 2, 12 )
			elif member == 'FaceCornerNodes':
			  value = Pad2DArray( m.groups()[1], 4, 6 )
			elif member == 'EdgesOfFaces':
			  value = Pad2DArray( m.groups()[1], 4, 6 )
			elif member == 'CornerParameterCoords':
				value = Pad2DArray( m.groups()[1], 3, 8, 0.0 )
			elif member == 'CellSimplices':
				value = Pad1DArray( m.groups()[1], 20 )
			else:
			  value = m.groups()[1]
			dm[ member ] = value
if len(dm) > 0:
	SaveMemberValues( (DomainShape, PolyInterpolant, PolyProdSpace), dm )

# OK, we've got cell specifications, it's time to write the output file.
# First, we'll write a preamble.
print >> cellout, '#include <math.h>\n'
print >> cellout, '#include <vtkCellOps.h>'
for sa in shapes:
	print >> cellout, '#include <Elements/%s.h>' % sa
print >> cellout, """
#include <Elements/Generic.h>
using namespace shoe;
// This file is autogenerated! Do NOT edit it or your changes will be lost.
// Instead, fix GenCellOps.py or CellOpsSpecifications.txt

"""


# Now output all the vtkCellOps structures
for pa in prodspaces:
	for ia in interps:
		for sa in shapes:
			k = (sa, ia, pa)
			if AllCellOps.has_key( k ):
				WriteOrderedMembers( sa, ia, pa, AllCellOps[ k ] )
print >> cellout, """
extern "C" int vtksnlRegisterElements()
{
  // Only register the finite elements once, no matter how many
  // times the function is called.
  static int beenHere = 0;
  if ( beenHere )
    {
    return 0;
    }
  beenHere = 1;

"""
for pa in prodspaces:
	for ia in interps:
		for sa in shapes:
			k = (sa, ia, pa)
			if AllCellOps.has_key( k ):
				WriteMemberInit( sa, ia, pa, AllCellOps[ k ] )
print >> cellout, """
  return 0;
}
"""

