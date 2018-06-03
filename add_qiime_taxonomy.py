
import sys
import cgi
#import cgitb
#cgitb.enable()
import csv
import numpy as np
import json
import os.path
import collections
from operator import itemgetter
import scipy.stats as stats
from scipy.spatial import distance
import math
import re
import skbio
import skbio.diversity as skdiv

def sorted_alnum( l ):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


def loaddata( fname ):
	reader = csv.reader(open( fname, "r"), delimiter='\t')
	data = list(reader)
	volumes = []
	for i in range( 0, len( data[0] ) ):
		if len( data[0][i] ) > 0:
			mn = i
			c0 = data[0][i][0]
			if c0 >= '0' and c0 <= '9':
				ml = i
			else:
				volumes.append( data[0][i] )
	return ( data, volumes, mn, ml )

def loadtaxonomy( data, ml, spfilter, ilevel ):
	kdict = {}
	kdnames = {}
	kgnames = {}
	cm = 0
	knorder = []
	kdata = {}
	if isinstance( spfilter, (list,)):
		excludelist = []
		reverse = 0
	else:
		excludelist = spfilter[0]
		reverse = spfilter[1]
	for d in data[1:]:
		if len( d ) < len( data[0] ):
			continue
		ckey = "".join( d[0:ilevel] )
		ckname = ""
		cgname = d[0]
		exflag = reverse
		if len( excludelist ) > 0:
			for k in range( ml ):
				if len( d[k] ) > 0 and d[k] in excludelist:
					exflag = 1 - reverse
					break
		if exflag == 1:
			continue
		for k in range( ilevel - 1, -1, -1 ):
			if len( d[k] ) > 0 and d[k] != "*":
				ckname = d[k]
				break
		if ml > 0 and len( d[1] ) > 0:
			cgname = d[1]
		if not ckey in kdict:
			kdict[ ckey ] = cm
			kdnames[ ckey ] = ckname
			kgnames[ ckey ] = cgname
			cm = cm + 1
			knorder.append( ckey )
			kdata[ ckey ] = d[0:ilevel] 
	return ( kdict, kdnames, kgnames, knorder, kdata )

def loadqtaxonomy( ifile ):
	rdic = {}
	maxtlevel = 9
	with open( ifile ) as f:
		for line in f:
			ls = line.strip().split( "\t" )
			ls2 = ls[1].split( ";" )
			ctax = [ s.strip()[ 3: ] for s in ls2 ]
			#maxtlevel = max( maxtlevel, min( len( ctax ), 9 ) )
			for k in range( len( ctax ), 9 ):
				ctax.append( "" )
			rdic[ ls[0] ] = ctax[ :9]
	return ( rdic, maxtlevel )


def processtags_m( volumes, tags, dfilter ):
	if dfilter == "none":
		findex = range( 0, len( volumes ) )
		mtags = tags
	else:
		findex = []
		for j in range( len( volumes ) ):
			if len( tags[ dfilter ][ j ] ):
				findex.append( j )
		mtags = {}
		for tkey in tags:
			mtval = []
			for j in range( len( volumes ) ):
				if j in findex: 				
					mtval.append( tags[tkey][ j ] ) 
			mtags[ tkey ] = mtval
	return ( findex, mtags )

def load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames ):
	edata = []
	site_ids = []
	species_ids = [ "" ] * len( kdict )
	dim2_ids = species_ids
	for i in range( ml + 1, mn + 1 ):
		if ( i - ml - 1 ) in findex:
			crow = [ 0 ] * len( kdict )
			site_ids.append( volumes[ i - ml - 1 ] )
			for d in data[1:]:
				ckey = "".join( d[0:ilevel] )
				if ckey in kdict:
					cind = kdict[ ckey ]
					crow[ cind ] += int( float( d[i] ) )
					if len( site_ids ) == 1:
						species_ids[ cind ] = kdnames[ ckey ] 
			edata.append( crow )
	return ( edata, site_ids, species_ids )


if len( sys.argv ) < 2:
	print "arguments: emap tax"
	sys.exit( 1 )
	
( data, volumes, mn, ml ) = loaddata( sys.argv[1] )
ilevel = ml + 1
dfilter = ""
( kdict, kdnames, kgnames, knorder, kdata ) = loadtaxonomy( data, ml, [], ilevel )
#( findex, mtags ) = processtags_m( volumes, tags, dfilter )
findex = range( 0, len( volumes ) )
( edata, site_ids, species_ids ) = load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )
( qtax, maxtlevel ) = loadqtaxonomy( sys.argv[2] )

#print >> sys.stderr, kdata

ostr = ""
for i in range( ilevel + maxtlevel ):
	ostr += str( i + 1 ) + "\t"
for vol in site_ids:
	ostr += vol + "\t"
print ostr
for knum in range( len( knorder ) ):
	key = knorder[ knum ]
	ostr = ""
	for ct in qtax[ kdata[ key ][ 0 ] ]:
		ostr += ct + "\t"
	for ct in kdata[ key ]:
		ostr += ct + "\t"
	for vnum in range( len( site_ids ) ):
		ostr += str( edata[ vnum ][ knum ] ) + "\t"
	print ostr
		
		
	



