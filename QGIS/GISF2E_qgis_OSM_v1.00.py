# -*- coding: utf-8 -*-
'''
GISF2E QGIS OSM - v1.00 - Note that the code was developed with pandas version 0.14.1 and it runs much faster for this version than with later versions.

To cite, plese use: Karduni, A., Kermanshah, A., & Derrible, S., 2016 “A protocol to convert spatial polyline data to network formats and applications to world urban road networks”, Scientific Data, 3:160046
The article is available at: http://www.nature.com/articles/sdata201646

Main authors: Felipe Macena, Rodrigo Marinho
'''

#import all libraries required
import shapefile 
import pandas as pd
from pyproj import Proj
import osgeo.ogr, osgeo.osr
from osgeo import ogr
import sys 
import os
import time
from math import pi, sin, cos, tan, sqrt
from qgis.core import *
def LLtoUTM(ReferenceEllipsoid, Lat, Long, zone = None):
    """converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees
    Written by Chuck Gantz- chuck.gantz@globalstar.com"""

    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    k0 = 0.9996

    #Make sure the longitude is between -180.00 .. 179.9
    LongTemp = (Long+180)-int((Long+180)/360)*360-180 # -180.00 .. 179.9

    LatRad = Lat*_deg2rad
    LongRad = LongTemp*_deg2rad

    if zone is None:
        ZoneNumber = int((LongTemp + 180)/6) + 1
    else:
        ZoneNumber = zone

    if Lat >= 56.0 and Lat < 64.0 and LongTemp >= 3.0 and LongTemp < 12.0:
        ZoneNumber = 32

    # Special zones for Svalbard
    if Lat >= 72.0 and Lat < 84.0:
        if  LongTemp >= 0.0  and LongTemp <  9.0:ZoneNumber = 31
        elif LongTemp >= 9.0  and LongTemp < 21.0: ZoneNumber = 33
        elif LongTemp >= 21.0 and LongTemp < 33.0: ZoneNumber = 35
        elif LongTemp >= 33.0 and LongTemp < 42.0: ZoneNumber = 37

    LongOrigin = (ZoneNumber - 1)*6 - 180 + 3 #+3 puts origin in middle of zone
    LongOriginRad = LongOrigin * _deg2rad

    #compute the UTM Zone from the latitude and longitude
    UTMZone = "%d%c" % (ZoneNumber, _UTMLetterDesignator(Lat))

    eccPrimeSquared = (eccSquared)/(1-eccSquared)
    N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))
    T = tan(LatRad)*tan(LatRad)
    C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
    A = cos(LatRad)*(LongRad-LongOriginRad)

    M = a*((1
            - eccSquared/4
            - 3*eccSquared*eccSquared/64
            - 5*eccSquared*eccSquared*eccSquared/256)*LatRad
           - (3*eccSquared/8
              + 3*eccSquared*eccSquared/32
              + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
           + (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad)
           - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad))

    UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
                        + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
                  + 500000.0)

    UTMNorthing = (k0*(M+N*tan(LatRad)*(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                                        + (61
                                           -58*T
                                           +T*T
                                           +600*C
                                           -330*eccPrimeSquared)*A*A*A*A*A*A/720)))
    if Lat < 0:
        UTMNorthing = UTMNorthing + 10000000.0; #10000000 meter offset for southern hemisphere
    return (UTMZone, UTMEasting, UTMNorthing)


def _UTMLetterDesignator(Lat):
    """This routine determines the correct UTM letter designator for the given
    latitude returns 'Z' if latitude is outside the UTM limits of 84N to 80S
    Written by Chuck Gantz- chuck.gantz@globalstar.com"""
    if Lat >= 0: return 'N'
    else: return 'S'

def add_link(start,end,cont):
    pontos = []
    pontos.append(osgeo.ogr.Geometry(osgeo.ogr.wkbLineString))
    lista = list([testando['coordinates_utm'][start:end+1]][0])
    for a in lista:        
        pontos[-1].AddPoint(a[0],a[1])  
    featureIndex = cont 
    feature = osgeo.ogr.Feature(layer_defn)
    feature.SetGeometry(pontos[-1])
    feature.SetFID(featureIndex)
    layer.CreateFeature(feature)
    for name in range(len(field_names)):
        if field_names[name] != 'Edge' and field_names[name] != 'Length':
            feature.SetField(name, str(testando[field_names[name]][start_node]))
        elif field_names[name] == 'Edge':
            feature.SetField(name, cont+1) 
        elif field_names[name] == 'Length':
            feature.SetField(name, distance(start,end))             
        layer.SetFeature(feature)
def add_field(field_names,featureIndex):
    for name in range(len(field_names)):
            if field_names[name] == 'coord_x' or field_names[name] == 'coord_y':
                name_aux = 'coordinates_utm'
            else:
                name_aux = field_names[name]
            feature = layer.GetFeature(featureIndex) #lets get the first feature (FID=='0')
            if field_names[name] == 'coord_x':       
                feature.SetField(name,testando_nodes[name_aux][featureIndex][0])
            elif field_names[name] == 'coord_y':
                feature.SetField(name,testando_nodes[name_aux][featureIndex][1])
            else:
                feature.SetField(name,int(testando_nodes[name_aux][featureIndex]))
            layer.SetFeature(feature) #now make the change permanent
def distance(start,end):
    dist = 0
    for u in range(start,end):
            dist = dist + ((testando['coordinates_utm'][u][0] - testando['coordinates_utm'][u+1][0])**2 + (testando['coordinates_utm'][u][1] - testando['coordinates_utm'][u+1][1])**2)**0.5
    return dist
def Edgelist(start,end,edge,length):
    global testando, columns
    datasheet = pd.DataFrame([[testando['coordinates_utm'][start][0],testando['coordinates_utm'][start][1],testando['coordinates_utm'][end][0],testando['coordinates_utm'][end][1],testando['Node'][start],testando['Node'][end],edge,length]],columns = columns, index = None)
    return datasheet

_deg2rad, _rad2deg, _EquatorialRadius, _eccentricitySquared  = pi / 180.0, 180.0 / pi, 2, 3

_ellipsoid = [
#  id, Ellipsoid name, Equatorial Radius, square of eccentricity
# first once is a placeholder only, To allow array indices to match id numbers
	[ -1, "Placeholder", 0, 0],
	[ 1, "Airy", 6377563, 0.00667054],
	[ 2, "Australian National", 6378160, 0.006694542],
	[ 3, "Bessel 1841", 6377397, 0.006674372],
	[ 4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],
	[ 5, "Clarke 1866", 6378206, 0.006768658],
	[ 6, "Clarke 1880", 6378249, 0.006803511],
	[ 7, "Everest", 6377276, 0.006637847],
	[ 8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],
	[ 9, "Fischer 1968", 6378150, 0.006693422],
	[ 10, "GRS 1967", 6378160, 0.006694605],
	[ 11, "GRS 1980", 6378137, 0.00669438],
	[ 12, "Helmert 1906", 6378200, 0.006693422],
	[ 13, "Hough", 6378270, 0.00672267],
	[ 14, "International", 6378388, 0.00672267],
	[ 15, "Krassovsky", 6378245, 0.006693422],
	[ 16, "Modified Airy", 6377340, 0.00667054],
	[ 17, "Modified Everest", 6377304, 0.006637847],
	[ 18, "Modified Fischer 1960", 6378155, 0.006693422],
	[ 19, "South American 1969", 6378160, 0.006694542],
	[ 20, "WGS 60", 6378165, 0.006693422],
	[ 21, "WGS 66", 6378145, 0.006694542],
	[ 22, "WGS-72", 6378135, 0.006694318],
	[ 23, "WGS-84", 6378137, 0.00669438]
]



#dividing each link by their osm_id information (osm_maps) or id information (non-osm_maps)

sf = shapefile.Reader("C:\Users\Rodrigo\Desktop\GIStool\input\Chicago_cut.shp")
shapes = sf.shapes()
fields=sf.fields
fields = [field[0] for field in sf.fields[1:]]
start = time.time()
print '\n'+str(time.time() - start)+'s - Reading and arranging the root file information...'
sys.stdout.flush()
#here you can change the utm zone depending on where the map is from
#p = Proj(proj='utm',zone=int(utm_zone[:-1]),ellps='WGS84')

if 'TYPE' in fields or 'Type' in fields:
    try:
        fields[fields.index('TYPE')] = 'type'
    except:
        fields[fields.index('Type')] = 'type'
      
heads = [x.upper() if x in ['Group','coordinates','primeiro','coordinates_utm','Node'] else x for x in [sf.shapeRecords()[0].shape.__geo_interface__.keys()[0]]+fields] + ['Group','coordinates','primeiro','coordinates_utm','Node']
type_ = False
if 'type' in map(lambda x: x.lower(),fields):
    heads.remove('type')
    heads.insert(3,'type 1')
    type_ = True

dicti = {}
for i in range(len(heads)-1):
    dicti[heads[i]] = []

testando = pd.DataFrame()

n = 0 
for feature in sf.shapeRecords():
    geom = feature.shape.__geo_interface__
    atr = dict(zip(fields, feature.record))
    if type_:
        atr['type 1'] = atr.pop('type')
    geom.update(atr)
    geom['Group'] = n    
    j = geom['coordinates']     
    for i in range(len(heads)-4):
        u = [geom[heads[i]]]
        for y in range(len(j)):
            dicti[heads[i]].append(u[0])      
    dicti['coordinates'] = list(dicti['coordinates']) + list(j)
    dicti['primeiro'] = dicti['primeiro'] + [1] + [0]*(len(j) - 2) + [1]
    n += 1
    if ((n-1) % 800) == 0:
        dicti['coordinates_utm'] = None
        dicti['Node'] = None
        testando = pd.concat([testando,pd.DataFrame.from_dict(dicti)])
        for i in range(len(heads)-1):
            dicti[heads[i]] = []
if ((n-1) % 800) != 0:
    dicti['coordinates_utm'] = None
    dicti['Node'] = None
    testando = pd.concat([testando,pd.DataFrame.from_dict(dicti)],ignore_index=True)     

utm_zone = LLtoUTM(23, testando['coordinates'][0][1], testando['coordinates'][0][0])[0]
p = Proj(proj='utm',zone=int(utm_zone[:-1]),ellps='WGS84')

testando['coordinates_utm'] = [p(testando['coordinates'][k][0],testando['coordinates'][k][1]) for k in range(len(testando['coordinates']))]
testando['Node'] = None

#testando = pd.DataFrame.from_dict(dicti)
dicti = None
geom = None  
atr = None
#test before creating nodes for OSM maps    
a= testando.groupby('coordinates_utm').groups
testando['continuidade'] = 1
testando.ix[testando.primeiro == 1,'continuidade'] = 0; testando
b1 = map(lambda x: (testando['coordinates_utm'][x],[x]),testando.groupby('continuidade').groups[0])

b = filter(lambda x: (testando['bridge'][x[1][0]] == testando['bridge'][int(x[1][1])]) or (testando['tunnel'][x[1][0]] == testando['tunnel'][int(x[1][1])]), filter(lambda x: len(x[1]) >= 2, zip(a.keys(),a.values())))

for h in (b):
    for k in range(len(h[1])):
        testando['continuidade'][h[1][k]] = 2
b = b + b1
points = sorted(testando.groupby('continuidade').groups[0] + testando.groupby('continuidade').groups[2])
testando['Node'][points] = range(1,len(points)+1)

print '\n'+str(time.time() - start)+'s - Reading and arranging are done!'
sys.stdout.flush()

#creating nodes shapefile
path='C:\Users\Rodrigo\Desktop\GIStool\output\Chicago_cut_nodes.shp'
spatialReference = osgeo.osr.SpatialReference() #will create a spatial reference locally to tell the system what the reference will be
spatialReference.ImportFromProj4('+proj=utm +zone='+utm_zone+' ellps=WGS84 +datum=WGS84 +units=m') #here we define this reference to be utm Zone 16N with wgs84
driver = osgeo.ogr.GetDriverByName('ESRI Shapefile') # will select the driver for our shp-file creation.
shapeData = driver.CreateDataSource(path) #so there we will store our data
layer = shapeData.CreateLayer('Nodes', spatialReference, osgeo.ogr.wkbPoint) #this will create a corresponding layer for our data with given spatial information
layer_defn = layer.GetLayerDefn() # gets parameters of the current shapefile
point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)

new_field = ogr.FieldDefn('Node', ogr.OFTInteger) #and a forth field 'oneway' stored as integer
layer.CreateField(new_field) #self explaining
new_field = ogr.FieldDefn('coord_x', ogr.OFTReal) #and a forth field 'oneway' stored as integer
layer.CreateField(new_field) #self explaining
new_field = ogr.FieldDefn('coord_y', ogr.OFTReal) #and a forth field 'oneway' stored as integer
layer.CreateField(new_field) #self explaining

new_field = ogr.FieldDefn('bridge', ogr.OFTInteger) #we will create a new field called bridge as integer
layer.CreateField(new_field) #self explaining
new_field = ogr.FieldDefn('tunnel', ogr.OFTInteger) #and a sixth field 'tunnel' stored as integer
layer.CreateField(new_field) #self explaining
field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
aux_b, featureIndex = list(zip(*b)[0]), 0

testando_nodes = testando[testando.continuidade != 1].drop_duplicates(['coordinates_utm'])
testando_nodes.index = range(len(testando_nodes))

t, aux_b, a, b, b1, points = None, None, None, None, None, None

print '\n'+str(time.time() - start)+'s - Adding nodes...'
sys.stdout.flush()   

for featureIndex in range(len(testando_nodes)):
        point.AddPoint(testando_nodes['coordinates_utm'][featureIndex][0],testando_nodes['coordinates_utm'][featureIndex][1]) #create a new point at given coordinates
        feature = osgeo.ogr.Feature(layer_defn)
        feature.SetGeometry(point)
        feature.SetFID(featureIndex)
        layer.CreateFeature(feature)
        add_field(field_names,featureIndex)    
shapeData = None

new_b = None
testando_nodes = None

print '\n'+str(time.time() - start)+'s - Nodes are done!'
sys.stdout.flush()

#Cut the raw links and Creat the edge shapefile
columns = ['X_Coord_Start','Y_Coord_Start','X_Coord_End','Y_Coord_End','Start_Node','End_Node','Edge','Length']
edgelist = pd.DataFrame(columns = columns, index=None)

#creating edges shapefile
path='C:\Users\Rodrigo\Desktop\GIStool\output\Chicago_cut_edges.shp'
spatialReference = osgeo.osr.SpatialReference() #will create a spatial reference locally to tell the system what the reference will be
#here you can change the utm zone depending on where the map is from
spatialReference.ImportFromProj4('+proj=utm +zone='+utm_zone+' ellps=WGS84 +datum=WGS84 +units=m') #here we define this reference to be utm Zone 16N with wgs84
driver = osgeo.ogr.GetDriverByName('ESRI Shapefile') # will select the driver for our shp-file creation.
shapeData = driver.CreateDataSource(path) #so there we will store our data
layer = shapeData.CreateLayer('linestrings', spatialReference, osgeo.ogr.wkbLineString) #this will create a corresponding layer for our data with given spatial information.
layer_defn = layer.GetLayerDefn() # gets parameters of the current shapefile

for f in ['Edge']+heads[:-5]+['Length']: 
    fieldDefn = ogr.FieldDefn(f, ogr.OFTString)
    layer.CreateField(fieldDefn)
field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]

print '\n'+str(time.time() - start)+'s - Adding links and creating the Edgelist...'
sys.stdout.flush()

links = (testando.groupby(['Group']).groups)
links1 = (testando.groupby(['Group','continuidade']).groups)
cont = 0   
for link in range(len(links)):
    try:
        length_0 =  len(links1[link,0])
    except:
        length_0 = -1
    try:
        length_1 =  len(links1[link,1])
    except:
        length_1 = -1
    try:
        length_2 =  len(links1[link,2])
    except:
        length_2 = -1      
    if length_0 == 1 and length_2 == 1:
        start_node, end_node = min(links1[link,0][0],links1[link,2][0]), max(links1[link,0][0],links1[link,2][0])
        add_link(start_node,end_node,cont)
        edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
        cont += 1
    elif length_0 == -1:
        for i in range(length_2-1):
            start_node, end_node = links1[link,2][i], links1[link,2][i+1]
            add_link(start_node,end_node,cont)
            edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
            cont += 1
    elif length_2 == -1:
        start_node, end_node = links1[link,0][0], links1[link,0][1]
        add_link(start_node,end_node,cont)
        edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
        cont += 1
    else:
        if length_0 > 1:
            start_node, end_node = links1[link,0][0], links1[link,2][0]
            add_link(start_node,end_node,cont)
            edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
            cont += 1
            for i in range(length_2-1):
                start_node, end_node = links1[link,2][i], links1[link,2][i+1]
                add_link(start_node,end_node,cont)
                edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
                cont += 1
            start_node, end_node = links1[link,2][-1], links1[link,0][1]
            add_link(start_node,end_node,cont)
            edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
            cont += 1                  
        else:
            if links1[link,0][0] < links1[link,2][0]:
                start_node, end_node = links1[link,0][0], links1[link,2][0]
                add_link(start_node,end_node,cont)
                edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
                cont += 1 
                for i in range(length_2 - 1):
                    start_node, end_node = links1[link,2][i], links1[link,2][i+1]
                    add_link(start_node,end_node,cont)
                    edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True) 
                    cont += 1
            else:
                for i in range(length_2 - 1):
                    start_node, end_node = links1[link,2][i], links1[link,2][i+1] 
                    add_link(start_node,end_node,cont)  
                    edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
                    cont += 1                  
                start_node, end_node = links1[link,2][-1], links1[link,0][0]
                add_link(start_node,end_node,cont)
                edgelist = edgelist.append(Edgelist(start_node,end_node,cont+1,distance(start_node,end_node)),ignore_index=True)
                cont += 1                 
shapeData = None       

print '\n'+str(time.time() - start)+'s - Links and Edgelist are done!'
sys.stdout.flush()
    
print '\n'+str(time.time() - start)+'s - Saving Edgelist...'     
sys.stdout.flush()

edgelist.to_csv('C:\Users\Rodrigo\Desktop\GIStool\output\Chicago_cut_edgelist.csv',index=False)

wb = QgsVectorLayer('C:\Users\Rodrigo\Desktop\GIStool\output\Chicago_cut_edges.shp', 'Chicago_cut_edges', 'ogr')
QgsMapLayerRegistry.instance().addMapLayer(wb)
wb = QgsVectorLayer('C:\Users\Rodrigo\Desktop\GIStool\output\Chicago_cut_nodes.shp', 'Chicago_cut_nodes', 'ogr')
QgsMapLayerRegistry.instance().addMapLayer(wb)
print '\n'+str(time.time() - start)+'s - You are all set.'
sys.stdout.flush()
