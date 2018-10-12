# -*- coding: utf-8 -*-
'''
GISF2E ArcGIS - v2.00

To cite, please use: Karduni, A., Kermanshah, A., & Derrible, S., 2016 “A protocol to convert spatial polyline data to network formats and applications to world urban road networks”, Scientific Data, 3:160046
The article is available at: http://www.nature.com/articles/sdata201646

Learn more at https://csun.uic.edu/codes/GISF2E.html
'''

#import all libraries required
import arcpy
import csv

#data_workspace = arcpy.GetParameterAsText(0)
#input_streets = arcpy.GetParameterAsText(1)
data_workspace = "C:/Users/kermanah/Desktop/SD paper/Tamrin1/"
path = data_workspace + "\\"
arcpy.AddMessage(path)
input_streets = path + "Sofla_Olia.shp"
arcpy.AddMessage(input_streets)
arcpy.env.overwriteOutput = True

elevFields = ["bridge","tunnel"]

def class_calc(bridge, tunnel):
    if (bridge == 1 or bridge == "T") and (tunnel==0 or tunnel == "F"):
        return 1
    elif (tunnel==1 or tunnel == "T") and (bridge == 0 or bridge == "F"):
        return 2
    elif (bridge == 0 or bridge == "F") and (tunnel==0 or tunnel == "F"):
        return 0
    elif (bridge == 1 or bridge == "T") and (tunnel==1 or tunnel == "T") :
        return 3


fieldnames = [field.name for field in arcpy.ListFields(input_streets)]

def unique_values(table, field):
    with arcpy.da.SearchCursor(table, [field]) as cursor:
        return sorted({row[0] for row in cursor})

# feature_to_lines = arcpy.FeatureToLine_management(input_streets,path+ "feature_to_lines.shp")

arcpy.AddField_management(input_streets,"class","String")

input_streets_cursor = arcpy.UpdateCursor(input_streets)

for row in input_streets_cursor:
    if "tunnel" in fieldnames and "bridge" in fieldnames:
        bridge = row.getValue(elevFields[0])
        tunnel = row.getValue(elevFields[1])
        roadClass = class_calc(bridge,tunnel)
        row.setValue("class",roadClass)
        input_streets_cursor.updateRow(row)
    else:
        row.setValue("class",0)
        input_streets_cursor.updateRow(row)

classValues = unique_values(input_streets,"class")
streetClasses = []
nodeClasses = []
for classValue in classValues:
    print classValue
    streetClass = streetsLayer = arcpy.MakeFeatureLayer_management (input_streets, "streetsLayer", where_clause="class = '"+ str(classValue)+"'")
    streetClass = arcpy.Dissolve_management(streetClass)
    streetClassToLine = streetClass
    if classValue == "0":
        streetClassToLine = arcpy.FeatureToLine_management(streetClass,path+classValue+"_lines")
        streetClassNodes = arcpy.FeatureVerticesToPoints_management(streetClassToLine,path+classValue+"_nodes","BOTH_ENDS")
    else:
        streetClassToLine = arcpy.CopyFeatures_management(streetClass,path+classValue+"_lines")
        streetClassNodes = arcpy.FeatureVerticesToPoints_management(streetClassToLine,path+classValue+"_nodes","BOTH_ENDS")
    nodeClasses.append(streetClassNodes)
    streetClasses.append(streetClassToLine)
    arcpy.Delete_management(streetClass)
print("created nodes and edges");
mergedEdges = arcpy.Merge_management(streetClasses,path+"edges.shp")
arcpy.AddGeometryAttributes_management(mergedEdges,"LENGTH_GEODESIC","METERS")
arcpy.AddField_management(mergedEdges,"LENGTH")
arcpy.CalculateField_management(mergedEdges,"LENGTH","[LENGTH_GEO]")
mergedNodes = arcpy.Merge_management(nodeClasses,path+"nodes.shp")

for edge in streetClasses:
    arcpy.Delete_management(edge)
for node in nodeClasses:
    arcpy.Delete_management(node)

arcpy.AddMessage("Created nodes and edges shapefiles")

mergedNodes = arcpy.DeleteIdentical_management(mergedNodes,"Shape")

EdgeNodeSpatialJoin = arcpy.SpatialJoin_analysis(mergedEdges,mergedNodes,path+"edgeNodeSpatialJoin.shp","JOIN_ONE_TO_MANY","KEEP_ALL",match_option="INTERSECT",)

NodeEdgeSpatialJoin = arcpy.SpatialJoin_analysis(mergedNodes,EdgeNodeSpatialJoin,path+"NodeEdgeSpatialJoin.shp","JOIN_ONE_TO_MANY","KEEP_ALL",match_option="INTERSECT")

arcpy.AddField_management(NodeEdgeSpatialJoin,"Start_Node")
arcpy.AddField_management(NodeEdgeSpatialJoin,"End_Node")
arcpy.AddField_management(NodeEdgeSpatialJoin,"Edge")
arcpy.CalculateField_management(NodeEdgeSpatialJoin,"Start_Node","[TARGET_F_1]")
arcpy.CalculateField_management(NodeEdgeSpatialJoin,"End_Node","[JOIN_FID]")
arcpy.CalculateField_management(NodeEdgeSpatialJoin,"Edge","[TARGET_FID]")

arcpy.AddMessage("created Edge list")
arcpy.AddMessage("Writing edge list to file")

edgesCheck = {}
outputEdgelist = []

outCursor = arcpy.SearchCursor(NodeEdgeSpatialJoin)

for row in outCursor:
    startNode = row.getValue("Start_Node")
    endNode = row.getValue("End_Node")
    edge = row.getValue("Edge")
    length = row.getValue("LENGTH")
    if startNode != endNode and edge not in edgesCheck:
        edgesCheck[edge] = 1
        edgelist = (startNode,endNode,edge,length)
        outputEdgelist.append(edgelist)

# outputEdgelist = set(outputEdgelist)
with open(path+"edgelist.csv","w") as csvFile:
    csvWriter = csv.writer(csvFile, lineterminator="\n")
    csvWriter.writerow(["Start_Node","End_Node","Edge","Length"])
    arcpy.AddMessage(outputEdgelist)
    csvWriter.writerows(outputEdgelist)

arcpy.Delete_management(EdgeNodeSpatialJoin)
arcpy.Delete_management(NodeEdgeSpatialJoin)
