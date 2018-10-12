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
import datetime

data_workspace = arcpy.GetParameterAsText(0)
input_streets = arcpy.GetParameterAsText(1)
print (data_workspace)
print (input_streets)
arcpy.AddMessage("starting GISF2E!")
arcpy.AddMessage(datetime.datetime.now())


test = False
if not data_workspace and not input_streets:
    data_workspace= raw_input("Please enter path output folder: ")
    input_streets = raw_input("Please enter path for input line shapefile: ")
    path = data_workspace + "\\"

arcpy.AddMessage("creating nodes, edges, and shapefiles at:")
arcpy.AddMessage(path)
arcpy.AddMessage("using centerline firle:")
arcpy.AddMessage(input_streets)
arcpy.env.overwriteOutput = True

elevFields = ["bridge","tunnel"]

#function to get a line class attribute from tunnel and bridge fields
def class_calc(bridge, tunnel):
    if (bridge == 1 or bridge == "T") and (tunnel==0 or tunnel == "F"):
        return 1
    elif (tunnel==1 or tunnel == "T") and (bridge == 0 or bridge == "F"):
        return 2
    elif (bridge == 0 or bridge == "F") and (tunnel==0 or tunnel == "F"):
        return 0
    elif (bridge == 1 or bridge == "T") and (tunnel==1 or tunnel == "T") :
        return 3

#function to get unique values from a field
def unique_values(table, field):
    with arcpy.da.SearchCursor(table, [field]) as cursor:
        return sorted({row[0] for row in cursor})


#getting a list of field names from shape files
fieldnames = [field.name for field in arcpy.ListFields(input_streets)]


# adding class attribute to line shapefile
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
#seperating line by class, disolving the lines and then seperating the edges based on the intersections
for classValue in classValues:
    streetClass = streetsLayer = arcpy.MakeFeatureLayer_management (input_streets, "streetsLayer", where_clause="class = '"+ str(classValue)+"'")
    streetClass = arcpy.Dissolve_management(streetClass)
    streetClassToLine = arcpy.FeatureToLine_management(streetClass,path+classValue+"_lines")
    streetClassNodes = arcpy.FeatureVerticesToPoints_management(streetClassToLine,path+classValue+"_nodes","BOTH_ENDS")
    nodeClasses.append(streetClassNodes)
    streetClasses.append(streetClassToLine)
    arcpy.Delete_management(streetClass)

arcpy.AddMessage("Different elevation layers seperated")

#merging the results, deleting the seperate line files
mergedEdges = arcpy.Merge_management(streetClasses,path+"edges.shp")
arcpy.AddGeometryAttributes_management(mergedEdges,"LENGTH_GEODESIC","METERS")
arcpy.AddField_management(mergedEdges,"LENGTH")
arcpy.CalculateField_management(mergedEdges,"LENGTH","[LENGTH_GEO]")
mergedNodes = arcpy.Merge_management(nodeClasses,path+"nodes.shp")
for edge in streetClasses:
    arcpy.Delete_management(edge)
for node in nodeClasses:
    arcpy.Delete_management(node)
    
mergedNodes = arcpy.DeleteIdentical_management(mergedNodes,"Shape")

arcpy.AddMessage("Created nodes and edges shapefiles")

# spatial join between edges file and nodes file to find start nodes and end nodes for each edge
EdgeNodeSpatialJoin = arcpy.SpatialJoin_analysis(mergedEdges,mergedNodes,path+"edgeNodeSpatialJoin.shp","JOIN_ONE_TO_MANY","KEEP_ALL",match_option="INTERSECT",)

arcpy.AddMessage("edge node spatial join completed")

# spatial join nodes to the previous spatial join to get start node, end node and edges
NodeEdgeSpatialJoin = arcpy.SpatialJoin_analysis(mergedNodes,EdgeNodeSpatialJoin,path+"NodeEdgeSpatialJoin.shp","JOIN_ONE_TO_MANY","KEEP_ALL",match_option="INTERSECT")

arcpy.AddMessage("node edge spatial join completed")


#deleting duplocate edges
outputEdgelist = []
arcpy.DeleteIdentical_management(NodeEdgeSpatialJoin,"TARGET_F_1")
outCursor = arcpy.SearchCursor(NodeEdgeSpatialJoin)
arcpy.AddMessage("creating Edge list")

#creating edgelist
for row in outCursor:
    startNode = row.getValue("TARGET_F_1")
    endNode = row.getValue("JOIN_FID")
    edge = row.getValue("TARGET_FID")
    length = row.getValue("LENGTH")
    if startNode != endNode:
        edgelist = (startNode,endNode,edge,length)
        outputEdgelist.append(edgelist)

arcpy.AddMessage("created Edge list")
arcpy.AddMessage("Writing edge list to file")

#writing results to edgelist csv
with open(path+"edgelist.csv","w") as csvFile:
    csvWriter = csv.writer(csvFile, lineterminator="\n")
    csvWriter.writerow(["Start_Node","End_Node","Edge","Length"])
    csvWriter.writerows(outputEdgelist)
    
arcpy.AddMessage("Deleting extra files")
arcpy.Delete_management(EdgeNodeSpatialJoin)
arcpy.Delete_management(NodeEdgeSpatialJoin)
arcpy.AddMessage(datetime.datetime.now())
arcpy.AddMessage("Done!")
