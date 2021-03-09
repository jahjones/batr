# Barrier Analysis Tool for R (batr)
GIS tool in R for fragmenting rivers (polylines) with instream barriers (points).

Inputs: 
topologically consistent river or drainage network as polylines in sf friendly format e.g. ESRI shapefile
instream barriers or points aligned to the above network in sf friendly format e.g. ESRI shapefile

Outputs:
river network split into fragments between barriers 
points attributed with various river network attributes e.g. distance to mouth and upstream habitat
river network source and root (river mouth)
