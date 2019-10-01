# Bathymetry base data

p = aegis.bathymetry::bathymetry_parameters( DS="bathymetry" )
bathymetry.db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015

M = bathymetry.db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
