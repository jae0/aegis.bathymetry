
# Bathymetry base data
p = aegis.bathymetry::bathymetry_parameters( DS="bathymetry" )
bathymetry.db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015


# The rest are for underlying areal units and data for modelling of bathymetry via carstm lattices
for ( areal_units_resolution_km in c(10, 20, 25) ) {
  # for ( spatial_domain in c("snowcrab", "SSE")) {
   for ( spatial_domain in c("snowcrab", "SSE")) {
    p = aegis.bathymetry::bathymetry_parameters(
      project_class = "carstm", # defines which parameter class / set to load
      project_name = "bathymetry",
      variabletomodel = "z",
      spatial_domain = spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      areal_units_resolution_km = areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
    )
    sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    M = bathymetry_carstm( p=p, DS="aggregated_data", redo=TRUE )  # already done in 01.bathymetry_data.R will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
    M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  }
}
