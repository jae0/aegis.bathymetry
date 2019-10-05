



# ------------------------------------------------
# Create unerlying areal units and data for modelling of bathymetry, etc using lattice methods via carstm

for ( areal_units_resolution_km in c(10, 20, 25) ) {
  # for ( spatial_domain in c("snowcrab", "SSE")) {
   for ( spatial_domain in c("snowcrab", "SSE")) {
    p = aegis.bathymetry::bathymetry_parameters(
      project_class = "carstm", # defines which parameter class / set to load
      spatial_domain = spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      areal_units_resolution_km = areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
    )
    # sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    # M = bathymetry.db( p=p, DS="aggregated_data", redo=TRUE )  # already done in 01.bathymetry_data.R will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
    M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  }
}



# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p = aegis.bathymetry::bathymetry_parameters(
  project_class = "carstm", # defines which parameter class / set to load
  spatial_domain = "snowcrab",  # defines spatial area, currenty: "snowcrab" or "SSE"
  inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
  areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
)



if (0) {

  # example sequence to force creating of input data for modelling

  sppoly = areal_units( p=p, redo=TRUE ); plot(sppoly) # or: spplot( sppoly, "StrataID", main="StrataID", sp.layout=p$coastLayout )
  M = bathymetry.db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  str(M)

  res = bathymetry_carstm( p=p, DS="carstm_modelled", redo=TRUE ) # run model and obtain predictions

  # loading saved results
  res = bathymetry_carstm( p=p, DS="carstm_modelled" ) # to load currently saved sppoly
  fit = bathymetry_carstm( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  plot(fit)
  plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

  s = summary(fit)
  s$dic$dic
  s$dic$p.eff

  # maps of some of the results
  carstm_plot( p=p, res=res, vn="z.predicted" )
  carstm_plot( p=p, res=res, vn="z.random_strata_nonspatial" )
  carstm_plot( p=p, res=res, vn="z.random_strata_spatial" )

}

# end
