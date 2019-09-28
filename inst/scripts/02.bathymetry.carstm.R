



# ------------------------------------------------
# Create unerlying areal units and data for modelling of bathymetry, etc using lattice methods via carstm

for ( areal_units_resolution_km in c(10, 20, 25) ) {
  # for ( spatial_domain in c("snowcrab", "SSE")) {
   for ( spatial_domain in c("snowcrab", "SSE")) {
    areal_units_overlay = "snowcrab_managementareas"
    if ( spatial_domain=="SSE") areal_units_overlay = "groundfish_strata"
    p = aegis.bathymetry::bathymetry_parameters(
      project_class = "carstm", # defines which parameter set to load
      id = paste("bathymetry", areal_units_overlay, sep="_"),
      inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      spatial_domain = spatial_domain,  # defines spatial area
      areal_units_resolution_km = areal_units_resolution_km, # km dim of lattice
      areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
      areal_units_strata_type = "lattice", # "aegis_lattice" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_overlay = areal_units_overlay, # additional polygon layers for subsequent analysis such as management area: "snowcrab" or "groundfish"  # for now ..
      areal_units_constraint="none", # set[, c("lon", "lat")],  # to limit to sppoly to only those with data that fall into them
       libs = RLibrary ( "sp", "spdep", "rgeos", "spatialreg", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry", "carstm" )
    )

    sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    M = bathymetry_carstm( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
    p$constant_offset = 2500
    M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  }
}




project = "snowcrab"


# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p = aegis.bathymetry::bathymetry_parameters(
  project_class = "carstm", # defines which parameter class / set to load
  id = paste("bathymetry", project, sep="_"),  # label to tag the results
  inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
  spatial_domain = "snowcrab",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_strata_type = "lattice", # "aegis_lattice" to use ageis fields instead of carstm fields ... note variables are not the same
  areal_units_overlay = "snowcrab_managementareas", # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
  # areal_units_resolution_km = 10, # km dim of lattice ~ 16 hrs
  areal_units_resolution_km = 20, # km dim of lattice ~ 1 hr
  # areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  # areal_units_proj4string_planar_km = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  carstm_modelengine = "inla.default",  # {model engine}.{label to use to store}
  carstm_modelcall = '
    inla(
      formula = z ~ 1
        + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
        + f(iid_error, model="iid", hyper=H$iid),
      family = "lognormal", # "zeroinflatedpoisson0",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),  # config=TRUE if doing posterior simulations
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      control.fixed=H$fixed,  # priors for fixed effects, generic is ok
      control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
      # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
      num.threads=4,
      blas.num.threads=4,
      verbose=TRUE
    ) ',
  # carstm_modelcall = 'glm( formula = z ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ] ) ',  # for modelengine='glm'
  # carstm_modelcall = 'gam( formula = z ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ] ) ',  # for modelengine='gam'
  constant_offset = 2500, # pre-modelling transformations
  libs = RLibrary ( "sp", "spdep", "rgeos", "spatialreg", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry", "carstm" )
)





if (0) {
  # example sequence to force creating of input data for modelling
  p = c(p, aegis.coastline::coastline_layout( p=p, redo=TRUE ) )  # set up default map projection
  sppoly = areal_units( p=p, redo=TRUE )
  plot(sppoly) # or: spplot( sppoly, "StrataID", main="StrataID", sp.layout=p$coastLayout )
  M = bathymetry_carstm( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  str(M)
}

# run model and obtain predictions
sppoly = bathymetry_carstm( p=p, DS="carstm_modelled", redo=TRUE ) # extract predictions in sppoly

if (0) {
  sppoly = bathymetry_carstm( p=p, DS="carstm_modelled" ) # to load currently saved sppoly
  fit =    bathymetry_carstm( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit

  plot(fit)

  s = summary(fit)

  s$dic$dic
  s$dic$p.eff

  # maps of some of the results
  p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
  p$mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")
  p = c(p, aegis.coastline::coastline_layout( p=p ) )  # set up default map projection

  vn = "z.predicted"
  dev.new();
  spplot( sppoly, vn, main=vn,
    col.regions=p$mypalette,
    at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
    sp.layout=p$coastLayout,
    col="transparent"
  )

  vn = "z.random_strata_nonspatial"
  dev.new();
  spplot( sppoly, vn, main=vn,
    col.regions=p$mypalette,
    at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
    sp.layout=p$coastLayout,
    col="transparent"
  )

  vn = "z.random_strata_spatial"
  dev.new();
  spplot( sppoly, vn, main=vn,
    col.regions=p$mypalette,
    at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
    sp.layout=p$coastLayout,
    col="transparent"
  )


}

# end
