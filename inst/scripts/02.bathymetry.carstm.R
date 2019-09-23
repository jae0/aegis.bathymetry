
# ------------------------------------------------
# Areal unit modelling of bathymetry


reset_input_data = FALSE
if (0) reset_input_data = TRUE # choose this if we are redoing input data "views"

# subproject defines spatial bounds and lattive areal_units_overlays
spatial_domain = "snowcrab"
subproject = "snowcrab"

if (0) {
  # alternatively:
  spatial_domain="SSE"
  subproject = "groundfish"
}

# --------------------------------
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p = aegis.bathymetry::bathymetry_parameters(
  project_class = "carstm", # defines which parameter set to load
  id = paste("bathymetry", subproject, sep="_"),
  inputdata_spatial_discretization_planar_km = 0.05,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  spatial_domain = spatial_domain,  # defines spatial area
  aegis_proj4string_planar_km = projection_proj4string("lonlat_wgs84"),  # coord system of the data lon/lats
  areal_units_resolution_km = 100, # km dim of lattice
  areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  # areal_units_proj4string_planar_km = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  areal_units_strata_type = "lattice", # "aegis_lattice" to use ageis fields instead of carstm fields ... note variables are not the same
  areal_units_overlay = subproject, # additional polygon layers for subsequent analysis such as management area: "snowcrab" or "groundfish"  # for now ..
  areal_units_constraint="none", # set[, c("lon", "lat")],  # to limit to sppoly to only those with data that fall into them
  carstm_modelengine = "inla",  # gam is also supported for now
  carstm_family = "lognormal",
  carstm_formula = formula( paste(
    'z ~ 1 +
      + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(iid_error, model="iid", hyper=H$iid)'
  )),
  constant_offset = 2500,
  libs = RLibrary ( "sp", "spdep", "rgeos", "spatialreg", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry", "carstm" )
)

p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
p$mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")
p = c(p, aegis.coastline::coastline_layout( p=p, redo=reset_input_data ) )  # set up default map projection



# --------------------------------
# ensure if polys exist and create if required
# for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygons_managementarea( species="snowcrab", au))

if (0) {
  # force creating of input data for modelling
  sppoly = bathymetry_carstm( p=p, DS="areal_units", redo=TRUE )  # will redo if not found
  plot(sppoly)

  M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  str(M)
}

# run model and obtain predictions
sppoly = bathymetry_carstm( p=p, DS="carstm_modelled", redo=TRUE ) # extract predictions in sppoly

if (0) {
  sppoly = bathymetry_carstm( p=p, DS="carstm_modelled" ) # to load currently saved sppoly
  fit =  bathymetry_carstm( p=p, DS="carstm_modelled" )  # extract currently saved model fit
}

vn = "z.predicted"
dev.new();
spplot( sppoly, vn, main=vn,
  col.regions=p$mypalette,
  at=interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile"),
  sp.layout=p$coastLayout,
  col="transparent"
)


# end
