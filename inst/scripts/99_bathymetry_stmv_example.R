
# model testing example


# testing "stmv" on a coarse grid and lower res data

# model testing
p0 = aegis::spatial_parameters(
  spatial_domain="bathymetry_example",
  aegis_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2
)
# or:
p0 = stmv_test_data( "aegis.test.parameters")


input = stmv::stmv_test_data( datasource="aegis.space", p=p0)
input = sf::st_as_sf( input, coords=c("lon","lat"), crs=st_crs(projection_proj4string("lonlat_wgs84")) )
input = sf::st_transform( input, crs=st_crs(p0$aegis_proj4string_planar_km) )
input = as.data.frame( cbind( input$z, st_coordinates(input) ) )
names(input) = c("z", "plon", "plat")
input = input[ which(is.finite(input$z)), ]
output = list( LOCS = spatial_grid(p0) )
DATA = list( input = input,  output = output )

input = output = NULL
gc()


p = bathymetry_parameters(
  p=p0,  # start with spatial settings of input data
  project_class="stmv",
  stmv_model_label="testing_statsgrid10km",
  data_root = file.path(tempdir(), "bathymetry_example"),
  DATA = DATA,
  spatial_domain = p0$spatial_domain,
  spatial_domain_subareas =NULL,
  inputdata_spatial_discretization_planar_km = p0$pres,  # pres = 0.5
  aegis_dimensionality="space",
  stmv_variables = list(Y="z"),  # required as fft has no formulae
  stmv_global_modelengine = "none",  # too much data to use glm as an entry into link space ... use a direct transformation
  stmv_local_modelengine="fft",
  stmv_fft_filter = "matern tapered lowpass modelled fast_predictions", #  matern with taper, fast predictions are sufficient as data density is high
  stmv_distance_statsgrid = 10, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_prediction_limits =c( 3, 25 ), # range of permissible predictions km (i.e 1/2 stats grid to upper
  stmv_runmode = list(
    scale = rep("localhost", 1),
    interpolate = list(
      c1 = rep("localhost", 1),
      c2 = rep("localhost", 1),
      c3 = rep("localhost", 1),
      c4 = rep("localhost", 1),
      c5 = rep("localhost", 1)
    ),
    interpolate_predictions = list(
      c1 = rep("localhost", 1),
      c2 = rep("localhost", 1),
      c3 = rep("localhost", 1),
      c4 = rep("localhost", 1),
      c5 = rep("localhost", 1),
      c6 = rep("localhost", 1),
      c7 = rep("localhost", 1) 
    ),
    globalmodel = FALSE,
    save_intermediate_results = TRUE,
    save_completed_data = TRUE
  ) 
)


if (0) {
  # to force parallel mode
   stmv_runmode = list(
    scale = rep("localhost", scale_ncpus),
    interpolate = list(
      c1 = rep("localhost", interpolate_ncpus),  # ncpus for each runmode
      c2 = rep("localhost", interpolate_ncpus),  # ncpus for each runmode
      c3 = rep("localhost", max(1, interpolate_ncpus-1)),
      c4 = rep("localhost", max(1, interpolate_ncpus-1)),
      c5 = rep("localhost", max(1, interpolate_ncpus-2))
    ),
    interpolate_predictions = list(
      c1 = rep("localhost", max(1, interpolate_ncpus-1)),  # ncpus for each runmode
      c2 = rep("localhost", max(1, interpolate_ncpus-1)),  # ncpus for each runmode
      c3 = rep("localhost", max(1, interpolate_ncpus-2)),
      c4 = rep("localhost", max(1, interpolate_ncpus-3)),
      c5 = rep("localhost", max(1, interpolate_ncpus-4)),
      c6 = rep("localhost", max(1, interpolate_ncpus-4)),
      c7 = rep("localhost", max(1, interpolate_ncpus-5))
    ),
    globalmodel = FALSE,
    # restart_load = "interpolate_correlation_basis_0.01" ,  # only needed if this is restarting from some saved instance
    save_intermediate_results = TRUE,
    save_completed_data = TRUE

  )  # ncpus for each runmode
}


# quick view
predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
statistics  = stmv_db( p=p, DS="stmv.stats" )
locations   = spatial_grid( p )

# comparison
dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )


statsvars = dimnames(statistics)[[2]]

# statsvars = c("sdTotal", "ndata", "fixed_mean", "fixed_sd", "dic", "dic_p_eff",
#   "waic", "waic_p_eff", "mlik", "Expected_number_of_parameters",
#   "Stdev_of_the_number_of_parameters", "Number_of_equivalent_replicates",
#   "Precision_for_the_Gaussian_observations", "Precision_for_aui",
#   "Phi_for_aui", "Precision_for_the_Gaussian_observations_sd", "Precision_for_aui_sd", "Phi_for_aui_sd"
# )

# statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
dev.new(); levelplot( predictions[] ~ locations[,1] + locations[,2], aspect="iso" )
dev.new(); levelplot( statistics[,match("Phi_for_auid", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdTotal", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("rsquared", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange
