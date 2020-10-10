


### map area unit problem example:

  maup = map_area_unit_problem( just_return_results=TRUE )  #default is bathymetry data


### stmv/carstm hybrid:

p = aegis.bathymetry::bathymetry_parameters(
  project_class="stmv",
  data_root = project.datadirectory( "aegis", "bathymetry" ),
  DATA = 'bathymetry_db( p=p, DS="stmv_inputs_highres" )',  # _highres
  stmv_variables = list(Y="z"),  # required as fft has no formulae
  spatial_domain = "canada.east.superhighres",
  aegis_dimensionality="space",
  stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
  stmv_local_modelengine="carstm",
  stmv_local_covariates_carstm = "",  # only model covariates
  stmv_local_all_carstm = "",  # ignoring aui
  stmv_local_modelcall = paste(
          'inla(
            formula = z ~ 1
              + f(aui, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "normal",
            data= dat,
            control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE),  # config=TRUE if doing posterior simulations
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            verbose=TRUE
          ) '
    ),
  stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
  stmv_Y_transform =list(
    transf = function(x) {log10(x + 3000)} ,
    invers = function(x) {10^(x) - 3000}
  ), # data range is from -1667 to 5467 m: make all positive valued
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
#  stmv_interpolation_basis_distance = 5,   # fixed distance 2 x statsgrid
  stmv_interpolation_basis_distance_choices = c(5),
  stmv_distance_prediction_limits =c( 3, 25 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
  global_sppoly = NULL, # force local lattice grid of pres ... bathymetry has enough data for this
  stmv_nmin = 100, # min number of data points req before attempting to model in a localized space
  stmv_nmax = 5000, # no real upper bound.. just speed /RAM
  stmv_force_complete_method = "linear_interp",
  stmv_runmode = list(
    carstm = rep("localhost", 2),
    globalmodel = FALSE,
    # restart_load = "interpolate_correlation_basis_0.01" ,  # only needed if this is restarting from some saved instance
    save_intermediate_results = TRUE,
    save_completed_data = TRUE
  )  # ncpus for each runmode
)

if (0) {  # model testing
  bathymetry_db( p=p, DS="stmv_inputs_highres_redo" )  # recreate fields for .. requires 60GB+
}


stmv( p=p )  # This will take from 40-70 hrs, depending upon system





if (0) {  # model testing
  # if resetting data for input to stmv run this or if altering discretization resolution
  bathymetry_db( p=p, DS="stmv_inputs_redo" )  # recreate fields for .. requires 60GB+
  # p$restart_load = paste("interpolate_correlation_basis_", p$stmv_autocorrelation_basis_interpolation[length(p$stmv_autocorrelation_basis_interpolation)], sep="")  # to choose the last save
}

stmv( p=p )  # This will take from 40-70 hrs, depending upon system


# quick view
  predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
  statistics  = stmv_db( p=p, DS="stmv.stats" )
  locations   = spatial_grid( p )

  # comparisons
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

  # stats
  # p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
  dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )

  dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
  dev.new(); levelplot( statistics[,match("sdSpatial", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
  dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange

# water only
o = which( predictions>0 & predictions <1000)
#levelplot( log( statistics[o,match("sdTotal", p$statsvars)] ) ~ locations[o,1] + locations[o,2], aspect="iso" ) #sd total
levelplot( log(predictions[o]) ~ locations[o,1] + locations[o,2], aspect="iso" )



# bring together stats and predictions and any other required computations: slope and curvature
# and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# .. still uses about 30-40 GB as the base layer is "superhighres" ..
# if parallelizing .. use different servers than local nodes


bathymetry_db( p=p, DS="complete.redo", spatial_domain_subareas = c( "canada.east.highres", "canada.east",  "SSE", "SSE.mpa" , "snowcrab") ) # finalise at diff resolutions 15 min ..


# ----





### WARNING::  running carstm directly... < 5 km resokution becomes rather slow ...

## SLOW :: about 1 hr for each config x 25 configs .. ie., 24 hrs .. consider removing configs if no need for posterior samples

# -- might need to run in a shell if number of polygons are large:  ulimit -s 16384


# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
# --- look inside "parameters_production" and define alternates based upon it
  p = aegis.bathymetry::bathymetry_carstm( DS = "parameters_production" )
    # DS = "parameters_production"; areal_units_resolution_km=5 ... takes 79 Hrs!


  if (0) {

    x = maup$resolution
    # x = log(  maup$n / maup$resolution^2  )  # data density
    yrange = range( maup$min, maup$max )
    plot(mean ~ x, maup, pch=20, ylim=yrange)
    lines( median ~ x, maup, col="green", lwd=1)
    lines(min ~ x, maup, col="red", lwd=4)
    lines(max ~ x, maup, col="blue", lwd=4)
    lines( I(mean+sd) ~ x, maup, col="gray", lwd=2)
    lines( I(mean-sd) ~ x, maup, col="gray", lwd=2)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

    plot( sd~ x, maup)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

  }

# example sequence to force creating of input data for modelling
  sppoly = areal_units( p=p, redo=TRUE ); plot(sppoly) # or: spplot( sppoly, "AUID", main="AUID", sp.layout=p$coastLayout )
  M = bathymetry_db( p=p, DS="aggregated_data" , redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  str(M)
  M = NULL; gc()

# run the model ... about 24 hrs
  fit = carstm_model( p=p, M='bathymetry_carstm( p=p, DS="carstm_inputs" )' ) # run model and obtain predictions

# loading saved results
  fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  res = carstm_summary( p=p ) # to load currently saved sppoly

  plot(fit)
  plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  s = summary(fit)
  s$dic$dic # 1404489
  s$dic$p.eff # 151412

  # maps of some of the results
  vn = paste(p$variabletomodel, "predicted", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

# end
