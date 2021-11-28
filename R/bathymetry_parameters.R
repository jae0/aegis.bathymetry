
bathymetry_parameters = function( p=list(), project_name="bathymetry", project_class="core", ... ) {

  p = parameters_add( p, list(...) ) # add passed args to parameter list, priority to args

  # ---------------------

  # create/update library list
  # p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
  #   "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "spdep", "splancs", "GADMTools", "INLA" ) )
  p$libs = c( p$libs, RLibrary ( "colorspace",  "lubridate",  "lattice",
    "parallel", "sf", "GADMTools", "INLA", "data.table" ) )

  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry",  "aegis.polygons", "aegis.coastline", "aegis.survey" ) )


  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )  # all unprocessed inputs (and simple manipulations)
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )  # all outputs

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )


  p = parameters_add_without_overwriting( p,
    variabletomodel = "z",  
    spatial_domain = "canada.east.superhighres",
    spatial_domain_subareas = c( "canada.east.highres", "canada.east", "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    aegis_dimensionality="space"
  )

  p$quantile_bounds =c(0, 0.95) # trim upper bounds

  p = spatial_parameters( p=p )  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change

  p$discretization = discretizations(p=p$discretization)  # key for discretization levels

  p = parameters_add_without_overwriting( p, inputdata_spatial_discretization_planar_km = 0.5 ) #  p$pres/2 is a bit too slow ..; controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)


  # ---------------------

  if (project_class=="core") {
    p$project_class = "core"
    return(p)  # minimal specifications
  }

  # ---------------------


  if (project_class %in% c("carstm") ) {
    # simple run of carstm. There are two types:
    #   one global, run directly from  polygons defined in aegis.bathymetry/inst/scripts/99.bathymetry.carstm.R
    #   and one that is called secondarily specific to a local project's polygons (eg. snow crab)
    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )
    p$project_class = "carstm"

    # defaults in case not provided ...
    p = parameters_add_without_overwriting( p,
      data_transformation=list( forward=function(x){ x+2500 }, backward=function(x) {x-2500} ),
      areal_units_xydata = "bathymetry_db(p=p, DS='areal_units_input')",
      areal_units_type = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_resolution_km = 5, # default in case not provided ... 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
      areal_units_proj4string_planar_km = p$aegis_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      areal_units_overlay = "none",
      areal_units_timeperiod = "none",
      areal_units_constraint_ntarget = 500,
      areal_units_constraint_nmin = 30 ,
      tus="none",
      fraction_todrop = 1/5,
      fraction_cv = 1.0,
      fraction_good_bad = 0.9,
      nAU_min = 30,
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "default",  
      carstm_inputs_prefilter = "aggregated",
      carstm_inputs_prefilter_n = 100  # used only if "sampled" 10
    )

    if ( grepl("inla", p$carstm_modelengine) ) {
      if ( !exists("formula", p)  ) {
        p$formula = as.formula( paste(
          p$variabletomodel, ' ~ 1',
             ' + f(space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, hyper=H$bym2) '
          ))
      }
      if ( !exists("family", p)  )  p$family = "lognormal"
    }

  
    return(p)

  }


  # ---------------------


  if (project_class %in% c("stmv", "default") ) {

    p$libs = c( p$libs, project.library ( "stmv" ) )
    p$project_class = "stmv"

    p = parameters_add_without_overwriting( p,
      DATA = 'bathymetry_db( p=p, DS="stmv_inputs" )',
      stmv_model_label="default",
      stmv_variables = list(Y="z"),  # required as fft has no formulae
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="fft",
      stmv_variogram_method = "fft",
      stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
      stmv_Y_transform =list(
        transf = function(x) {log10(x + 2500)} ,
        invers = function(x) {10^(x) - 2500}
      ), # data range is from -1667 to 5467 m: make all positive valued
      stmv_rsquared_threshold = 0.01, # lower threshold  .. ignore
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_force_complete_method = "linear_interp"
    )



    p = parameters_add_without_overwriting( p,
      stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1/2, 5 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
      stmv_distance_scale = p$stmv_distance_statsgrid * c( 1, 2, 3, 4, 5, 10, 15), # km ... approx guesses of 95% AC range
      stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 5, 10, 15 ),  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 8) # finalizing preds using linear interpolation
    )


    # tweaked override the defaults of aegis_parameters( p=p, DS="stmv") :
    if ( p$stmv_local_modelengine %in% c("krige" )) {
      # nothing to do  .. this is faster than "gstat" .. do not use for bathymetry as it is oversmoothed
    }

    if ( p$stmv_local_modelengine =="gaussianprocess2Dt" ) {
      p$libs = unique( c( p$libs, RLibrary ("fields")) )
      # too slow to use right now
      if (!exists("fields.cov.function", p)) p$fields.cov.function = "stationary.taper.cov"  # Wendland tapering; also "stationary.cov"  #
      if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS
      if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
        if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)
        p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance
      }
    }

    if (p$stmv_local_modelengine == "fft") {
      nu = 0.5  # exponential smoothing
      ac_local = 0.1  # ac at which to designate "effective range"
      p = parameters_add_without_overwriting( p,
        stmv_fft_filter = "matern tapered lowpass modelled fast_predictions", #  act as a low pass filter first before matern with taper
        stmv_autocorrelation_fft_taper = 0.9,  # benchmark from which to taper
        stmv_autocorrelation_localrange = ac_local,  # for output to stats
        stmv_autocorrelation_interpolation = c(0.25, 0.1, 0.05, 0.01),
        stmv_lowpass_nu = nu, # exp
        stmv_lowpass_phi = stmv::matern_distance2phi( distance=p$pres, nu=nu, cor=ac_local )
      )
    }

    if (p$stmv_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      # timings:
      # 14 hrs on hyperion with 100 knots
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      p$libs = unique( c( p$libs, RLibrary ("mgcv")))
      p = parameters_add_without_overwriting( p,
        stmv_local_modelformula = formula( paste(
          p$variabletomodel, ' ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts")') ),
        stmv_local_model_distanceweighted = TRUE,
        stmv_gam_optimizer = c("outer", "bfgs")
      )
    }

    if ( p$stmv_local_modelengine == "bayesx" ) {
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      p$libs = unique( c( p$libs, RLibrary ("bayesx")) )
      p = parameters_add_without_overwriting( p,
        stmv_local_model_bayesxmethod="MCMC",  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
        stmv_local_model_distanceweighted = TRUE,
        stmv_local_modelformula = formula( paste(
          p$variabletomodel, ' ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te")') )   # more detail than "gs" .. "te" is preferred
      )
    }

    if (p$stmv_local_modelengine == "inla" ) {
      # old method .. took a month to finish .. results look good but very slow
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      p$libs = unique( c( p$libs, RLibrary ("INLA")) )

      p = parameters_add_without_overwriting( p,
        inla.alpha = 0.5, # bessel function curviness .. ie "nu"
        stmv_local_modelformula = formula( paste(
          p$variabletomodel, ' ~ -1 + intercept + f( spatial.field, model=SPDE )' ) ), # SPDE is the spatial covaria0nce model .. defined in stmv___inla
        stmv.posterior.extract = function(s, rnm) {
          # rnm are the rownames that will contain info about the indices ..
          # optimally the grep search should only be done once but doing so would
          # make it difficult to implement in a simple structure/manner ...
          # the overhead is minimal relative to the speed of modelling and posterior sampling
          i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
          i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
          return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
        }
      )
    }

    # default to serial mode
    p = parameters_add_without_overwriting( p,
      stmv_runmode = list(
        globalmodel = FALSE,
        scale = list(
          cor_0.25 = rep("localhost", 1),
          cor_0.1  = rep("localhost", 1),
          cor_0.05 = rep("localhost", 1),
          cor_0.01 = rep("localhost", 1)
        ),
        interpolate_correlation_basis = list(
          cor_0.25 = rep("localhost", 1),
          cor_0.1  = rep("localhost", 1),
          cor_0.05 = rep("localhost", 1),
          cor_0.01 = rep("localhost", 1)
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
        save_intermediate_results = TRUE,
        save_completed_data = TRUE # just a dummy variable with the correct name
      )
    )

    p = aegis_parameters( p=p, DS="stmv" )

    return(p)
  }


  # ---------------------


  if (project_class %in% c("hybrid")  ) {
    p$project_class = "hybrid"

    p = parameters_add_without_overwriting( p,
      DATA = 'bathymetry_db( p=p, DS="stmv_inputs_highres" )',  # _highres
      stmv_variables = list(Y="z"),  # required as fft has no formulae
      stmv_model_label="default",
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="carstm",
      stmv_local_covariates_carstm = "",  # only model covariates
      stmv_local_all_carstm = "",  # ignoring au
      stmv_local_modelcall = paste(
        'inla(
          formula = z ~ 1
            + f(space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE),
          family = "gaussian",
          data = dat,
          inla.mode="experimental",
          control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE, return.marginals.predictor=TRUE),  # config=TRUE if doing posterior simulations
          control.predictor=list(compute=FALSE, link=1 ),
          verbose=FALSE
        ) '
      ),   # NOTE:: this is a local model call
      stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
      stmv_Y_transform =list(
        transf = function(x) {log10(x + 2500)} ,
        invers = function(x) {10^(x) - 2500}
      ), # data range is from -1667 to 5467 m: make all positive valued
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_au_distance_reference = "none", # additional filters upon polygons relative to windowsize: "centroid", "inside_or_touches_boundary", completely_inside_boundary"
      stmv_au_buffer_links = 0, # number of additional neighbours to extend beyond initial solution
#      pres = 1  # this governs resolution of lattice predictions
      stmv_nmin = 100, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000 # no real upper bound.. just speed /RAM
    )


    p = parameters_add_without_overwriting( p,
      stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1, 2 ), # range of permissible predictions km (i.e  stats grid to upper limit based upon data density)
      stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2 ),  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2) # finalizing preds using linear interpolation
    )


    p = parameters_add_without_overwriting( p,
      stmv_runmode = list(
        carstm = list(
          car1 = rep("localhost", 1),
          car2 = rep("localhost", 1),
          car3 = rep("localhost", 1),
          car4 = rep("localhost", 1)
        ),
        globalmodel = FALSE,
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )
    )

    p = aegis_parameters( p=p, DS="stmv" )  # get defaults

    if ( p$inputdata_spatial_discretization_planar_km >= p$pres ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$pres " )
    }
    message ("p$stmv_distance_statsgrid: ", p$stmv_distance_statsgrid)

    return(p)
  }

}
