
bathymetry_parameters = function( p=list(), project_name="bathymetry", project_class="default", ... ) {

  p = parameters_add( p, list(...) ) # add passed args to parameter list, priority to args

  # ---------------------

  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "spdep", "splancs", "GADMTools", "INLA" ) )
  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry",  "aegis.polygons", "aegis.coastline" ) )

  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )

  p = parameters_add_without_overwriting( p,
    variabletomodel = "z",
    spatial_domain = "canada.east.superhighres",
    spatial_domain_subareas = c( "canada.east",  "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    aegis_dimensionality="space"
  )
  p = spatial_parameters( p=p )  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change

  p = parameters_add_without_overwriting( p,
      inputdata_spatial_discretization_planar_km = p$pres/10  #  controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
  )

  # ---------------------

  if (project_class=="default") return(p)  # minimal specifications

  # ---------------------


  if (project_class=="carstm") {
    # simple run of carstm. There are two types:
    #   one global, run directly from  polygons defined in aegis.bathymetry/inst/scripts/99.bathymetry.carstm.R
    #   and one that is called secondarily specific to a local project's polygons (eg. snow crab)
    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )

    # defaults in case not provided ...
    p = parameters_add_without_overwriting( p,
      data_transformation=list( forward=function(x){ x+2500 }, backward=function(x) {x-2500} ),
      areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_resolution_km = 5, # default in case not provided ... 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      areal_units_overlay = "none",
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "default",
      carstm_inputs_aggregated = FALSE
    )


    if ( !exists("carstm_model_call", p)  ) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$carstm_model_call = paste(
          'inla(
            formula = ', p$variabletomodel, ' ~ 1
              + f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "lognormal",
            data= M,
            control.compute=list(dic=TRUE, waic=TRUE, config=TRUE),  # config=TRUE if doing posterior simulations
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            control.inla = list(h=1e-4, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            verbose=TRUE
          ) ' )
      }
    }

    p = carstm_parameters( p=p )  # fill in anything missing with defaults and do some checks

    return(p)

  }


  # ---------------------


  if (project_class=="stmv") {
    p = parameters_add_without_overwriting( p,
      project_class="stmv",
      DATA = 'bathymetry_db( p=p, DS="stmv_inputs" )',
      stmv_variables = list(Y="z"),  # required as fft has no formulae
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="fft",
      stmv_fft_filter = "matern tapered lowpass modelled fast_predictions", #  act as a low pass filter first before matern with taper .. depth has enough data for this. Otherwise, use:
      stmv_lowpass_nu = 0.5, # exp
      stmv_lowpass_phi = stmv::matern_distance2phi( distance=0.1, nu=0.5, cor=0.1 ),
      stmv_autocorrelation_fft_taper = 0.9,  # benchmark from which to taper
      stmv_autocorrelation_localrange = 0.1,  # for output to stats
      stmv_autocorrelation_basis_interpolation = c(0.25, 0.1, 0.05, 0.01),
      stmv_variogram_method = "fft",
      stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
      stmv_Y_transform =list(
        transf = function(x) {log10(x + 2500)} ,
        invers = function(x) {10^(x) - 2500}
      ), # data range is from -1667 to 5467 m: make all positive valued
      stmv_rsquared_threshold = 0.01, # lower threshold  .. ignore
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_distance_prediction_limits =c( 3, 25 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
      stmv_distance_scale = c( 5, 10, 20, 25, 40, 80, 150, 200), # km ... approx guesses of 95% AC range
      stmv_distance_basis_interpolation = c(  2.5 , 5, 10, 15, 20, 40, 80, 150, 200 ) , # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_force_complete_method = "linear_interp"
    )
    p$stmv_variables = parameters_add_without_overwriting( p$stmv_variables, LOCS = c("plon", "plat") )
    p$libs = c( p$libs, project.library ( "stmv" ) )


    # tweaked override the defaults of aegis_parameters( p=p, DS="stmv") :
    if ( p$stmv_local_modelengine %in% c("krige" )) {
      # nothing to do  .. this is faster than "gstat" .. do not use for bathymetry as it is oversmoothed
    }

    if ( p$stmv_local_modelengine =="gaussianprocess2Dt" ) {
      # too slow to use right now
      if (!exists("fields.cov.function", p)) p$fields.cov.function = "stationary.taper.cov"  # Wendland tapering; also "stationary.cov"  #
      if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS
      if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
        if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)
        p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance
      }
    }

    if (p$stmv_local_modelengine == "fft") {
      # ~ 3.25 days hr with 68, 3 Ghz cpus on beowulf using fft method, bigmemory-filebacked jc: 2016
      # ~ 14 hrs with 8, 3.2 Ghz cpus on thoth; 1 GB per process and a total of 6 GB usage;  method RAM based jc: 2016
      # 12 hrs to complete stage 1 on hyperion
      # ~ 5.5 hr on hyperion
      # definitely a cleaner (not overly smoothed) image than a GAM
      # NOTE that  p$stmv_lowpass_phi and  p$stmv_lowpass_nu are very critical choices
      p = parameters_add_without_overwriting( p,
        stmv_fft_filter = "matern_tapered_modelled", # only act as a low pass filter .. depth has enough data for this. Otherwise, use:
        stmv_lowpass_phi = 0.5, # low pass FFT filter range .. 0.5 seems to be optimal (by visual inspection)
        stmv_lowpass_nu = 0.5 # this is exponential covar
      )
    }

    if (p$stmv_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      # timings:
      # 14 hrs on hyperion with 100 knots
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      p = parameters_add_without_overwriting( p,
        stmv_local_modelformula = formula( paste(
          p$variabletomodel, ' ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts")') ),
        stmv_local_model_distanceweighted = TRUE,
        stmv_gam_optimizer = c("outer", "bfgs")
      )
    }

    if ( p$stmv_local_modelengine == "bayesx" ) {
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      p = parameters_add_without_overwriting( p,
        stmv_local_model_bayesxmethod="MCMC",  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
        stmv_local_model_distanceweighted = TRUE,
        stmv_local_modelformula = formula( paste(
          p$variabletomodel, ' ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te")') )   # more detail than "gs" .. "te" is preferred
      )
    }

    if (p$stmv_local_modelengine == "inla" ){
      # old method .. took a month to finish .. results look good but very slow
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
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
    p = aegis_parameters( p=p, DS="stmv" )
    return(p)
  }


  # ---------------------


  if (project_class=="production") {

    # "best" interpolation method .. in this case, due to data density ... stmv_carstm hybrid

    p = parameters_add_without_overwriting( p,
      project_class="stmv",
      DATA = 'bathymetry_db( p=p, DS="stmv_inputs_highres" )',  # _highres
      stmv_variables = list(Y="z"),  # required as fft has no formulae
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="carstm",
      stmv_local_covariates_carstm = "",  # only model covariates
      stmv_local_all_carstm = "",  # ignoring aui
      stmv_local_modelcall = paste(
        'inla(
          formula = z ~ 1
            + f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
          family = "normal",
          data= dat,
          control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE),  # config=TRUE if doing posterior simulations
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list(compute=FALSE, link=1 ),
          control.fixed=H$fixed,  # priors for fixed effects, generic is ok
          verbose=TRUE
        ) '
      ),   # NOTE:: this is a local model call
      stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
      stmv_Y_transform =list(
        transf = function(x) {log10(x + 2500)} ,
        invers = function(x) {10^(x) - 2500}
      ), # data range is from -1667 to 5467 m: make all positive valued
      stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    #  stmv_interpolation_basis_distance = 5,   # fixed distance 2 x statsgrid
      stmv_distance_prediction_limits =c( 2, 25 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
      stmv_interpolation_basis_distance_choices = c(2, 4, 8),
      stmv_nmin = 10, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 10000, # no real upper bound.. just speed /RAM
      stmv_force_complete_method = "linear_interp",
      stmv_runmode = list(
        carstm = rep("localhost", 2),
        globalmodel = FALSE,
        # restart_load = "interpolate_correlation_basis_0.01" ,  # only needed if this is restarting from some saved instance
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )  # ncpus for each runmode
    )

    return(p)
  }



}
