bathymetry_parameters = function( p=NULL, project_name=NULL, project_class="default", ... ) {

  # ---------------------

  if (is.null(p)) p=list()

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "spdep", "splancs", "GADMTools", "INLA" ) )
  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry",  "aegis.polygons", "aegis.coastline") )

  p$project_name = ifelse ( !is.null(project_name), project_name, "bathymetry" )

  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )

  p = parameters_add_without_overwriting( p,
    datadir  = file.path( p$data_root, "data" ),
    modeldir = file.path( p$data_root, "modelled" ),
    variabletomodel = "z",
    spatial_domain = "canada.east.superhighres",
    spatial_domain_subareas = c( "canada.east", "SSE", "snowcrab", "SSE.mpa" ),
    aegis_dimensionality="space"
  )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=F, recursive=T )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=F, recursive=T )


  p = spatial_parameters( p=p )  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change

  if (project_class=="default") return(p)

  if (project_class=="model") {

    p = parameters_add_without_overwriting( p,
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
        transf = function(x) {log10(x + 2500)} ,
        invers = function(x) {10^(x) - 2500}
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

    return(p)
  }

  # ---------------------


  if (project_class=="carstm") {

      # default "production" / best/optimal/etc param list for global analysis in aegis.bathymetry/inst/scripts/02.bathymetry.carstm.R
      p = parameters_add_without_overwriting( p,
        project_name = "bathymetry",
        spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
        variabletomodel ="z",
        carstm_model_label = "production",
        inputdata_spatial_discretization_planar_km = 1,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
        areal_units_resolution_km = 2, # 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
        areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
        areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
        areal_units_overlay = "none",
        carstm_inputs_aggregated = FALSE
      )
      return(p)
  }


  if (project_class=="carstm_force_default") {

    # reset a few project specific params, forcing the use of defaults
    p$data_root = NULL
    p$datadir  = NULL
    p$data_transformation=list( forward=function(x){ x+2500 }, backward=function(x) {x-2500} )
    p$carstm_modelcall = NULL  # defaults to generic
    p$carstm_model_tag = NULL

    p = aegis_parameters( p=p, DS="default" )

    p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry", "carstm"  ) )

    p = parameters_add_without_overwriting( p,
      areal_units_source = "lattice", # "stmv_lattice" to use aegis fields instead of carstm fields ... note variables are not the same
      carstm_modelengine = "inla"  # {model engine}.{label to use to store}
    )


    if ( p$spatial_domain == "SSE" ) {
      p = parameters_add_without_overwriting( p,
        areal_units_overlay = "groundfish_strata", #.. additional polygon layers for subsequent analysis for now ..
        areal_units_resolution_km = 25,  # km dim of lattice ~ 1 hr
        areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
        inputdata_spatial_discretization_planar_km = 0.5  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      )
    }

    if ( p$spatial_domain == "snowcrab" ) {
      p = parameters_add_without_overwriting( p,
        areal_units_overlay = "snowcrab_managementareas", # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
        areal_units_resolution_km = 25 , # km dim of lattice ~ 1 hr
        areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      #   areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
        inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      )
    }


    if ( !exists("carstm_modelcall", p)  ) {

      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )

        p = parameters_add_without_overwriting( p, carstm_model_label = "production" )

        p$carstm_modelcall = paste(
          'inla(
            formula = ', p$variabletomodel, ' ~ 1
              + f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "lognormal",
            data= M,
            control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE, config=TRUE),  # config=TRUE if doing posterior simulations
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),  # extra work to get tails
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            # control.inla = list(cmin = 0 ),
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            control.inla = list(h=1e-4, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            verbose=TRUE
          ) ' )
      }

      if ( grepl("glm", p$carstm_modelengine) ) {
        p = parameters_add_without_overwriting( p, carstm_model_label = "default_glm" )
        p$carstm_modelcall = paste( 'glm( formula = ', p$variabletomodel, '  ~ 1 + AUID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) ' ) # for modelengine='glm'
      }

      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv" ) ) )
        p = parameters_add_without_overwriting( p, carstm_model_label = "default_gam" )
        p$carstm_modelcall = paste( 'gam( formula = ', p$variabletomodel, '  ~ 1 + AUID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) ' )  # for modelengine='gam'
      }

      p = carstm_parameters( p=p )  # fill in anything missing and some checks

    }

    return(p)
  }



  if (project_class=="stmv") {
    p = parameters_add_without_overwriting(
      project_class="stmv",
      data_root = project.datadirectory( "aegis", "bathymetry" ),
      DATA = 'bathymetry_db( p=p, DS="stmv_inputs" )',
      stmv_variables = list(Y="z"),  # required as fft has no formulae
      inputdata_spatial_discretization_planar_km = pres,  #  controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
      spatial_domain = "canada.east.superhighres",
      spatial_domain_subareas = c( "canada.east.highres", "canada.east",  "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
      aegis_dimensionality="space",
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

    p$libs = c( p$libs, project.library ( "stmv" ) )

    if ( !exists("DATA", p) ) p$DATA = 'bathymetry_db( p=p, DS="stmv_inputs" )'
    if ( !exists("stmv_variables", p)) p$stmv_variables = list()
    if ( !exists("LOCS", p$stmv_variables)) p$stmv_variables$LOCS = c("plon", "plat")
    if ( !exists("inputdata_spatial_discretization_planar_km", p) ) p$inputdata_spatial_discretization_planar_km = p$pres / 10 # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)

    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine="fft"  # fft is the perferred approach for bathymetry

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
      if (!exists("stmv_fft_filter", p))  p$stmv_fft_filter = "matern_tapered_modelled" # only act as a low pass filter .. depth has enough data for this. Otherwise, use:
      # p$stmv_fft_filter = "matern" to ~ krige
      if (!exists("stmv_lowpass_phi", p))  p$stmv_lowpass_phi = 0.5 # low pass FFT filter range .. 0.5 seems to be optimal (by visual inspection)
      if (!exists("stmv_lowpass_nu", p))  p$stmv_lowpass_nu = 0.5 # this is exponential covar
    }

    if (p$stmv_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      # timings:
      # 14 hrs on hyperion with 100 knots
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$variabletomodel, ' ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts")') )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      if (!exists("stmv_gam_optimizer", p))  p$stmv_gam_optimizer = c("outer", "bfgs")
    }

    if ( p$stmv_local_modelengine == "bayesx" ) {
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      if (!exists("stmv_local_model_bayesxmethod", p)) p$stmv_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      if (!exists("stmv_local_modelformula", p)) p$stmv_local_modelformula = formula( paste(
        p$variabletomodel, ' ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te")') )   # more detail than "gs" .. "te" is preferred
    }

    if (p$stmv_local_modelengine == "inla" ){
      # old method .. took a month to finish .. results look good but very slow
      ## data range is from -1667 to 5467 m .. 2000 shifts all to positive valued by one order of magnitude
      if (!exists("inla.alpha", p))  p$inla.alpha = 0.5 # bessel function curviness .. ie "nu"
      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$variabletomodel, ' ~ -1 + intercept + f( spatial.field, model=SPDE )' ) ) # SPDE is the spatial covaria0nce model .. defined in stmv___inla
      if (!exists("stmv.posterior.extract", p)) p$stmv.posterior.extract = function(s, rnm) {
        # rnm are the rownames that will contain info about the indices ..
        # optimally the grep search should only be done once but doing so would
        # make it difficult to implement in a simple structure/manner ...
        # the overhead is minimal relative to the speed of modelling and posterior sampling
        i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
        i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
        return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
      }

      p = aegis_parameters( p=p, DS="stmv" )

    }
    return(p)
  }


}
