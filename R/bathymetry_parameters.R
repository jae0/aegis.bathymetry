bathymetry_parameters = function( p=NULL, project_name=NULL, project_class="default", DS="none", ... ) {

  # ---------------------
  p = parameters_control(p, list(...), control="add") # add passed args to parameter list, priority to args

  if (DS=="default") {

    p = bathymetry_parameters(
      p=p,
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
    return(p)
  }

  # ---------------------

  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "spdep", "splancs", "GADMTools", "INLA" ) )
  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry",  "aegis.polygons", "aegis.coastline") )

  p$project_name = ifelse ( !is.null(project_name), project_name, "bathymetry" )

  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=F, recursive=T )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=F, recursive=T )

  if ( !exists("variabletomodel", p)) p$variabletomodel = "z"

  if (!exists("spatial_domain", p) ) p$spatial_domain = "canada.east.superhighres"
  if (!exists("spatial_domain_subareas", p)) p$spatial_domain_subareas = c( "canada.east", "SSE", "snowcrab", "SSE.mpa" )

  p$aegis_dimensionality="space"

  p = spatial_parameters( p=p)  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change


  if (project_class=="default") return(p)

  if (project_class=="carstm")  return(p)

  if (project_class=="stmv") {
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
