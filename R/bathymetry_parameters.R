bathymetry_parameters = function( p=NULL, project_name=NULL, project_class="default", ... ) {

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  # ---------------------

  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "spdep", "splancs", "GADMTools" ) )
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

  if (!exists("aegis_dimensionality", p)) p$aegis_dimensionality="space"

  p = spatial_parameters( p=p)  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change


  if (project_class=="default") {
    return(p)
  }

  if (project_class=="stmv") {
    p$libs = c( p$libs, project.library ( "stmv" ) )

    if ( !exists("DATA", p) ) p$DATA = 'bathymetry.db( p=p, DS="stmv_inputs" )'
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
