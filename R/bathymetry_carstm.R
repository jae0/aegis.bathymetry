
  bathymetry_carstm = function( p=NULL, DS="parameters", redo=FALSE, ... ) {

    #\\ Note inverted convention: depths are positive valued
    #\\ i.e., negative valued for above sea level and positive valued for below sea level

  # ------------------


  if (DS=="parameters_production") {

      # copy of param list for global analysis in aegis.bathymetry/inst/scripts/02.bathymetry.carstm.R
      p = aegis.bathymetry::bathymetry_carstm(
        DS = "parameters",
        project_name = "bathymetry",
        spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
        variabletomodel ="z",
        carstm_model_label = "production",
        inputdata_spatial_discretization_planar_km = 1,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
        areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
        areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
        areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
        areal_units_overlay = "none"
      )
      return(p)
  }


  if (DS=="parameters") {

    if ( is.null(p)) {

      p = aegis.bathymetry::bathymetry_parameters(...)

    } else {

      # if p is passed, assume it is a secondary call ... overwrite nonrelevent params to force use of project defaults
      p$data_root = NULL
      p$datadir  = NULL
      p$data_transformation=list( forward=function(x){ x+2500 }, backward=function(x) {x-2500} )
      p$carstm_modelcall = NULL  # defaults to generic
      p$carstm_model_tag = NULL

      p = aegis.bathymetry::bathymetry_parameters(p=p, ...)
    }

    p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry", "carstm"  ) )

    if ( !exists("areal_units_source", p)) p$areal_units_source = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( p$spatial_domain == "SSE" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "groundfish_strata" #.. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 0.5  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( p$spatial_domain == "snowcrab" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    }



    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)  ) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )

        if ( !exists("carstm_model_label", p)) p$carstm_model_label = "production"

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
        if ( !exists("carstm_model_label", p))  p$carstm_model_label = "default_glm"
        p$carstm_modelcall = paste( 'glm( formula = ', p$variabletomodel, '  ~ 1 + AUID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) ' ) # for modelengine='glm'
      }
      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv" ) ) )
        if ( !exists("carstm_model_label", p)) p$carstm_model_label = "default_gam"
        p$carstm_modelcall = paste( 'gam( formula = ', p$variabletomodel, '  ~ 1 + AUID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ], family=gaussian(link="log")  ) ' )  # for modelengine='gam'
      }

      p = carstm_parameters( p=p )  # fill in anything missing and some checks

      p$carstm_inputs_aggregated = TRUE

    }

    return(p)
  }


  # --------------------


  if ( DS=="carstm_inputs") {

    # prediction surface
    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))
    sppoly = areal_units( p=p )  # will redo if not found
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    if (p$carstm_inputs_aggregated) {
      fn = carstm_filenames( p=p, projectname="bathymetry", projecttype="carstm_inputs", areal_units_fn=areal_units_fn )
    } else {
      fn = paste( "bathymetry", "carstm_inputs", areal_units_fn, "rawdata", "rdata", sep=".")
    }

    fn = file.path( p$modeldir, fn)

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }

    # reduce size
    if (p$carstm_inputs_aggregated) {
      M = bathymetry.db ( p=p, DS="aggregated_data"   )  # 16 GB in RAM just to store!
      names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel

    } else {
      M = bathymetry.db ( p=p, DS="z.lonlat.rawdata"  )  # 16 GB in RAM just to store!
      names(M)[which(names(M)=="z") ] = p$variabletomodel
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")
          # p$quantile_bounds_data = c(0.0005, 0.9995)
      if (exists("quantile_bounds_data", p)) {
        TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds_data, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
        keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
        # this was -1.7, 21.8 in 2015
      }
    }

    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    M = lonlat2planar(M, p$aegis_proj4string_planar_km)  # should not be required but to make sure
    # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")

    if( exists("spatial_domain", p)) M = geo_subset( spatial_domain=p$spatial_domain, Z=M ) # need to be careful with extrapolation ...  filter depths

    M$AUID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(!is.na(M$AUID)),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    M$tag = "observations"

    sppoly_df = as.data.frame(sppoly)
    sppoly_df[, p$variabletomodel] = NA
    sppoly_df$AUID = as.character( sppoly_df$AUID )
    sppoly_df$tag ="predictions"

    vn = c("z", "tag", "AUID")

    M = rbind( M[, vn], sppoly_df[, vn] )
    sppoly_df = NULL

    M$auid  = as.numeric( factor(M$AUID) )

    save( M, file=fn, compress=TRUE )
    return( M )
  }

}  # end
