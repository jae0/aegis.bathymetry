
  bathymetry_carstm = function( p=NULL, DS="aggregated_data", sppoly=NULL, id=NULL, redo=FALSE, map=TRUE, ... ) {

    #\\ Note inverted convention: depths are positive valued
    #\\ i.e., negative valued for above sea level and positive valued for below sea level
    if ( is.null(p)) p = bathymetry_parameters(...)

    if ( !exists("project.name", p)) p$project.name = "bathymetry"
    if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project.name )
    if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
    if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

    if (is.null(id)) id = paste( p$resolution, p$spatial.domain, p$overlay, sep="_" )

    if ( DS=="aggregated_data") {

      fn = file.path( p$modeldir, paste( "bathymetry", "aggregated_data", id, "rdata", sep=".") )
      if (!redo)  {
        print( "Warning: aggregated_data is loading from a saved instance ... add redo=TRUE if data needs to be refresh" )
        if (file.exists(fn)) {
          load( fn)
          return( sppoly )
        }
        print( "Warning: aggregated_data load from saved instance failed ... " )
      }

      print( "Warning: aggregated_data is being recreated ... " )
      print( "Warning: this needs a lot of RAM .. ~XX GB depending upon resolution of discretization .. a few hours " )
      B = bathymetry.db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
      crs_lonlat = sp::CRS("+proj=longlat +datum=WGS84")
      spp = sppoly["StrataID"]
      B$StrataID = over( SpatialPoints( B[, c("lon", "lat")], crs_lonlat ), spTransform(spp, crs_lonlat ) )$StrataID # match each datum to an area
      o = NULL
      B$lon = NULL
      B$lat = NULL
      gc()
      bb = as.data.frame( t( simplify2array(
        tapply( X=B$z, INDEX=list(paste( B$StrataID) ),
          FUN = function(w) { c(
            mean(w, na.rm=TRUE),
            sd(w, na.rm=TRUE),
            length( which(is.finite(w)) )
          ) }, simplify=TRUE )
      )))
      B = NULL
      colnames(bb) = c("z.mean", "z.sd", "z.n")
      bb$StrataID = rownames(bb)
      sppoly$z = NA
      sppoly$z.sd = NA
      sppoly$z.n = NA
      j = match( bb$StrataID, sppoly$StrataID )
      if (length(j) > 0)  {
        sppoly$z[j] = bb$z.mean
        sppoly$z.sd[j] = bb$z.sd
        sppoly$z.n[j] = bb$z.n
      }
      save( sppoly, file=fn, compress=TRUE )
      return( sppoly )
    }




    if ( DS %in% c("carstm_modelled", "carstm_modelled_fit") ) {

      fn = file.path( p$modeldir, paste( "bathymetry", "carstm_modelled", id, "rdata", sep=".") )
      fn_fit = file.path( p$modeldir, paste( "bathymetry", "carstm_modelled_fit", id, "rdata", sep=".") )

      if (!redo)  {
        print( "Warning: carstm_modelled is loading from a saved instance ... add redo=TRUE if data needs to be refresh" )
        if (DS=="carstm_modelled") {
          if (file.exists(fn)) {
            load( fn)
            return( sppoly )
          }
        }
        if (DS=="carstm_modelled_fit") {
          if (file.exists(fn_fit)) {
            load( fn_fit )
            return( fit )
          }
        }
        print( "Warning: carstm_modelled load from saved instance failed ... " )
      }

      print( "Warning: carstm_modelled is being recreated ... " )
      print( "Warning: this needs a lot of RAM .. ~XX GB depending upon resolution of discretization .. a few hours " )
      M = bathymetry.db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
      crs_lonlat = sp::CRS("+proj=longlat +datum=WGS84")
      spp = sppoly["StrataID"]

      # do this immediately to reduce storage for spp (before adding other variables)
      M$StrataID = as.character( over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(spp, crs_lonlat ) )$StrataID) # match each datum to an area
      o = NULL
      M$lon = NULL
      M$lat = NULL
      M$z = M$z + 2500 # make all positive
      M$tag = "observations"

      spp = as.data.frame(spp)
      spp$z = NA
      spp$StrataID = as.character( spp$StrataID )
      spp$tag ="predictions"

      M = rbind( M, spp[, names(M)] )
      M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
      M$strata  = as.numeric( M$StrataID)
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      spp = NULL
      gc()

      H = carstm_hyperparameters( sd(log(M$z), na.rm=TRUE), alpha=0.5, median( log(M$z), na.rm=TRUE) )
      fit = inla(
        formula =
          z ~ 1 +
            + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
            + f(iid_error, model="iid", hyper=H$iid)
          ,
        family = "lognormal",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
        # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        num.threads=2,
        blas.num.threads=2,
        verbose=TRUE
      )
      save( fit, file=fn_fit, compress=TRUE )
      s = summary(fit)
      s$dic$dic  # 31225
      s$dic$p.eff # 5200

      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

      # reformat predictions into matrix form
      ii = which(M$tag=="predictions")
      jj = match(M$StrataID[ii], sppoly$StrataID)
      sppoly@data$z.predicted = exp( fit$summary.fitted.values[ ii[jj], "mean" ]) - 2500
      sppoly@data$z.predicted_lb = exp( fit$summary.fitted.values[ ii[jj], "0.025quant" ]) - 2500
      sppoly@data$z.predicted_ub = exp( fit$summary.fitted.values[ ii[jj], "0.975quant" ]) - 2500
      sppoly@data$z.random_strata_nonspatial = exp( fit$summary.random$strata[ jj, "mean" ])
      sppoly@data$z.random_strata_spatial = exp( fit$summary.random$strata[ jj+max(jj), "mean" ])
      sppoly@data$z.random_sample_iid = exp( fit$summary.random$iid_error[ ii[jj], "mean" ])
      save( spplot, file=fn, compress=TRUE )

      if (map) {
        vn = "z.predicted"
        brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
        dev.new()
        spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )
      }

      return( sppoly )

    }

  }
