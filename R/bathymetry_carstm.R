
bathymetry_carstm = function(p=NULL, DS=NULL, sppoly=NULL, id=NULL, redo=FALSE, map=FALSE, ...) {

  #\\ Note inverted convention: depths are positive valued
  #\\ i.e., negative valued for above sea level and positive valued for below sea level
  if ( is.null(p)) p = bathymetry_parameters(...)

  if ( !exists("project_name", p)) p$project_name = "bathymetry"
  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )


  if (is.null(id)) id = paste( p$spatial_domain, p$areal_units_overlay, p$areal_units_resolution_km, p$areal_units_strata_type, sep="_" )


  # -----------------

  if (DS=="areal_units") {

    fn = file.path( p$modeldir, paste( "areal_units", id, "rdata", sep=".") )
    sppoly = NULL
    if (!redo) {
      if (file.exists(fn)) load(fn)
      return(sppoly)
    }

    sppoly = areal_units(
      spatial_domain=p$spatial_domain,
      areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km,
      areal_units_strata_type=p$areal_units_strata_type,
      areal_units_resolution_km=p$areal_units_resolution_km,
      areal_units_overlay= ifelse(!exists("areal_units_overlay", p) || !is.finite(p$areal_units_overlay) || !p$areal_units_overlay, "none", p$areal_units_overlay),
      areal_units_constraint=ifelse(!exists("areal_units_", p) || !is.finite(p$areal_units_) || !p$areal_units_, "none", p$areal_units_),
      redo=TRUE
    )

    W.nb = poly2nb(sppoly, row.names=sppoly$StrataID, queen=TRUE)  # slow .. ~1hr?
    W.remove = which(card(W.nb) == 0)

    if ( length(W.remove) > 0 ) {
      # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
      W.keep = which(card(W.nb) > 0)
      W.nb = nb_remove( W.nb, W.remove )
      sppoly = sppoly[W.keep,]
      row.names(sppoly) = as.character(sppoly$StrataID)
      sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
      sppoly$StrataID = factor( as.character(sppoly$StrataID) )
      sppoly$strata = as.numeric( sppoly$StrataID )
      sppoly = sppoly[order(sppoly$strata),]
    }

    attr(sppoly, "nb") = W.nb  # adding neighbourhood as an attribute to sppoly
    save(sppoly, file=fn, compress=TRUE)
    return( sppoly )
  }


  # ------------

  if ( DS=="carstm_inputs") {


    fn = file.path( p$modeldir, paste( "bathymetry", "carstm_inputs", id, "rdata", sep=".") )
    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }
    M = bathymetry.db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]

    if (exists("inputdata_spatial_discretization_planar_km", p)) {
      # thin data a bit ... remove potential duplicates and robustify
      M = lonlat2planar( M, proj.type=p$aegis_proj4string_planar_km )
      M$plon = round(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      M$plat = round(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      keep = which( !duplicated(paste( M$plon, M$plat )) )
      M = M[ keep , ]
      keep = NULL
      gc()
      M$plon = NULL
      M$plat = NULL
      M = M[ which( is.finite( M$z )) ,]
    }

    # prediction surface
    sppoly = bathymetry_carstm( p=p, DS="areal_units" )  # will redo if not found
    sppoly = sppoly["StrataID"]

    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))

    # do this immediately to reduce storage for sppoly (before adding other variables)
    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
    M$lon = NULL
    M$lat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area
    M$z = M$z + p$constant_offset # make all positive
    M$tag = "observations"

    sppoly_df = as.data.frame(sppoly)
    sppoly_df$z = NA
    sppoly_df$StrataID = as.character( sppoly_df$StrataID )
    sppoly_df$tag ="predictions"

    M = rbind( M, sppoly_df[, names(M)] )
    sppoly_df = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    sppoly = NULL
    M$strata  = as.numeric( M$StrataID)
    M$iid_error = 1:nrow(M) # for inla indexing for set level variation

    save( M, file=fn, compress=TRUE )
    return( M )
  }

  # -----------------------

  if ( DS=="aggregated_data") {
    #\\ not used (yet) .. jusr empirical averages over sppoly grids

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

    sppoly = bathymetry_carstm( p=p, DS="areal_units"  )  # will redo if not found
    sppoly = sppoly["StrataID"]

    # B$StrataID = over( SpatialPoints( B[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area

    B$lon = NULL
    B$lat = NULL
    B = B[ which(is.finite(B$StrataID))]

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


  # ------------


  if ( DS %in% c("carstm_modelled", "carstm_modelled_fit") ) {

    fn = file.path( p$modeldir, paste( "bathymetry", "carstm_modelled", id, p$carstm_modelengine, p$carstm_family, "rdata", sep=".") )
    fn_fit = file.path( p$modeldir, paste( "bathymetry", "carstm_modelled_fit", id, p$carstm_modelengine, "rdata", sep=".") )

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

    gc()

    # prediction surface
    sppoly = bathymetry_carstm( p=p, DS="areal_units" )  # will redo if not found
#    sppoly = sppoly["StrataID"]

    M = bathymetry_carstm( p=p, DS="carstm_inputs" )  # will redo if not found

    if (p$carstm_modelengine == "glm" ) {
      fit = glm( formula = z ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ] )
      s = summary(fit)
      AIC(fit)  # 104487274
      # reformat predictions into matrix form
      ii = which( M$tag=="predictions" & M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      # out = reformat_to_matrix(
      #   input = preds$fit,
      #   matchfrom = list( StrataID=M$StrataID[ii] ),
      #   matchto   = list( StrataID=sppoly$StrataID  )
      # )
      # iy = match( as.character(sppoly$StrataID), aps$StrataID )
        sppoly@data[,"z.predicted"] = exp( preds$fit) - p$constant_offset
        sppoly@data[,"z.predicted_se"] = exp( preds$se.fit)
        sppoly@data[,"z.predicted_lb"] = exp( preds$fit - preds$se.fit ) - p$constant_offset
        sppoly@data[,"z.predicted_ub"] = exp( preds$fit + preds$se.fit ) - p$constant_offset
      }

      if ( p$carstm_modelengine == "gam"  ) {
        fit = gam( formula = z ~ 1 + StrataID,  family = gaussian(link="log"), data= M[ which(M$tag=="observations"), ] )
        s = summary(fit)
        AIC(fit)  # 104487274
        # reformat predictions into matrix form
        ii = which( M$tag=="predictions" & M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] )
        preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
        sppoly@data[,"z.predicted"] = exp( preds$fit) - p$constant_offset
        sppoly@data[,"z.predicted_se"] = exp( preds$se.fit)
        sppoly@data[,"z.predicted_lb"] = exp( preds$fit - preds$se.fit ) - p$constant_offset
        sppoly@data[,"z.predicted_ub"] = exp( preds$fit + preds$se.fit ) - p$constant_offset
      }


      # out[ out>1e10] = NA
      # convert numbers/km to biomass/strata (kg)..
      # RES$glm = colSums( {out * sppoly$sa_strata_km2}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
      # RES$glm_cfanorth = colSums( {out * sppoly$cfanorth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
      # RES$glm_cfasouth = colSums( {out * sppoly$cfasouth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
      # RES$glm_cfa4x = colSums( {out * sppoly$cfa4x_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km

      # plot( glm ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")
      # plot( glm_cfanorth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
      # plot( glm_cfasouth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
      # plot( glm_cfa4x ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")

    }

    if (p$carstm_modelengine == "inla") {

      H = carstm_hyperparameters( sd(log(M$z), na.rm=TRUE), alpha=0.5, median( log(M$z), na.rm=TRUE) )

      fit = inla(
        formula = p$carstm_formula,
        family = p$carstm_family,
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
      sppoly@data$z.predicted = exp( fit$summary.fitted.values[ ii[jj], "mean" ]) - p$constant_offset
      sppoly@data$z.predicted_lb = exp( fit$summary.fitted.values[ ii[jj], "0.025quant" ]) - p$constant_offset
      sppoly@data$z.predicted_ub = exp( fit$summary.fitted.values[ ii[jj], "0.975quant" ]) - p$constant_offset
      sppoly@data$z.random_strata_nonspatial = exp( fit$summary.random$strata[ jj, "mean" ])
      sppoly@data$z.random_strata_spatial = exp( fit$summary.random$strata[ jj+max(jj), "mean" ])
      sppoly@data$z.random_sample_iid = exp( fit$summary.random$iid_error[ ii[jj], "mean" ])
      save( spplot, file=fn, compress=TRUE )
    }


    if (map) {
      vn = "z.predicted"
      brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
      dev.new();  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )
    }

    return( sppoly )

  }

}