
bathymetry_carstm = function(p=NULL, DS=NULL, sppoly=NULL, redo=FALSE, map_results=FALSE, ...) {

  #\\ Note inverted convention: depths are positive valued
  #\\ i.e., negative valued for above sea level and positive valued for below sea level

  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable



  # ----------------------


  if ( DS=="carstm_inputs") {

    fn = file.path( p$modeldir, paste( "bathymetry", "carstm_inputs", p$auid, "rdata", sep=".") )
    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }

    # prediction surface
    if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found
    sppoly = sppoly["StrataID"]

    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))

    # do this immediately to reduce storage for sppoly (before adding other variables)

    M = bathymetry.db ( p=p, DS="aggregated_data" )  # 16 GB in RAM just to store!

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area

    names(M)[which(names(M)=="z.mean" )] = "z"

    M$tag = "observations"

    sppoly_df = as.data.frame(sppoly)
    sppoly_df$z = NA
    sppoly_df$StrataID = as.character( sppoly_df$StrataID )
    sppoly_df$tag ="predictions"

    vn = c("z", "tag", "StrataID")

    M = rbind( M[, vn], sppoly_df[, vn] )
    sppoly_df = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    sppoly = NULL
    M$strata  = as.numeric( M$StrataID)
    M$iid_error = 1:nrow(M) # for inla indexing for set level variation

    save( M, file=fn, compress=TRUE )
    return( M )
  }


  # ------------


  if ( DS %in% c("carstm_modelled", "carstm_modelled_fit") ) {

    fn = file.path( p$modeldir, paste( "bathymetry", "carstm_modelled", p$auid, p$carstm_modelengine, "rdata", sep=".") )
    fn_fit = file.path( p$modeldir, paste( "bathymetry", "carstm_modelled_fit", p$auid, p$carstm_modelengine, "rdata", sep=".") )

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

    # prediction surface
    if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found

    M = bathymetry_carstm( p=p, DS="carstm_inputs" )  # will redo if not found
    M$z = M$z + p$constant_offset # make all positive

    fit  = NULL

    if ( grepl("glm", p$carstm_modelengine) ) {

      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) error("model fit error")
      if ("try-error" %in% class(fit) ) error("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      # s = summary(fit)
      # AIC(fit)  # 104487274
      # reformat predictions into matrix form
      ii = which( M$tag=="predictions" & M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      sppoly@data[,"z.predicted"] = exp( preds$fit) - p$constant_offset
      sppoly@data[,"z.predicted_se"] = exp( preds$se.fit)
      sppoly@data[,"z.predicted_lb"] = exp( preds$fit - preds$se.fit ) - p$constant_offset
      sppoly@data[,"z.predicted_ub"] = exp( preds$fit + preds$se.fit ) - p$constant_offset
      save( sppoly, file=fn, compress=TRUE )
    }

    if ( grepl("gam", p$carstm_modelengine) ) {
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) error("model fit error")
      if ("try-error" %in% class(fit) ) error("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      s = summary(fit)
      AIC(fit)  # 104487274
      # reformat predictions into matrix form
      ii = which( M$tag=="predictions" & M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      sppoly@data[,"z.predicted"] = exp( preds$fit) - p$constant_offset
      sppoly@data[,"z.predicted_se"] = exp( preds$se.fit)
      sppoly@data[,"z.predicted_lb"] = exp( preds$fit - preds$se.fit ) - p$constant_offset
      sppoly@data[,"z.predicted_ub"] = exp( preds$fit + preds$se.fit ) - p$constant_offset
      save( sppoly, file=fn, compress=TRUE )
    }


    if ( grepl("inla", p$carstm_modelengine) ) {

      H = carstm_hyperparameters( sd(log(M$z), na.rm=TRUE), alpha=0.5, median( log(M$z), na.rm=TRUE) )
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) error("model fit error")
      if ("try-error" %in% class(fit) ) error("model fit error")
      save( fit, file=fn_fit, compress=TRUE )


      # plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

      # reformat predictions into matrix form
      ii = which(M$tag=="predictions")
      jj = match(M$StrataID[ii], sppoly$StrataID)
      sppoly@data$z.predicted = exp( fit$summary.fitted.values[ ii[jj], "mean" ]) - p$constant_offset
      sppoly@data$z.predicted_lb = exp( fit$summary.fitted.values[ ii[jj], "0.025quant" ]) - p$constant_offset
      sppoly@data$z.predicted_ub = exp( fit$summary.fitted.values[ ii[jj], "0.975quant" ]) - p$constant_offset
      sppoly@data$z.random_strata_nonspatial = exp( fit$summary.random$strata[ jj, "mean" ])
      sppoly@data$z.random_strata_spatial = exp( fit$summary.random$strata[ jj+max(jj), "mean" ])
      sppoly@data$z.random_sample_iid = exp( fit$summary.random$iid_error[ ii[jj], "mean" ])
      save( sppoly, file=fn, compress=TRUE )
    }


    if (map_results) {
      vn = "z.predicted"
      brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
      dev.new();  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )
    }

    return( sppoly )

  }

}