
isobath_db = function( 
  p = NULL,
  depths=c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 ),
  DS="isobath",
  project_to=projection_proj4string("lonlat_wgs84"),
  data_dir=project.datadirectory( "aegis", "bathymetry" ),
  add_missing=FALSE,
  p_max_resolution="canada.east.superhighres",
  aRange = 3  # # pixels to approx 1 SD ,
   ) {

  #\\ create or return isobaths and coastlines/coast polygons
  #\\ isobaths come from aggregated data (resolution of pres) which is then locally smoothed through a guassian kernal process 
  # require(stmv)
  if (DS %in% c( "isobath", "isobath.redo" )) {

    require (fields)
    
    if (is.null(p)) p = aegis.bathymetry::bathymetry_parameters() 
    pDomain = aegis.bathymetry::bathymetry_parameters( spatial_domain=p_max_resolution )

    fn.iso = file.path( data_dir, "isobaths", paste("isobaths", p$spatial_domain, "rdata", sep=".") )  # in case there is an alternate project

    isobaths = NULL
    notfound = NULL

    if ( DS == "isobath" ) {
      if (file.exists(fn.iso)) {
        load(fn.iso)
        isobaths = as( isobaths, "sf")
        st_crs(isobaths) = st_crs( p$aegis_proj4string_planar_km  ) 

        nn = row.names(isobaths)
        notfound = setdiff( as.character(depths), nn )
        if (length( notfound) > 0 ) {
          message( "Some isobaths not found ... add them:")
          print(notfound)
        }
        if (!add_missing) {
          if ( st_crs( isobaths ) != st_crs(project_to) ) isobaths = st_transform( isobaths, st_crs( project_to ) )
          return( isobaths[ which( nn %in% as.character(depths)), ] )
        } else {
          message( "adding missing depths to current data object .. " )
        }
      }
    }


    pres = pDomain$pres

    x = seq(min(pDomain$corners$plon), max(pDomain$corners$plon), by=pres)
    y = seq(min(pDomain$corners$plat), max(pDomain$corners$plat), by=pres)

    options( max.contour.segments=50000 )

    if ( !is.null(notfound) ) {

      depths = sort( unique( as.numeric( notfound ) ) )
      message("Adding isobaths ..")
      Zsmoothed = attributes( isobaths)$Zsmoothed 

    } else {
      
      depths = sort( unique( depths ) )
      Z = bathymetry_db( p=pDomain, DS="aggregated_data" )
      Zi = array_map( "xy->2", Z[, c("plon", "plat")], gridparams=pDomain$gridparams )

      # remove raw data outside of the bounding box
        good = which( Zi[,1] >= 1 & Zi[,1] <= pDomain$nplons & Zi[,2] >= 1 & Zi[,2] <= pDomain$nplats )
        Zi = Zi[good,]
        Z = Z[good,]

      Zmatrix = matrix(NA, nrow=pDomain$nplons, ncol=pDomain$nplats )
      Zmatrix[Zi] = Z$z.mean
      Zsmoothed = image.smooth( Zmatrix, aRange=aRange )

    }

    cl = contourLines( x=x, y=y, Zsmoothed$z, levels=depths )

    isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS( p$aegis_proj4string_planar_km ) )
    isobaths = as( isobaths, "sf")
    st_crs(isobaths) = st_crs( p$aegis_proj4string_planar_km  ) 

    isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
    row.names(isobaths) = as.character(isobaths$level)

    attr( isobaths, "Zsmoothed" ) = Zsmoothed
    attr( isobaths, "aRange" ) =  aRange

    attr( isobaths, "pres" ) =  pres
    attr( isobaths, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
    attr( isobaths, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

    if ( !is.null(notfound) ) {
      toadd = isobaths
      load(fn.iso)
      isobaths = rbind(isobaths, toadd)           
    }

    save( isobaths, file=fn.iso, compress=TRUE)

    if ( ! st_crs( isobaths ) == st_crs( project_to) ) isobaths = st_transform( isobaths, st_crs( project_to ) )

    return( isobaths )
  }

  # ------------------------

  if (DS %in% c( "coastLine", "coastLine.redo")) {
    #\\ synomym for coastline_db ... left for historical compatibility .. deprecated
    if (DS=="coastline") return( coastline_db( project_to = project_to   ) )

  }

  # ------------------------

  if (DS %in% c("coastPolygon", "coastPolygon.redo") ) {
    #\\ synomym for coastline_db ... left for historical compatibility .. deprecated
    if (DS=="coastPolygon") return( coastline_db( project_to = project_to   ) )
  }


}
