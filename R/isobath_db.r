
isobath_db = function( 
  depths=c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 ),
  DS="isobath",
  project_to=projection_proj4string("lonlat_wgs84"),
  data_dir=project.datadirectory( "aegis", "bathymetry" ),
  aRange = 3  # # pixels to approx 1 SD 
   ) {

  #\\ create or return isobaths and coastlines/coast polygons
  #\\ isobaths come from aggregated data (resolution of pres) which is then locally smoothed through a guassian kernal process 
  # require(stmv)
  if (DS %in% c( "isobath", "isobath.redo" )) {

    p0 = bathymetry_parameters()
    fn.iso = file.path( data_dir, "isobaths", paste("isobaths", p0$spatial_domain, "rdata", sep=".") )  # in case there is an alternate project

    isobaths = NULL
    notfound = NULL

    if ( DS != "isobath.redo" & file.exists(fn.iso) ) {
      load(fn.iso)
      notfound = setdiff( as.character(depths), names(isobaths) )
      if (length( notfound)==0) {
        if ( st_crs( isobaths ) != st_crs(project_to) ) isobaths = st_transform( isobaths, st_crs( project_to ) )
        return( isobaths[ which(isobaths$level %in% as.character(depths)), ] )
      }
    }

    depths = sort( unique( depths ) )
    x=seq(min(p0$corners$plon), max(p0$corners$plon), by=p0$pres)
    y=seq(min(p0$corners$plat), max(p0$corners$plat), by=p0$pres)

    Z = bathymetry_db( p=p0, DS="aggregated_data" )

    Zi = array_map( "xy->2", Z[, c("plon", "plat")], gridparams=p0$gridparams )

    # remove raw data outside of the bounding box
      good = which( Zi[,1] >= 1 & Zi[,1] <= p0$nplons & Zi[,2] >= 1 & Zi[,2] <= p0$nplats )
      Zi = Zi[good,]
      Z = Z[good,]

    Zmatrix = matrix(NA, nrow=p0$nplons, ncol=p0$nplats )
    Zmatrix[Zi] = Z$z.mean

    Zsmoothed = image.smooth( Zmatrix, aRange=aRange )
  
    cl = contourLines( x=x, y=y, Zsmoothed$z, levels=depths )

    isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS( p0$aegis_proj4string_planar_km ) )
    isobaths = as( isobaths, "sf")
    st_crs(isobaths) = st_crs( p0$aegis_proj4string_planar_km  ) 

    isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
    row.names(isobaths) = as.character(isobaths$level)

    attr( isobaths, "Zmatrix" ) = Zmatrix
    attr( isobaths, "Zsmoothed" ) = Zsmoothed
    attr( isobaths, "aRange" ) =  aRange

    attr( isobaths, "pres" ) =  p0$pres
    attr( isobaths, "proj4string_planar" ) =  p0$aegis_proj4string_planar_km
    attr( isobaths, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")
         
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
