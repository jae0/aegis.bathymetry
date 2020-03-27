
isobath.db = function( ip=NULL, p=NULL, depths=c(100, 200), DS="isobath", crs=projection_proj4string("lonlat_wgs84", data_dir=project.datadirectory( "aegis", "bathymetry" ) ) ) {
  #\\ create or return isobaths and coastlines/coast polygons
  # require(stmv)
  if (DS %in% c( "isobath", "isobath.redo" )) {

    fn.iso = file.path( data_dir, "isobaths", paste("isobaths", p$spatial_domain, "rdata", sep=".") )  # in case there is an alternate project

    isobaths = NULL
    notfound = NULL

    if ( DS != "isobath.redo" & file.exists(fn.iso) ) {
      load(fn.iso)
      notfound = setdiff( as.character(depths), names(isobaths) )
      if (length( notfound)==0) {
        if ( proj4string( isobaths ) != as.character(crs) ) isobaths = spTransform( isobaths, sp::CRS( crs ) )
        return( isobaths[ as.character(depths) ] )
      }
    }

    p = spatial_parameters( p )
    depths = sort( unique(c(depths, notfound) ))
    x=seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    y=seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)

    Z = bathymetry.db( p=p, DS="complete", varnames=c("plon", "plat", "z") )
    Zi = array_map( "xy->2", Z[, c("plon", "plat")], gridparams=p$gridparams )
    Zm = matrix( NA, ncol=p$nplats, nrow=p$nplons )
    Zm[Zi] = Z$z
    rm(Z); gc()

    # Zm = fields::image.smooth( Zm, theta=p$pres, dx=p$pres, dy=p$pres ) # a little smoothed to make contours cleaner     .. too slow

    cl = contourLines( x=x, y=y, Zm, levels=depths )

    isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS( p$aegis_proj4string_planar_km ) )
    row.names(slot(isobaths, "data")) = as.character(depths)
    for (i in 1:length(depths)) slot( slot(isobaths, "lines")[[i]], "ID") = as.character(depths[i])
    isobaths = as.SpatialLines.SLDF( isobaths )
    sp::proj4string( isobaths ) =  p$aegis_proj4string_planar_km   # crs gets reset .. not sure why
    isobaths = spTransform( isobaths, sp::CRS(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format

    save( isobaths, file=fn.iso, compress=TRUE) # save spherical
    if ( ! proj4string( isobaths ) == as.character( crs) ) isobaths = spTransform( isobaths, sp::CRS( crs ) )
    return( isobaths )
  }

  # ------------------------

  if (DS %in% c( "coastLine", "coastLine.redo")) {
    #\\ synomym for coastline.db ... left for historical compatibility .. deprecated
    if (DS=="coastline") return( coastline.db( p=p, DS="mapdata.coastLine", crs=crs ) )
    if (DS=="coastline.redo") return( coastline.db( p=p, DS="mapdata.coastLine.redo", crs=crs ) )
  }

  # ------------------------

  if (DS %in% c("coastPolygon", "coastPolygon.redo") ) {
    #\\ synomym for coastline.db ... left for historical compatibility .. deprecated
    if (DS=="coastPolygon") return( coastline.db( p=p, DS="mapdata.coastPolygon", crs=crs ) )
    if (DS=="coastPolygon.redo") return( coastline.db( p=p, DS="mapdata.coastPolygon.redo", crs=crs ) )
  }


}
