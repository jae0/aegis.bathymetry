
isobath.db = function( ip=NULL, p=NULL, depths=c(100, 200), DS="isobath", project_to=projection_proj4string("lonlat_wgs84"), data_dir=project.datadirectory( "aegis", "bathymetry" ) ) {
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
        if ( proj4string( isobaths ) != as.character(project_to) ) isobaths = spTransform( isobaths, sp::CRS( project_to ) )
        return( isobaths[ as.character(depths) ] )
      }
    }

    p0 = spatial_parameters( spatial_domain=p$spatial_domain )
    depths = sort( unique(c(depths, notfound) ))
    x=seq(min(p0$corners$plon), max(p0$corners$plon), by=p0$pres)
    y=seq(min(p0$corners$plat), max(p0$corners$plat), by=p0$pres)

    Z = bathymetry.db( p=p0, DS="complete", varnames=c("plon", "plat", "z") )
    Zi = array_map( "xy->2", Z[, c("plon", "plat")], gridparams=p0$gridparams )
    Zm = matrix( NA, ncol=p0$nplats, nrow=p0$nplons )
    Zm[Zi] = Z$z
    rm(Z); gc()

    # Zm = fields::image.smooth( Zm, theta=p0$pres, dx=p0$pres, dy=p0$pres ) # a little smoothed to make contours cleaner     .. too slow

    cl = contourLines( x=x, y=y, Zm, levels=depths )

    isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS( p$aegis_proj4string_planar_km ) )
    row.names(slot(isobaths, "data")) = as.character(depths)
    for (i in 1:length(depths)) slot( slot(isobaths, "lines")[[i]], "ID") = as.character(depths[i])
    isobaths = as.SpatialLines.SLDF( isobaths )
    sp::proj4string( isobaths ) =  p$aegis_proj4string_planar_km   # project_to gets reset .. not sure why
    isobaths = spTransform( isobaths, sp::CRS(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format

    save( isobaths, file=fn.iso, compress=TRUE) # save spherical
    if ( ! proj4string( isobaths ) == as.character( project_to) ) isobaths = spTransform( isobaths, sp::CRS( project_to   ) )
    return( isobaths )
  }

  # ------------------------

  if (DS %in% c( "coastLine", "coastLine.redo")) {
    #\\ synomym for coastline_db ... left for historical compatibility .. deprecated
    if (DS=="coastline") return( coastline_db( p=p, DS="mapdata.coastLine", project_to = project_to   ) )
    if (DS=="coastline.redo") return( coastline_db( p=p, DS="mapdata.coastLine.redo", project_to = project_to   ) )
  }

  # ------------------------

  if (DS %in% c("coastPolygon", "coastPolygon.redo") ) {
    #\\ synomym for coastline_db ... left for historical compatibility .. deprecated
    if (DS=="coastPolygon") return( coastline_db( p=p, DS="mapdata.coastPolygon", project_to = project_to   ) )
    if (DS=="coastPolygon.redo") return( coastline_db( p=p, DS="mapdata.coastPolygon.redo", project_to = project_to   ) )
  }


}
