
isobath_db = function( 
  depths=c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 ),
  DS="isobath",
  project_to=projection_proj4string("lonlat_wgs84"),
  data_dir=project.datadirectory( "aegis", "bathymetry" ) ) {

  #\\ create or return isobaths and coastlines/coast polygons
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

    Zm = bathymetry_db( p=p0, DS="aggregated_data_as_matrix" )
    # Zm = fields::image.smooth( Zm, theta=p0$pres, dx=p0$pres, dy=p0$pres ) # a little smoothed to make contours cleaner     .. too slow
    cl = contourLines( x=x, y=y, Zm, levels=depths )

    isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS( p0$aegis_proj4string_planar_km ) )
    isobaths = as( isobaths, "sf")
    st_crs(isobaths) = st_crs( p0$aegis_proj4string_planar_km  ) 

    isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
    row.names(isobaths) = as.character(depths)

    save( isobaths, file=fn.iso, compress=TRUE)

    if ( ! st_crs( isobaths ) == st_crs( project_to) ) isobaths = st_transform( isobaths, st_crs( project_to ) )

    return( isobaths )
  }

  # ------------------------

  if (DS %in% c( "coastLine", "coastLine.redo")) {
    #\\ synomym for coastline_db ... left for historical compatibility .. deprecated
    if (DS=="coastline") return( coastline_db( project_to = project_to   ) )
    # if (DS=="coastline.redo") return( coastline_db( p=p, DS="mapdata.coastLine.redo", project_to = project_to   ) )
  }

  # ------------------------

  if (DS %in% c("coastPolygon", "coastPolygon.redo") ) {
    #\\ synomym for coastline_db ... left for historical compatibility .. deprecated
    if (DS=="coastPolygon") return( coastline_db( project_to = project_to   ) )
    # if (DS=="coastPolygon.redo") return( coastline_db( p=p, DS="mapdata.coastPolygon.redo", project_to = project_to   ) )
  }


}
