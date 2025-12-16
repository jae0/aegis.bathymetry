
isobath_db = function( 
  p = NULL,
  spatial_domain="canada.east.superhighres",
  depths=c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 ),
  DS="isobath",
  project_to=projection_proj4string("lonlat_wgs84"),
  data_dir=project.datadirectory( "aegis", "bathymetry" )
   ) {

  #\\ create or return isobaths and coastlines/coast polygons
  #\\ isobaths come from aggregated data (resolution of pres) which is then locally smoothed through a guassian kernal process 
    
  if (DS %in% c( "isobath", "isobath.redo" )) {
    
    # behviour determmined by p. If passed then a p-specific Zsmooth is created otherwise use the highest resoltuion for the region
    require (fields)
    
    isobaths = NULL

    options( max.contour.segments=50000 )

    depths = sort( unique( depths ) )

    if ( is.null(p)) {
    
      if ( DS == "isobath" ) {
    
        fn.iso = file.path( data_dir, "isobaths", paste("isobaths", spatial_domain, "rdz", sep=".") )  # in case there is an alternate project
        
        if (file.exists(fn.iso)) {
          isobaths = read_write_fast(fn.iso)
          # isobaths = as( isobaths, "sf")  # in case an old file from sp*

          nn = row.names(isobaths)
          if ( st_crs( isobaths ) != st_crs(project_to) ) isobaths = st_transform( isobaths, st_crs( project_to ) )

          notfound = setdiff( as.character(depths), nn )

          if (length( notfound) > 0 ) {
            message( "matching isobaths not found, computing on the fly  .. " )
            
            Z =  st_as_sf(Z, coords=c("plon", "plat"), crs=st_crs(p$aegis_proj4string_planar_km) ) 

            Z = stars::st_rasterize( Z["z.mean"], dx=p$pres, dy=p$pres )

            isobaths = st_contour( Z, contour_lines = FALSE, na.rm = TRUE, breaks = depths )
            isobaths = as( isobaths, "sf") 
            isobaths = isobaths[, "Min"]
            names(isobaths) = c("levels" , "geometry")
            st_crs(isobaths) = st_crs( p$aegis_proj4string_planar_km  ) 

            isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
            row.names(isobaths) = as.character(isobaths$level)
  
          }
          
          retain = which( row.names(isobaths) %in% as.character(depths) )
          return( isobaths[retain,]  )
        }
      }
    
    } 
    
    # here is redoing or p is passed and a lower (alt) resolution, p-specific isobath is desired 
    if (is.null(p)) p = aegis.bathymetry::bathymetry_parameters( spatial_domain=spatial_domain ) 

    fn.iso = file.path( data_dir, "isobaths", paste("isobaths", spatial_domain, "rdz", sep=".") )  # in case there is an alternate project

    Z = bathymetry_db( p=p, DS="aggregated_data" )
    
    Z =  st_as_sf(Z, coords=c("plon", "plat"), crs=st_crs(p$aegis_proj4string_planar_km) ) 

    Z = stars::st_rasterize( Z["z.mean"], dx=p$pres, dy=p$pres )

    isobaths = st_contour( Z, contour_lines = FALSE, na.rm = TRUE, breaks = depths )
    isobaths = as( isobaths, "sf") 
    isobaths = isobaths[, "Min"]
    names(isobaths) = c("levels" , "geometry")
    st_crs(isobaths) = st_crs( p$aegis_proj4string_planar_km  ) 

    isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
    row.names(isobaths) = as.character(isobaths$level)

    attr( isobaths, "corners" ) =  p$corners
    attr( isobaths, "pres" ) =  p$pres
    attr( isobaths, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
    attr( isobaths, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

    read_write_fast( isobaths, file=fn.iso)

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
