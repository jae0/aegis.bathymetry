
isobath_db = function( 
  p = NULL,
  spatial_domain="canada.east.superhighres",
  depths=c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 ),
  DS="isobath",
  project_to=projection_proj4string("lonlat_wgs84"),
  data_dir=project.datadirectory( "aegis", "bathymetry" ),
  aRange = 3  # # pixels to approx 1 SD ,
   ) {

  #\\ create or return isobaths and coastlines/coast polygons
  #\\ isobaths come from aggregated data (resolution of pres) which is then locally smoothed through a guassian kernal process 
    
  if (DS %in% c( "isobath", "isobath.redo" )) {
    
    # behviour determmined by p. If passed then a p-specific Zsmooth is created otherwise use the highest resoltuion for the region
    require (fields)
    
    if ( spatial_domain %in% c( "canada.east.superhighres", "canada.east.highres", "canada.east", "SSE", "SSE.mpa" , "snowcrab" ) ) spatial_domain_input = "canada.east.superhighres"
        

    isobaths = NULL

    options( max.contour.segments=50000 )

    depths = sort( unique( depths ) )

    if ( is.null(p)) {
    
      if ( DS == "isobath" ) {
    
        fn.iso = file.path( data_dir, "isobaths", paste("isobaths", spatial_domain_input, "rdata", sep=".") )  # in case there is an alternate project
        
        if (file.exists(fn.iso)) {
          load(fn.iso)
          isobaths = as( isobaths, "sf")  # in case an old file from sp*

          nn = row.names(isobaths)
          if ( st_crs( isobaths ) != st_crs(project_to) ) isobaths = st_transform( isobaths, st_crs( project_to ) )

          notfound = setdiff( as.character(depths), nn )

          if (length( notfound) > 0 ) {
            message( "matching isobaths not found, computing on the fly  .. " )

            Zsmoothed = attributes( isobaths)$Zsmoothed 

            x = seq(min(attributes( isobaths)$corners$plon), max(attributes( isobaths)$corners$plon), by=attributes( isobaths)$pres)
            y = seq(min(attributes( isobaths)$corners$plat), max(attributes( isobaths)$corners$plat), by=attributes( isobaths)$pres)
            
            {
            message( "FIX ME:: maptools depreciated use stars::st_contour")
            # cl = stars::st_contour( stars, breaks= depths )
            cl = contourLines( x=x, y=y, Zsmoothed$z, levels=depths )

            iso_crs = attributes( isobaths)$proj4string_planar 
            
            isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS(iso_crs) )
            }

            isobaths = as( isobaths, "sf")
            st_crs(isobaths) = st_crs( iso_crs )  

            isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
            row.names(isobaths) = as.character(isobaths$level)
  
          }
          retain = which( row.names(isobaths) %in% as.character(depths) )
          return( isobaths[retain,]  )
        }
      }
    
    } 
    
    # here is redoing or p is passed and a lower (alt) resolution, p-specific isobath is desired 
    if (is.null(p)) p = aegis.bathymetry::bathymetry_parameters( spatial_domain=spatial_domain_input ) 

    fn.iso = file.path( data_dir, "isobaths", paste("isobaths", spatial_domain, "rdata", sep=".") )  # in case there is an alternate project

    Z = bathymetry_db( p=p, DS="aggregated_data" )
    Zi = array_map( "xy->2", Z[, c("plon", "plat")], gridparams=p$gridparams )

    # remove raw data outside of the bounding box
    good = which( Zi[,1] >= 1 & Zi[,1] <= p$nplons & Zi[,2] >= 1 & Zi[,2] <= p$nplats )
    Zi = Zi[good,]
    Z = Z[good,]

    Zmatrix = matrix(NA, nrow=p$nplons, ncol=p$nplats )
    Zmatrix[Zi] = Z$z.mean
    Zsmoothed = image.smooth( Zmatrix, aRange=aRange )

    x = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    y = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)

{
    message( "FIX ME !!! :: maptools depreciated use stars::st_contour")
    # ?stars::st_as_stars
    # cl = stars::st_contour( stars, breaks= depths )

    cl = contourLines( x=x, y=y, Zsmoothed$z, levels=depths )

    isobaths = maptools::ContourLines2SLDF(cl, proj4string=sp::CRS( p$aegis_proj4string_planar_km ) )
}
    isobaths = as( isobaths, "sf")
    st_crs(isobaths) = st_crs( p$aegis_proj4string_planar_km  ) 

    isobaths = st_transform( isobaths, st_crs(projection_proj4string("lonlat_wgs84")) )  ## longlat  as storage format
    row.names(isobaths) = as.character(isobaths$level)

    attr( isobaths, "Zsmoothed" ) = Zsmoothed
    attr( isobaths, "aRange" ) =  aRange
    attr( isobaths, "corners" ) =  p$corners
    attr( isobaths, "pres" ) =  p$pres
    attr( isobaths, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
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
