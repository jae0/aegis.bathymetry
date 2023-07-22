
  landmask = function( lonlat=NULL, db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="index", crs=NULL,
    data_root=file.path( project.datadirectory("aegis"), "bathymetry"), ... ) {
    #\\ Using the world coastline data base:
    #\\ return.value determines what is returned: "coast.lonat", "coast.polygon" (sp),
    #\\ or indices of"not.land", "land"

    message( "defunct ... here for posterity")

    if (is.null(crs)) stop("crs is required")

message("FIXE ME::: deprecated libs, use sf/stars")

    # deprecated: use sf
    # require(maps)
    # require(mapdata)
    # require(maptools)
    # require(rgdal)
    # require(sp)

    if (is.null(crs)) crs = sp::CRS( sp::proj4string(lonlat) )
    if (is.null(crs)) crs = sp::CRS( projection_proj4string("lonlat_wgs84") )

    fno = paste( tag, db, paste0(regions, collapse=""), paste0(crs, collapse=""), "rdata", sep=".")

    defaultdir = project.datadirectory( "aegis", "bathymetry" )

    fn = file.path( data_root, "landmask", fno )
    if (!file.exists( fn  )) fn  = file.path( defaultdir, "landmask", fno )  # in case there is an alternate project

    dir.create( dirname(fn), recursive=TRUE, showWarnings=FALSE  )

    if (is.null( lonlat)) {
      #\\ When lonlat is NULL, this is a flag to return a previous generated and saved version
      #\\ found in bio.data/bathymetry/landmask/ ...
      land = NULL
      if (file.exists( fn)) load(fn)

      if ( return.value=="test" )     return( land)
      if ( return.value=="not.land")  return( which ( is.na(land)) )
      if ( return.value=="land")      return( which ( !is.na(land)) )
    }

{
   message("FIX ME!!! :: map/maptools deprecated " ) 
    # use sf

    coastline = maps::map( database=db, regions=regions, fill=TRUE, plot=FALSE, ...)
    if ( return.value=="coast.lonlat") return (coastline)

    coastlineSp = maptools::map2SpatialPolygons( coastline, IDs=coastline$names, proj4string=crs  )
    if ( return.value=="coast.polygon") return (coastlineSp)

    land = sp::over( SpatialPoints( lonlat, crs ), coastlineSp )

}
    save( land, file=fn, compress=TRUE )
    return(fn)

    if (0) {

        coastline_source="eastcoast_gadm"
        crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
        coast = st_transform( coastline_db( p=p, DS=coastline_source ), crs_lonlat )
        coast$inside = TRUE
        inside = st_points_in_polygons(
          pts = st_as_sf( set[, c("lon", "lat")], coords=c("lon","lat"), crs=crs_lonlat ),
          polys = coast,
          varname = "inside"
        )
    }

  }
