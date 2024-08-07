
# Bathymetry base data
p = aegis.bathymetry::bathymetry_parameters()  # default params

bathymetry_db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015

for ( dom in p$spatial_domain_subareas ) {
  Z = bathymetry_db( 
    p = aegis.bathymetry::bathymetry_parameters( spatial_domain=dom ), 
    DS ="aggregated_data", 
    redo=TRUE 
  )  
  str(Z)
  Z = NULL; gc()
}


### -----------------------------------------------------------------
# to update/recreate new polygons, run the following:  aegis.bathymetry/inst/scripts/01.bathymetry_data.R

bathyclines.redo = FALSE

if( bathyclines.redo ) {
  # note these polygons are created at the resolution specified in spatial_domain ..
  # which by default is very high ("canada.east.highres" = 0.5 km .. p$pres ).
  # For lower, specify an appropriate p$spatial_domain, specify the p=p, where the latter is specific to the spatial_domain
  # options(max.contour.segments=1000) # might be required if superhighres is being used

  depths = c( 0, 10, 20, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
              1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 )
  
  p = aegis.bathymetry::bathymetry_parameters()  # default params
  doms = c( p$spatial_domain, p$spatial_domain_subareas )
  for ( dom in doms ) {
    plygn = isobath_db( 
      spatial_domain=dom,
      DS="isobath.redo", 
      depths=depths, 
      project_to=projection_proj4string("lonlat_wgs84")  
    )
  }

  
  (plygn)

  plygn = NULL; gc()
  
  depths = c(25,150)
  plygn = isobath_db( DS="isobath", depths=depths, project_to=projection_proj4string("lonlat_wgs84" ),add_missing=TRUE  )

}


if (0) {
    ### -----------------------------------------------------------------
    # some test plots
    RLibrary( "aegis.bathymetry" , "aegis.coastline", "aegis.polygons")

    p = bathymetry_parameters( spatial_domain="canada.east" ) # reset to lower resolution
    depths = c( 100, 200, 300, 500, 1000)
    plygn = isobath_db(  depths=depths  )

    coast = coastline_db( xlim=c(-75,-52), ylim=c(41,50), no.clip=TRUE )  # no.clip is an option for maptools::getRgshhsMap
    plot( coast, col="transparent", border="steelblue2" , xlim=c(-68,-52), ylim=c(41,50),  xaxs="i", yaxs="i", axes=TRUE )  # ie. coastline
    plot( plygn[  as.character(c( 100, 200, 300 ))  ,], col="gray90", add=TRUE ) # for multiple polygons
    plot( plygn[  as.character(c( 500, 1000)) , ], col="gray80", add=TRUE ) # for multiple polygons
    # plot( plygn, xlim=c(-68,-52), ylim=c(41,50))  # all isobaths commented as it is slow ..


    # or to get in projected (planar) coords as defined by p$spatial_domain
    plygn = isobath_db( DS="isobath", depths=c(100) , project_to=p$aegis_proj4string_planar_km ) # as SpatialLines
    plot(plygn)

    plygn_aslist = st_coordinates( plygn)
    plot( 0,0, type="n", xlim=c(-200,200), ylim=c(-200,200)  )
    lapply( plygn_aslist[[1]], points, pch="." )

    plygn_as_xypoints = st_coordinates( plygn  )# ... etc...
    plot(plygn_as_xypoints, pch=".",  xaxs="i", yaxs="i", axes=TRUE)
}
 