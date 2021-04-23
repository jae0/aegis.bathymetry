
# Bathymetry base data
p = aegis.bathymetry::bathymetry_parameters()  # default params

bathymetry_db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015


Z = bathymetry_db( p=p, DS="aggregated_data", redo=TRUE )  
Z = NULL; gc()

Z = bathymetry_db( p=p, DS="aggregated_data_as_matrix", redo=TRUE )  
Z = NULL; gc()


### -----------------------------------------------------------------
# to update/recreate new polygons, run the following:
bathyclines.redo = FALSE

if( bathyclines.redo ) {
  # note these polygons are created at the resolution specified in p$spatial_domain ..
  # which by default is very high ("canada.east.highres" = 0.5 km .. p$pres ).
  # For lower one specify an appropriate p$spatial_domain
  # options(max.contour.segments=1000) # might be required if superhighres is being used

  depths = c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
              1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 )
  plygn = isobath_db( DS="isobath.redo", depths=depths, project_to=projection_proj4string("lonlat_wgs84")  )

}


if (0) {
    ### -----------------------------------------------------------------
    # some test plots
    RLibrary( "aegis.bathymetry" , "aegis.coastline", "aegis.polygons")

    p = bathymetry_parameters( spatial_domain="canada.east" ) # reset to lower resolution
    depths = c( 100, 200, 300, 500, 1000)
    plygn = isobath_db( DS="isobath", depths=depths  )

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


# The rest are for underlying areal units and data for modelling of bathymetry via carstm lattices
require(carstm)
require(aegis)
require(aegis.bathymetry)

for ( areal_units_resolution_km in c(100, 50, 25, 20, 15, 10, 5, 1) ) {
  # for ( spatial_domain in c("snowcrab", "SSE")) {
   for ( spatial_domain in c( "SSE", "snowcrab", "canada.east.superhighres", "canada.east.highres", "canada.east" )) {
    p = bathymetry_parameters(
      spatial_domain = spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_resolution_km = areal_units_resolution_km # km dim of lattice ~ 1 hr
      # areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm  --  uses p$aegis_proj4string_planar_km if not set ..

    )
    sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    M = bathymetry_db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  }
}
