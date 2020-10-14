

### stmv/carstm hybrid:
### generate a modelled surface using areal units placed on a lattice system (1km x 1 km grid)
### CARSTM-based, however, run within a range of influence of dimension defined by : p$stmv_interpolation_basis_distance_choices = 5
### i.e., the same size as the stats grid (this is 1/2 of the window so 5 km surround each stats node )

p = aegis.bathymetry::bathymetry_parameters( project_class="model" )

if (0)  bathymetry_db( p=p, DS="stmv_inputs_highres_redo" )  # recreate fields for .. requires 60GB+

stmv( p=p )  # This will take from 40-70 hrs, depending upon system

# quick view
  predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
  statistics  = stmv_db( p=p, DS="stmv.stats" )
  locations   = spatial_grid( p )

  # comparisons
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

  # stats
  # p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
  dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )

  dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
  dev.new(); levelplot( statistics[,match("sdSpatial", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
  dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange

# water only
o = which( predictions>0 & predictions <1000)
#levelplot( log( statistics[o,match("sdTotal", p$statsvars)] ) ~ locations[o,1] + locations[o,2], aspect="iso" ) #sd total
levelplot( log(predictions[o]) ~ locations[o,1] + locations[o,2], aspect="iso" )


# bring together stats and predictions and any other required computations: slope and curvature
# and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# .. still uses about 30-40 GB as the base layer is "superhighres" ..
# if parallelizing .. use different servers than local nodes


bathymetry_db( p=p, DS="complete.redo", spatial_domain_subareas = c( "canada.east.highres", "canada.east",  "SSE", "SSE.mpa" , "snowcrab") ) # finalise at diff resolutions 15 min ..


# ----
### demonstration of the map area unit problem example:

maup = map_area_unit_problem( just_return_results=TRUE )  #default is bathymetry data

  if (0) {

    x = maup$resolution
    # x = log(  maup$n / maup$resolution^2  )  # data density
    yrange = range( maup$min, maup$max )
    plot(mean ~ x, maup, pch=20, ylim=yrange)
    lines( median ~ x, maup, col="green", lwd=1)
    lines(min ~ x, maup, col="red", lwd=4)
    lines(max ~ x, maup, col="blue", lwd=4)
    lines( I(mean+sd) ~ x, maup, col="gray", lwd=2)
    lines( I(mean-sd) ~ x, maup, col="gray", lwd=2)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

    plot( sd~ x, maup)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

  }




# end
