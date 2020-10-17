

### stmv/carstm hybrid:  seems to be the "best" interpolation method .. in this case, due to data density

### generate a modelled surface using areal units placed on a lattice system (1km x 1 km grid)
### CARSTM-based, however, run within a range of influence of dimension defined by : p$stmv_interpolation_basis_distance_choices = 5
### i.e., the same size as the stats grid (this is 1/2 of the window so 5 km surround each stats node )

p = aegis.bathymetry::bathymetry_parameters( project_class="production" )  # "production" uses the "best" interpolation method given data constraints

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

p$statsvars = dimnames(statistics)[[2]]


dev.new(); levelplot( statistics[,match("Phi_for_aui", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("rsquared", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu

dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdSpatial", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange

# water only
o = which( predictions>5 & predictions <1000)
#levelplot( log( statistics[o,match("sdTotal", p$statsvars)] ) ~ locations[o,1] + locations[o,2], aspect="iso" ) #sd total
levelplot( (predictions[o]) ~ locations[o,1] + locations[o,2], aspect="iso" )


# bring together stats and predictions and any other required computations: slope and curvature
# and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# .. still uses about 30-40 GB as the base layer is "superhighres" ..
# if parallelizing .. use different servers than local nodes


bathymetry_db( p=p, DS="complete.redo", spatial_domain_subareas = c( "canada.east.highres", "canada.east",  "SSE", "SSE.mpa" , "snowcrab") ) # finalise at diff resolutions 15 min ..


# end
