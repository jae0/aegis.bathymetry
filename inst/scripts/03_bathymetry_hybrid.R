

  ### stmv/carstm hybrid:  slow
  
  ### generate a modelled surface using areal units placed on a lattice system (1km x 1 km grid)
  ### CARSTM-based, however, run within a range of influence of dimension defined by : p$stmv_distance_interpolation = 5
  ### i.e., the same size as the stats grid (this is 1/2 of the window so 5 km surround each stats node )

  p = aegis.bathymetry::bathymetry_parameters( project_class="hybrid" )  # "hybrid" uses the "best" interpolation method given data constraints

    # stmv_au_distance_reference="completely_inside_boundary",
    # stmv_au_buffer_links=1,
    # stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
    # stmv_distance_statsgrid = 1, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    # stmv_distance_prediction_limits =c( 5, 10 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
    # stmv_nmin = 50,  # min number of data points req before attempting to model in a localized space
    # stmv_nmax = 600, # no real upper bound.. just speed /RAM


  ncores = ram_local( "ncores", ram_main=?, ram_process=? ) # in GB; 

  p$stmv_runmode$carstm = rep("localhost", ncores)


  if (redo_inouts) {
      bathymetry_db( p=p, DS="stmv_inputs_redo" )  # recreate fields for .. requires 60GB+
      bathymetry_db( p=p, DS="stmv_inputs_highres_redo" )  # recreate fields for .. requires 60GB+
  }

  stmv( p=p )  # This will take from 40-70 hrs, depending upon system
  # stmv_db( p=p, DS="cleanup.all" )
      
  # quick view
  predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
  statistics  = stmv_db( p=p, DS="stmv.stats" )
  locations   = spatial_grid( p )

  predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
  statistics  = stmv_db( p=p, DS="stmv.stats" )
  locations =  spatial_grid( p )

  # comparison
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )





  # bring together stats and predictions and any other required computations: slope and curvature
  # and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
  # .. still uses about 30-40 GB as the base layer is "superhighres" ..
  # if parallelizing .. use different servers than local nodes

  bathymetry_db( p=p, DS="complete.redo", spatial_domain_subareas = c( "canada.east.highres", "canada.east",  "SSE", "SSE.mpa" , "snowcrab") ) # finalise at diff resolutions 15 min ..


  # some plots

  # comparisons
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

  # stats
  dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )

  # statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
  statsvars = dimnames(statistics)[[2]]


  dev.new(); levelplot( statistics[,match("Phi_for_aui", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
  dev.new(); levelplot( statistics[,match("rsquared", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu

  dev.new(); levelplot( statistics[,match("sdSpatial", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
  dev.new(); levelplot( statistics[,match("sdTotal", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total

  # water only
  o = which( predictions>5 & predictions <1000)
  #levelplot( log( statistics[o,match("sdTotal", statsvars)] ) ~ locations[o,1] + locations[o,2], aspect="iso" ) #sd total
  levelplot( (predictions[o]) ~ locations[o,1] + locations[o,2], aspect="iso" )




  # end
