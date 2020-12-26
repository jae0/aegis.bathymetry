
  # Bathymetry

  # Spatial interpolation using stmv on a lattice grid .. continuous form  -- fastest and best performing 

  # Total "superhighres": 2-5 GB/process and 4 GB in parent for fft
  # gam method requires more ~ 2X
  # boundary def takes too long .. too much data to process -- skip
  # "highres": ~ 20 hr with 8, 3.2 Ghz cpus on thoth using fft method jc: 2016 or 2~ 6 hr on hyperion
  # "superhighres" fft: looks to be the best in performance/quality; req ~5 GB per process req
  # FFT is the method of choice for speed and ability to capture the variability
  # krige method is a bit too oversmoothed, especially where rapid changes are occuring

  #  ~24 hrs to scale
  # ~18 + 10 + xxx =  hrs to interpolate
  
  p = aegis.bathymetry::bathymetry_parameters( project_class="stmv" )

  # p$DATA = 'bathymetry_db( p=p, DS="stmv_inputs" )'   # if using lower resolution data
  # p$stmv_distance_statsgrid = 25     # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  
    
  use_parallel_mode = FALSE
  # use_parallel_mode = TRUE
  if (use_parallel_mode) {
      # default is serial mode .. to enable parallel processing, pick and choose:
      scale_ncpus = ram_local( "ncores", ram_main=10, ram_process=4 ) # in GB; about 24  hr
      interpolate_ncpus = ram_local( "ncores", ram_main=2, ram_process=2 ) # nn hrs

      if (!exists("stmv_runmode", p)) p$stmv_runmode = list()
      
      p$stmv_runmode$globalmodel = FALSE
      
      p$stmv_runmode$scale = list(
        cor_0.25 = rep("localhost", scale_ncpus),
        cor_0.1  = rep("localhost", scale_ncpus),
        cor_0.05 = rep("localhost", scale_ncpus),
        cor_0.01 = rep("localhost", scale_ncpus)
      ) 
      
      p$stmv_runmode$interpolate_correlation_basis = list(
        cor_0.25 = rep("localhost", interpolate_ncpus),
        cor_0.1  = rep("localhost", interpolate_ncpus),
        cor_0.05 = rep("localhost", interpolate_ncpus),
        cor_0.01 = rep("localhost", interpolate_ncpus)
      ) 
      
      # p$stmv_runmode$restart_load = "interpolate_correlation_basis"   # only needed if this is restarting from some saved instance

      # if a good idea of autocorrelation is missing, forcing via explicit distance limits is an option
      if (0) {
        p$stmv_runmode$interpolate_distance_basis = list(
          d1 = rep("localhost", interpolate_ncpus),
          d2 = rep("localhost", interpolate_ncpus),
          d3 = rep("localhost", interpolate_ncpus),
          d4 = rep("localhost", interpolate_ncpus),
          d5 = rep("localhost", interpolate_ncpus),
          d6 = rep("localhost", interpolate_ncpus)
        )
      }
      
      p$stmv_runmode$interpolate_predictions = list(
        c1 = rep("localhost", interpolate_ncpus),  
        c2 = rep("localhost", interpolate_ncpus),  
        c3 = rep("localhost", interpolate_ncpus),
        c4 = rep("localhost", interpolate_ncpus),
        c5 = rep("localhost", interpolate_ncpus),
        c6 = rep("localhost", interpolate_ncpus),
        c7 = rep("localhost", interpolate_ncpus)
      )

      p$stmv_runmode$save_intermediate_results = TRUE
      p$stmv_runmode$save_completed_data = TRUE      

  }


  # DATA inputs area created on-the-fly

  stmv( p=p )  # This will take from 40-70 hrs, depending upon system

# quick view
  predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
  statistics  = stmv_db( p=p, DS="stmv.stats" )
  locations   = spatial_grid( p )

  # comparisons
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

  # stats
  # statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
  statsvars = dimnames(statistics)[[2]]

  dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )

  dev.new(); levelplot( statistics[,match("nu", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
  dev.new(); levelplot( statistics[,match("sdSpatial", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
  dev.new(); levelplot( statistics[,match("localrange", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange

  # water only
  o = which( predictions>0 & predictions <500)
  #levelplot( log( statistics[o,match("sdTotal", statsvars)] ) ~ locations[o,1] + locations[o,2], aspect="iso" ) #sd total
   dev.new(); levelplot( log(predictions[o]) ~ locations[o,1] + locations[o,2], aspect="iso" )



# bring together stats and predictions and any other required computations: slope and curvature
# and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# .. still uses about 30-40 GB as the base layer is "superhighres" ..
# if parallelizing .. use different servers than local nodes
bathymetry_db( p=p, DS="complete.redo" ) # finalise at diff resolutions 15 min ..
bathymetry_db( p=p, DS="baseline.redo" )  # coords of areas of interest ..filtering of areas and or depth to reduce file size, in planar coords only




# a few plots :
pb = spatial_parameters( p=p, spatial_domain="canada.east" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )

# a few plots :
pb = spatial_parameters( p=p, spatial_domain="canada.east.highres" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )

pb = spatial_parameters( p=p, spatial_domain="canada.east.superhighres" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )


pb = spatial_parameters( p=p, spatial_domain="snowcrab" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )


pb = spatial_parameters( p=p, spatial_domain="SSE" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )



