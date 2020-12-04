
# Bathymetry

# Spatial interpolation using stmv on a lattice grid .. continuous form

# Total "superhighres": 2-5 GB/process and 4 GB in parent for fft
# gam method requires more ~ 2X
# boundary def takes too long .. too much data to process -- skip
# "highres": ~ 20 hr with 8, 3.2 Ghz cpus on thoth using fft method jc: 2016 or 2~ 6 hr on hyperion
# "superhighres" fft: looks to be the best in performance/quality; req ~5 GB per process req
# FFT is the method of choice for speed and ability to capture the variability
# krige method is a bit too oversmoothed, especially where rapid changes are occuring

#  ~40 hrs to scale
scale_ram_required_main_process = 20 # GB twostep / fft
scale_ram_required_per_process  = 6 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB)
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

# ~96 hrs
interpolate_ram_required_main_process = 10 # GB twostep / fft
interpolate_ram_required_per_process  = 8 # twostep / fft /fields vario ..
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )


p = aegis.bathymetry::bathymetry_parameters(
  project_class="stmv",
  stmv_runmode = list(
    scale = rep("localhost", scale_ncpus),
    interpolate = list(
      c1 = rep("localhost", interpolate_ncpus),  # ncpus for each runmode
      c2 = rep("localhost", interpolate_ncpus),  # ncpus for each runmode
      c3 = rep("localhost", max(1, interpolate_ncpus-1)),
      c4 = rep("localhost", max(1, interpolate_ncpus-1)),
      c5 = rep("localhost", max(1, interpolate_ncpus-2))
    ),
    # if a good idea of autocorrelation is missing, forcing via explicit distance limits is an option
    # interpolate_distance_basis = list(
    #   d1 = rep("localhost", interpolate_ncpus),
    #   d2 = rep("localhost", interpolate_ncpus),
    #   d3 = rep("localhost", max(1, interpolate_ncpus-1)),
    #   d4 = rep("localhost", max(1, interpolate_ncpus-1)),
    #   d5 = rep("localhost", max(1, interpolate_ncpus-2)),
    #   d6 = rep("localhost", max(1, interpolate_ncpus-2))
    # ),
    interpolate_predictions = list(
      c1 = rep("localhost", max(1, interpolate_ncpus-1)),  # ncpus for each runmode
      c2 = rep("localhost", max(1, interpolate_ncpus-1)),  # ncpus for each runmode
      c3 = rep("localhost", max(1, interpolate_ncpus-2)),
      c4 = rep("localhost", max(1, interpolate_ncpus-3)),
      c5 = rep("localhost", max(1, interpolate_ncpus-4)),
      c6 = rep("localhost", max(1, interpolate_ncpus-4)),
      c7 = rep("localhost", max(1, interpolate_ncpus-5))
     ),
    globalmodel = FALSE,
    # restart_load = "interpolate_correlation_basis_0.01" ,  # only needed if this is restarting from some saved instance
    save_intermediate_results = TRUE,
    save_completed_data = TRUE

  )  # ncpus for each runmode
)



if (0) {  # model testing
  # if resetting data for input to stmv run this or if altering discretization resolution
  bathymetry_db( p=p, DS="stmv_inputs_redo" )  # recreate fields for .. requires 60GB+
  # p$restart_load = paste("interpolate_correlation_basis_", p$stmv_autocorrelation_basis_interpolation[length(p$stmv_autocorrelation_basis_interpolation)], sep="")  # to choose the last save
}

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
o = which( predictions>0 & predictions <1000)
#levelplot( log( statistics[o,match("sdTotal", statsvars)] ) ~ locations[o,1] + locations[o,2], aspect="iso" ) #sd total
levelplot( log(predictions[o]) ~ locations[o,1] + locations[o,2], aspect="iso" )



# bring together stats and predictions and any other required computations: slope and curvature
# and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# .. still uses about 30-40 GB as the base layer is "superhighres" ..
# if parallelizing .. use different servers than local nodes
bathymetry_db( p=p, DS="complete.redo" ) # finalise at diff resolutions 15 min ..
bathymetry_db( p=p, DS="baseline.redo" )  # coords of areas of interest ..filtering of areas and or depth to reduce file size, in planar coords only




# a few plots :
pb = aegis.bathymetry::bathymetry_parameters( project_class="stmv", spatial_domain="canada.east.highres" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )

pb = aegis.bathymetry::bathymetry_parameters( project_class="stmv", spatial_domain="canada.east.superhighres" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )


pb = aegis.bathymetry::bathymetry_parameters( project_class="stmv", spatial_domain="snowcrab" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )


pb = aegis.bathymetry::bathymetry_parameters( project_class="stmv", spatial_domain="SSE" )
bathymetry_figures( p=pb, varnames=c("z", "dZ", "ddZ", "b.localrange"), logyvar=TRUE, savetofile="png" )
bathymetry_figures( p=pb, varnames=c("b.sdTotal", "b.sdSpatial", "b.sdObs"), logyvar=FALSE, savetofile="png" )



