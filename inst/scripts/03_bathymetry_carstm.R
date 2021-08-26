
### --------------------------------------------------------------------------------------------------------
### this works fine, but resolutions less than 5 km are very slow and occasionally fails to converge
### it is here as an example
### the stmv-based results are better in that spatial resolution is higher (0.5 km or better),
### permitting more useful higher order descriptions (slope, curviture) and variability stats.

### WARNING::  running carstm directly on  raw data and resolution < 5 km  becomes very slow ...
### instead use the aggregated data and larger resolution.

### IF number of polygons are large, this can fail (due to INLA representation) a potential solution
### is to run the following in a shell:  ulimit -s 16384
### --------------------------------------------------------------------------------------------------------



# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
  p = aegis.bathymetry::bathymetry_parameters( project_class="carstm", areal_units_resolution_km = 5 )  # defaults are hard coded

    # adjust based upon RAM requirements and ncores
    # inla.setOption(num.threads= floor( parallel::detectCores() / 3 ) )
    # inla.setOption(blas.num.threads= 3 )
 

    if (0) {
      # to recreate the underlying data:
      xydata=bathymetry_db(p=p, DS="areal_units_input", redo=TRUE)

      sppoly = areal_units( p=p , redo=T )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      plot( sppoly[ "AUID" ] )
      
      bathymetry_db( p=p, DS="carstm_inputs", redo=TRUE )
    }


# run the model ... about 24 hrs

  res = carstm_model( 
    p=p, 
    data='bathymetry_db( p=p, DS="carstm_inputs" )', 
    num.threads="4:2",
    compress="xz", 
    control.inla = list( strategy='adaptive', int.strategy='eb' ),
    redo_fit=TRUE, 
    verbose=TRUE   
  ) 

    # run model and obtain predictions, 0== no file compression
    # fixed and random effects are multiplicative effects 
    # quantile_bounds=c(0,1) means do not extrapolate
    
    if (0) {
      # loading saved fit and results
      # very large files .. slow 

      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }

# extract results and examine
  res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  res$summary  


# maps of some of the results
  outputdir = file.path( gsub( ".rdata", "", carstm_filenames(p, returntype="carstm_modelled_fit") ), "figures" )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  tmout = carstm_map( res=res, vn = "predictions",
    title="Bathymetry predicted (m)",
    palette="-Spectral",
    plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
    outfilename= file.path( outputdir, "bathymetry_predictions_carstm.png"),
    tmap_zoom= c((p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8, 6.5)
  )
  tmout

# random effects  ..i.e.,  deviation from lognormal model
  tmout = carstm_map( res=res, vn = c( "random", "space", "combined" ), 
    title="Bathymetry random spatial (m)",
    palette="-Spectral",
    plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
    outfilename= file.path( outputdir, "bathymetry_spatialeffect_carstm.png"),
    tmap_zoom= c((p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8, 6.5)
  )
  tmout




# end
 
