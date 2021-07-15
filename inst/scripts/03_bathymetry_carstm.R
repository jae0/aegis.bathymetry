
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
  p = aegis.bathymetry::bathymetry_parameters( project_class="carstm" )  # defaults are hard coded

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
# careful: any reasonable compression can increase save time by hours ..
  res = carstm_model( p=p, M='bathymetry_db( p=p, DS="carstm_inputs" )', compression_level=0, quantile_bounds=c(0,1), nposteriors=1000, redo_fit=TRUE ) 
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

  tmout = carstm_map(  
    res=res,    
    vn = "predictions" ,
    palette="-Spectral",
    plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
    outfilename= file.path( project.datadirectory("aegis", "bathymetry", "maps"), "bathymetry_predictions_carstm.png"),
    main="Bathymetry predicted" ,
    tmap_zoom=6.5
  )
  tmout

  tmout = carstm_map(  
    res=res, 
    vn = c( "random", "space", "combined" ), 
    palette="-Spectral",
    plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
    outfilename= file.path( project.datadirectory("aegis", "bathymetry", "maps"), "bathymetry_spatialeffect_carstm.png"),
    main="Bathymetry random spatial" ,
    tmap_zoom=6.5
  )
  


# end
 
