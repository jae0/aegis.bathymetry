
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
  p = aegis.bathymetry::bathymetry_parameters( project_class="carstm", areal_units_resolution_km = 5 )  # defaults are hard coded as a lattice .. anything else takes a very long time
 

    if (0) {
      # to recreate the underlying data:
      xydata=bathymetry_db(p=p, DS="areal_units_input", redo=TRUE)
    
      sppoly = areal_units( p=p , redo=TRUE )  # this is the same as  aegis.polygons::01 polygons.R  
      plot( sppoly[ "AUID" ] )

      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=areal_units( p=p ) , redo=TRUE )
      
      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=areal_units( p=p )  )
      str(M)

    }


# run the model ... about 24 hrs depending upon no posteriors to keep

  res = carstm_model( 
    p=p, 
    sppoly=areal_units( p=p ),
    data='bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly )', 
    nposteriors = 1000,  # do not need too many as stmv solutions are default
    # redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
    # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
    # debug = TRUE,
    control.mode = list( restart=FALSE, theta= c( 8.988, 3.704, -2.970 ) ) ,
    # control.inla = list( strategy='laplace'),
    # control.inla = list( strategy='adaptive', int.strategy="eb" ),
    # control.inla = list( strategy='adaptive', int.strategy="eb" ),
    num.threads="1:1",  # very memory intensive ... serial process
    verbose=TRUE   
  ) 

  

    # run model and obtain predictions, 0== no file compression
    # fixed and random effects are multiplicative effects 
    # quantile_bounds=c(0,1) means do not extrapolate
    
    if (0) {
      # loading saved fit and results
      # very large files .. slow 

      fit = carstm_model( p=p, DS="carstm_modelled_fit", sppoly=sppoly )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }

# extract results and examine
  
  sppoly = areal_units( p=p )

  res = carstm_model( p=p, DS="carstm_modelled_summary", sppoly=sppoly ) # to load currently saved results
  res$summary  


  # bbox = c(-71.5, 41, -52.5,  50.5 )
  additional_features = additional_features_tmap( 
      p=p, 
      isobaths=c( 10, 100, 200, 300, 500, 1000 ), 
      coastline =  c("canada"), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) 
  )

# maps of some of the results
  outputdir = file.path(p$data_root, "maps", p$carstm_model_label )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

  
  vn = "predictions"  
  brks = pretty(  quantile( carstm_results_unpack( res, vn )[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )
  outfilename = file.path( outputdir, "bathymetry_predictions_carstm.png")

  tmout = carstm_map( res=res, vn = vn,
    breaks = brks, 
    title="Bathymetry predicted (m)",
    palette="-Spectral",
    plot_elements=c(  "compass", "scale_bar", "legend" ),
    additional_features=additional_features,
    outfilename=outfilename
  )
  tmout


  

# random effects  ..i.e.,  deviation from lognormal model
  vn = c( "random", "space", "combined" )
  oo = carstm_results_unpack( res, vn )[,"mean"]
  oo = oo[oo< 5000]
  brks = pretty(  quantile( oo, probs=c(0,0.975), na.rm=TRUE )  )

  outfilename= file.path( outputdir, "bathymetry_spatialeffect_carstm.png")

  tmout = carstm_map( res=res, vn=vn, 
    breaks = brks, 
    title="Bathymetry random spatial (m)",
    palette="-Spectral",
    plot_elements=c(  "compass", "scale_bar", "legend" ),
    additional_features=additional_features,
    outfilename=outfilename
  )
  tmout

 

  predictions_errors_removed = posterior_summary( res$sims$predictions - res$sims$space$combined)

  outfilename= file.path( outputdir, "bathymetry_denoised_carstm.png")
  brks = pretty(  quantile( carstm_results_unpack( res, vn )[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )

  sppoly$z_denoised = predictions_errors_removed$mean
  tmout = carstm_map( vn="z_denoised", 
    sppoly = sppoly,
    breaks = brks, 
    title="Bathymetry denoised (m)",
    palette="-Spectral",
    plot_elements=c(  "compass", "scale_bar", "legend" ),
    additional_features=additional_features,
    outfilename=outfilename
  )
  tmout

# end
 
