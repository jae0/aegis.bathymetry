
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


set.seed(12345)


# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
  p = aegis.bathymetry::bathymetry_parameters( project_class="carstm", areal_units_resolution_km = 5 )  # defaults are hard coded as a lattice .. anything else takes a very long time
 

    if (0) {
      # to recreate the underlying data:
      xydata=bathymetry_db(p=p, DS="areal_units_input", redo=TRUE)
    
      sppoly = areal_units( p=p , redo=TRUE )  # this is the same as  aegis.polygons::01 polygons.R  
      plot( sppoly[ "AUID" ] )

      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly , redo=TRUE )
      
      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly  )
      str(M)

    }
 
  sppoly = areal_units( p=p )  # this is the same as  aegis.polygons::01 polygons.R  
 
  p$space_name = sppoly$AUID 
  p$space_id = 1:nrow(sppoly)  # numst match M$space

# run the model ... about 24 hrs depending upon number of posteriors to keep

  res = carstm_model( 
    p=p, 
    sppoly=areal_units( p=p ),
    data='bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly )', 
    nposteriors = 1000,  # do not need too many as stmv solutions are default, this is to show proof of concept
    # redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
    # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
    # debug = TRUE,      
    theta = c( 4.1828, -0.9646, 3.6723 ),
    toget = c("summary", "random_spatial", "predictions"),
    posterior_simulations_to_retain = c("predictions"),
    control.mode = list( restart=TRUE  ) ,
    # control.inla = list( strategy='laplace'),
    # control.inla = list( strategy='adaptive', int.strategy="eb" ),
    # control.inla = list( strategy='adaptive', int.strategy="eb" ),
    num.threads="4:2",  # very memory intensive ... serial process
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
  
  sppoly = areal_units( p=p )
 
  smmy = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and direct summary
  smmy$direct 


  # bbox = c(-71.5, 41, -52.5,  50.5 )
  additional_features = features_to_add( 
      p=p, 
      isobaths=c( 100, 200, 300, 400, 500  ), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) 
  )

# maps of some of the results
  outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )
  carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="random_spatial", probs=c(0.025, 0.975) ) 


  vn = "predictions"  
  # brks = pretty(  quantile( carstm_results_unpack( res, vn )[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )
  brks = seq(0,600, 100)
  
  outfilename = file.path( outputdir, "bathymetry_predictions_carstm.png")

  plt = carstm_map( res=predictions, vn = vn,
    modelinfo=modelinfo,
    breaks=brks,
    title="Bathymetry predicted (m)",
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
    additional_features=additional_features,
    outfilename=outfilename
  )
  plt

   
  
# end
 
