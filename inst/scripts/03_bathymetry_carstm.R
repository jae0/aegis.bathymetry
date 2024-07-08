
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
    data='bathymetry_db( p=p, DS="carstm_inputs" )', 
    nposteriors = 1000,  # do not need too many as stmv solutions are default, this is to show proof of concept
    # redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
    redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
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

      carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, vn="Precision for space", transf=FALSE )  # no conversion to SD 
      carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, vn="Phi for space" )  

      # posterior predictive check
      carstm_posterior_predictive_check(p=p, M=substrate_db( p=p, DS="carstm_inputs" ), transf=TRUE  )

    }

# extract results and examine
  
  sppoly = areal_units( p=p )
 
  smmy = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and direct summary
  smmy$direct 

Deviance Information Criterion (DIC) ...............: 24755267.86
Deviance Information Criterion (DIC, saturated) ....: 2843680.48
Effective number of parameters .....................: 2288.50

Watanabe-Akaike information criterion (WAIC) ...: 26358325.05
Effective number of parameters .................: 881372.01

Marginal log-Likelihood:  -12506191.81 

  # bbox = c(-71.5, 41, -52.5,  50.5 )
  additional_features = features_to_add( 
      p=p, 
      isobaths=c( 100, 200, 300, 400, 500  ), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) 
  )

  # maps of some of the results
  outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )

  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="random_spatial", probs=c(0.025, 0.975),  transf=log10, 
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) 

  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="predictions", colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) )

 
   
  
# end
 
