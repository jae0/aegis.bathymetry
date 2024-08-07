
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
    
      # if any new parameter settings are used for sppoly creation, 
      # then move them into temperature_parameters.R as the lookup mechanism 
      # uses these parameter settings for file lookup 
      sppoly = areal_units( p=p , redo=TRUE )  # this is the same as  aegis.polygons::01 polygons.R  
      plot( sppoly[ "AUID" ] )

      #NOTE:: as this is a lognormal model, all values < 1 is removed below
      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly , redo=TRUE )

      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly  )
      str(M)
    }
 
  sppoly = areal_units( p=p )  # this is the same as  aegis.polygons::01 polygons.R  
 
  p$space_name = sppoly$AUID 
  p$space_id = 1:nrow(sppoly)  # numst match M$space
  
 
# run the model ... about 24 hrs depending upon number of posteriors to keep

  carstm_model( 
    p=p, 
    sppoly=areal_units( p=p ),
    data='bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly )', 
    nposteriors = 1000,  # do not need too many as stmv solutions are default, this is to show proof of concept
    # redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
    # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
    # debug = TRUE,      
    theta = c( 8.988, 3.704, -2.970 ),
    toget = c("summary", "random_spatial", "predictions"),
    posterior_simulations_to_retain = c("predictions"),
    control.mode = list( restart=TRUE  ) ,
    control.inla = list( strategy="laplace", optimiser="gsl", restart=1 ),  # gsl = gsl::bfgs2
    # control.inla = list( strategy='laplace'),
    # control.inla = list( strategy='auto', int.strategy="eb" ),
    # control.inla = list( strategy='adaptive', int.strategy="eb" ),
    num.threads="4:2",  # very memory intensive ... serial process
    verbose=TRUE   
  ) 
  # return is modelinfo and input parameters p  (all saved to disk) for restarts
  
    # fixed and random effects are multiplicative effects 
    
    if (0) {
      # loading saved fit and results
      # very large files .. slow 

      fit = carstm_model( p=p, DS="modelled_fit" )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


      # EXAMINE POSTERIORS AND PRIORS
      all.hypers = INLA:::inla.all.hyper.postprocess(fit$all.hyper)
      hypers = fit$marginals.hyperpar
      names(hypers)

      carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, vn="Precision for space", transf=FALSE )  # no conversion to SD 
      carstm_prior_posterior_compare( hypers=hypers, all.hypers=all.hypers, vn="Phi for space" )  

      # posterior predictive check
      carstm_posterior_predictive_check(p=p, M=bathymetry_db( p=p, DS="carstm_inputs" )   )

    }

# extract results and examine
  
  sppoly = areal_units( p=p )
 
  smmy = carstm_model(  p=p, sppoly=sppoly, DS="carstm_modelled_summary" )  # parameters in p and direct summary
  smmy$direct
  Time used:
    Pre = 68.4, Running = 318, Post = 107, Total = 493 
Fixed effects:
            mean    sd 0.025quant 0.5quant 0.975quant mode kld
(Intercept) 7.93 0.001      7.928     7.93      7.931 7.93   0

Random effects:
  Name	  Model
    space BYM2 model

Model hyperparameters:
                                             mean    sd 0.025quant 0.5quant
Precision for the lognormal observations 9179.399 8.772   9159.719 9180.279
Precision for space                        48.832 0.438     47.852   48.874
Phi for space                               0.073 0.005      0.067    0.072
                                         0.975quant     mode
Precision for the lognormal observations   9193.556 9184.593
Precision for space                          49.541   49.091
Phi for space                                 0.085    0.068

Deviance Information Criterion (DIC) ...............: 27491869.50
Deviance Information Criterion (DIC, saturated) ....: 2903404.05
Effective number of parameters .....................: 29299.72

Watanabe-Akaike information criterion (WAIC) ...: 27286062.71
Effective number of parameters .................: 29314.60

Marginal log-Likelihood:  -13864228.09 
 is computed 

Model hyperparameters:
                                             mean    sd 0.025quant 0.5quant
Precision for the lognormal observations 9209.489 9.589   9187.494 9210.686
Precision for space                        19.824 0.211     19.495   19.799
Phi for space                               0.963 0.002      0.957    0.963
                                         0.975quant     mode
Precision for the lognormal observations   9224.153 9216.793
Precision for space                          20.305   19.686
Phi for space                                 0.966    0.965

Deviance Information Criterion (DIC) ...............: 27491807.87
Deviance Information Criterion (DIC, saturated) ....: 2912863.99
Effective number of parameters .....................: 29302.25

Watanabe-Akaike information criterion (WAIC) ...: 27284247.15
Effective number of parameters .................: 29419.03

Marginal log-Likelihood:  -13862283.54 
 
  # bbox = c(-71.5, 41, -52.5,  50.5 )
  additional_features = features_to_add( 
      p=p, 
      isobaths=c( 100, 200, 300, 400, 500  ), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) 
  )

  # maps of some of the results
  outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )

  
  # you might require the "grid"  library:  install.packages("grid")
  randomeffects = carstm_model(  p=p, sppoly=sppoly, DS="carstm_randomeffects" ) 
  modelinfo = carstm_model(  p=p, sppoly=sppoly, DS="carstm_modelinfo" ) 
  
  brks = seq(0, 600, 100)
  outfilename= file.path( outputdir, "bathymetry_spatialeffect_carstm.png")

  plt = carstm_map(  res=randomeffects, vn=c( "space", "re_total" ), 
#      breaks = brks, 
      modelinfo=modelinfo,
      title="Bathymetry random spatial (m)",
      colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
      additional_features=additional_features,
      outfilename=outfilename
    )
  plt


  predictions = carstm_model(  p=p, sppoly=sppoly, DS="carstm_predictions" ) 

  vn = "predictions"  
  # brks = pretty(  quantile( carstm_results_unpack( res, vn )[,"mean"], probs=c(0,0.975), na.rm=TRUE )  )
  brks = seq(0,600, 100)
  outfilename = file.path( outputdir, "bathymetry_predictions_carstm.png")

  plt = carstm_map( res=predictions, vn = vn,
    modelinfo=modelinfo,
    title="Bathymetry predicted (m)",
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
    additional_features=additional_features,
    outfilename=outfilename
  )
  plt

   
  
# end
 
