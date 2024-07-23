
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
  

  # REQUIRED: this is necessary as data spans zero (-XXX to +XXX) .. 
  # naive lognormal of this would break due to infinity at log( x<= 0) .. THIS KEEPS DATA POSITIVE VALUED
  # or subset to positive-valued data o rmannual trnasform before and after ..
  p$data_transformation=list( forward=function(x){ x+2500 }, backward=function(x) {x-2500} )  

# run the model ... about 24 hrs depending upon number of posteriors to keep

  res = carstm_model( 
    p=p, 
    sppoly=sppoly,
    data='bathymetry_db( p=p, DS="carstm_inputs" )', 
    nposteriors = 1000,  # do not need too many as stmv solutions are default, this is to show proof of concept
    # redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
    # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
    # debug = TRUE,      
    # debug ="random_spatial",
    theta = c( 9.1236, 2.9961, 3.7036 ),
    toget = c("summary", "random_spatial", "predictions"),
    posterior_simulations_to_retain = c("predictions"),
    family = "lognormal",
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

      fit = carstm_model( p=p, DS="modelled_fit" )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
 
    }

# extract results and examine
  
  sppoly = areal_units( p=p )
  
  # posterior predictive check
  M = speciescomposition_db( p=p, DS='carstm_inputs', sppoly=sppoly  )
  carstm_posterior_predictive_check(p=p, M=M  )

  # EXAMINE POSTERIORS AND PRIORS
  res = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and summary

  names(res$hypers)
  for (i in 1:length(names(res$hypers)) ){
    o = carstm_prior_posterior_compare( hypers=res$hypers, all.hypers=res$all.hypers, vn=names(res$hypers)[i] )  
    dev.new(); print(o)
  }     


Time used:
    Pre = 40.3, Running = 197, Post = 67.3, Total = 304 
Fixed effects:
            mean sd 0.025quant 0.5quant 0.975quant mode kld
(Intercept) 7.93  0      7.929     7.93       7.93 7.93   0

Random effects:
  Name	  Model
    space BYM2 model

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

  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="random_spatial", probs=c(0.025, 0.975), 
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) 

  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="predictions", colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
    brks=seq(1, 501, 100) )

 
   
  
# end
 
