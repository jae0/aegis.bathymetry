
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
options(sf_use_s2 = FALSE)  # seems to cause some problems ... in areal_units()

# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
  p = aegis.bathymetry::bathymetry_parameters( project_class="carstm", areal_units_resolution_km = 5 )  # defaults are hard coded as a lattice .. anything else takes a very long time
 
  # bbox = c(-71.5, 41, -52.5,  50.5 )
  additional_features = features_to_add( 
      p=p, 
      isobaths=c( 100, 200, 300, 400, 500  ), 
      xlim=c(-80,-40), 
      ylim=c(38, 60) 
  )


    if (0) {
      # to recreate the underlying data:
      xydata=bathymetry_db(p=p, DS="areal_units_input", redo=TRUE)
    
      # if any new parameter settings are used for sppoly creation, 
      # then move them into temperature_parameters.R as the lookup mechanism 
      # uses these parameter settings for file lookup 
      sppoly = areal_units( p=p , redo=TRUE )  # this is the same as  aegis.polygons::01 polygons.R  
    
      plt = areal_units( sppoly=sppoly, xydata=xydata, additional_features=additional_features, plotit=TRUE )
      
      (plt)

  
      #NOTE:: as this is a lognormal model, all values < 1 is removed below
      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly , redo=TRUE )

      M = bathymetry_db( p=p, DS="carstm_inputs", sppoly=sppoly  )
      str(M)
    }
 
  sppoly = areal_units( p=p )  # this is the same as  aegis.polygons::01 polygons.R  
 
  p$space_name = sppoly$AUID 
  p$space_id = 1:nrow(sppoly)  # numst match M$space
  
 
# run the model ... about 5 min (fit) depending upon number of posteriors to keep

  carstm_model( 
    p=p, 
    sppoly=sppoly,
    data='bathymetry_db( p=p, DS="carstm_inputs" )', 
    # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
    # debug = TRUE,      
    # debug ="random_spatial",
    theta = c( 3.9725, 2.1879, 4.1276 ),  #  3.9729 2.1851 3.8436   
    toget = c("summary", "random_spatial", "predictions"),
    # nposteriors = 0,  # do not need samples ut this is where you need to specify along with posterior_simulations_to_retain
    # posterior_simulations_to_retain = c("predictions"),
    family = "lognormal",
    # control.mode = list( restart=TRUE  ) ,
    # control.inla = list( strategy="laplace", optimiser="gsl", restart=1 ),  # gsl = gsl::bfgs2
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
 
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

			# posterior predictive check
			M = bathymetry_db( p=p, DS='carstm_inputs', sppoly=sppoly  )
			carstm_posterior_predictive_check(p=p, M=M  )


      # EXAMINE POSTERIORS AND PRIORS
      res = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and summary

      outputdir = file.path(p$modeldir, p$carstm_model_label)
      
      res_vars = c( names( res$hypers), names(res$fixed) )
      for (i in 1:length(res_vars) ) {
        o = carstm_prior_posterior_compare( res, vn=res_vars[i], outputdir=outputdir )  
        dev.new(); print(o)
      }     
 

    }

# extract results and examine
  
  sppoly = areal_units( p=p )
  

 
  # maps of some of the results
  outputdir = file.path(p$modeldir, p$carstm_model_label, "maps" )

  # random effects are on response scale 
  carstm_plot_map( p=p, outputdir=outputdir, additional_features=additional_features, 
    toplot="random_spatial", probs=c(0.025, 0.975), transf=log10,  # log10 transform for map
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) 
 
  # predictions are on response scale
  carstm_plot_map( p=p, outputdir=outputdir,  
    toplot="predictions", transf=log10,  # log10 transform for map
    additional_features=additional_features, 
    colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) )
 
 
  # more direct control over map
  # random effects  ..i.e.,  deviation from lognormal model ( pure spatial effect )
    res = carstm_model(  p=p, DS="carstm_randomeffects" )  

    outfilename= file.path( outputdir, paste("depth_spatialeffect_carstm", "png", sep=".") )
    plt = carstm_map(  res=res, vn= c(  "space", "re_total" ), 
        sppoly=sppoly,
        transformation=log10,
        title="Depth spatial error (log 10 m)",
        colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
        additional_features=additional_features,
        outfilename=outfilename
    )  
    plt
  
 
   
  
# end
 
