
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
    inla.setOption(num.threads= floor( parallel::detectCores() / 3 ) )
    inla.setOption(blas.num.threads= 3 )
 

    if (0) {
      # to recreate the underlying data:
      xydata=bathymetry_db(p=p, DS="areal_units_input", redo=TRUE)

      sppoly = areal_units( p=p , redo=T )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      plot( sppoly[ "AUID" ] )
    }


# run the model ... about 24 hrs
  fit = carstm_model( p=p, M='bathymetry_db( p=p, DS="carstm_inputs" )' ) # run model and obtain predictions

    if (0) {
      # loading saved fit and results
      # very large files .. slow 
      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }

# extract results and examine
  res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  
    res$summary$dic$dic
    res$summary$dic$p.eff
    res$dyear



  plot_crs = p$aegis_proj4string_planar_km
  coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
  isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs )


# maps of some of the results
  vn = paste(p$variabletomodel, "predicted", sep=".")

  carstm_map(  res=res, vn=vn, 
      breaks =pretty(p$discretization$z),
      palette="viridis",
      coastline=coastline,
      isobaths=isobaths,
      main="Bathymetry random unstructured" 
  )


  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_map(  res=res, vn=vn, 
     palette="viridis",
      coastline=coastline,
      isobaths=isobaths,
      main="Bathymetry random iid" 
  )


  vn = paste(p$variabletomodel, "random_space_nonspatial", sep=".")
  carstm_map(  res=res, vn=vn, 
      palette="viridis",
      coastline=coastline,
      isobaths=isobaths,
      main="Bathymetry random unstructured" 
  )


  vn = paste(p$variabletomodel, "random_space_spatial", sep=".")
  carstm_map(  res=res, vn=vn, 
      palette="viridis",
      coastline=coastline,
      isobaths=isobaths,
      main="Bathymetry spatially structured" 
  )

# end

  if (0) {
    # example work flow for testing alternate params .. large grid is for debugging
      p = aegis.bathymetry::bathymetry_parameters(
        project_class="carstm",
        areal_units_resolution_km=100,
        carstm_inputs_prefilter = "sampled",
        carstm_inputs_prefilter_n = 10
      )

    # example sequence to force creating of input data for modelling
      sppoly = areal_units( p=p, redo=TRUE );
      plot(sppoly[, "AUID"])
      for( i in 1:3) plot( as(p$coastLayout[[i]][[2]], "sf"), add=TRUE )

      # set up default map projection
      oo = aegis.coastline::coastline_layout( p=p )
      p = parameters_add_without_overwriting( p,
        coastLayout = oo[["coastLayout"]],
        bounding_domain = oo[["bounding_domain"]]
      )
      oo = NULL

      B = bathymetry_db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
      str(B)
      B= NULL ;  gc()

      M = bathymetry_db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
      str(M)
      fit = carstm_model( p=p, M=M, DS="redo" ) # run model and obtain predictions
      # fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit

      res = carstm_model( p=p, DS="carstm_modelled_summary" ) # to load currently saved sppoly

      vn = paste(p$variabletomodel, "predicted", sep=".")
      carstm_map(  res=res, vn=vn, time_match=time_match, 
            breaks =pretty(p$discretization$z),
            palette="viridis",
            coastline=coastline,
            isobaths=isobaths,
            main= "Depth" 
      )
    }


