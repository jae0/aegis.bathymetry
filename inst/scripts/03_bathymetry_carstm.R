
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

if(0) {
      p$fraction_todrop = 1/4 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
      p$fraction_cv = 1.0  #sd/mean no.
      p$fraction_good_bad = 0.9
      p$areal_units_constraint_nmin = 1000  # length(p$yrs)
      p$nAU_min = 100
}
  # to recreate the underlying data
  # xydata=bathymetry_db(p=p, DS="areal_units_input", redo=TRUE)

  sppoly = areal_units( p=p , redo=T )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
  plot( sppoly[ "AUID" ] )



# run the model ... about 24 hrs
  fit = carstm_model( p=p, M='bathymetry_db( p=p, DS="carstm_inputs" )' ) # run model and obtain predictions

# loading saved fit and results
  fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  res = carstm_summary( p=p ) # to load currently saved sppoly

  plot(fit)
  plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
  s = summary(fit)

  s$dic$dic # 1404489
  s$dic$p.eff # 151412

# maps of some of the results
  vn = paste(p$variabletomodel, "predicted", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

# end

  if (0) {
    # example work flow for testing alternate params .. large grid is for debugging
      p = aegis.bathymetry::bathymetry_parameters(
        project_class="carstm",
        areal_units_resolution_km=100,
        carstm_inputs_aggregated = TRUE
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
      fit = carstm_model( p=p, M=M ) # run model and obtain predictions
      # fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit

      res = carstm_summary( p=p ) # to load currently saved sppoly

      vn = paste(p$variabletomodel, "predicted", sep=".")
      carstm_plot( p=p, res=res, vn=vn )
  }


