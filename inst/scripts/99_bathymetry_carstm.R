

### this works fine, but resolutions less than 5 km are very slow and occasionally fails to converge
### it is here as an example


### WARNING::  running carstm directly... < 5 km resokution becomes rather slow ...

## SLOW :: about 1 hr for each config x 25 configs .. ie., 24 hrs .. consider removing configs if no need for posterior samples

# -- might need to run in a shell if number of polygons are large:  ulimit -s 16384


# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
# --- look inside "parameters_default" and define alternates based upon it

p = aegis.bathymetry::bathymetry_parameters( project_class="carstm" )  # carstm on 25 km lattice grid


    # DS = "parameters_default"; areal_units_resolution_km=5 ... takes 79 Hrs!

# example sequence to force creating of input data for modelling
  sppoly = areal_units( p=p, redo=TRUE ); plot(sppoly[, "AUID"]) # or: spplot( sppoly, "AUID", main="AUID", sp.layout=p$coastLayout )
  M = bathymetry_db( p=p, DS="aggregated_data" , redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = bathymetry_db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  str(M)
  M = NULL; gc()

# run the model ... about 24 hrs
  fit = carstm_model( p=p, M='bathymetry_db( p=p, DS="carstm_inputs" )' ) # run model and obtain predictions

# loading saved results
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
