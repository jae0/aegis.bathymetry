
## SLOW :: about 1 hr for each config x 25 configs .. ie., 24 hrs .. consider removing configs if no need for posterior samples

# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
  p = aegis.bathymetry::bathymetry_carstm(
    DS = "parameters",
    project_name = "bathymetry",
    project_class = "carstm", # defines which parameter class / set to load
    spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
    inputdata_spatial_discretization_planar_km = 1,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
    areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
  )

# example sequence to force creating of input data for modelling
  sppoly = areal_units( p=p, redo=TRUE ); plot(sppoly) # or: spplot( sppoly, "StrataID", main="StrataID", sp.layout=p$coastLayout )
  M = bathymetry_carstm( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M = bathymetry_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  str(M)

# run the model ... about 24 hrs
  res = carstm_model( p=p, M=bathymetry_cartsm( p=p, DS="carstm_inputs" ) ) # run model and obtain predictions

# loading saved results
  res = carstm_model( p=p, DS="carstm_modelled" ) # to load currently saved sppoly
  fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit

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

  vn = paste(p$variabletomodel, "random_strata_nonspatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

  vn = paste(p$variabletomodel, "random_strata_spatial", sep=".")
  carstm_plot( p=p, res=res, vn=vn )

# end
