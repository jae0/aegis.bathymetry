
bathymetry_lookup = function( p, locs, vnames=NULL, output_data_class="points", source_data_class="aggregated_rawdata" ){

  # locs = output point locations or sp_polygon

  require(aegis.bathymetry)

  if ( output_data_class=="points") {

    pAg = spatial_parameters( spatial_domain=p$spatial_domain )

    if (!exists("inputdata_spatial_discretization_planar_km", pAg)) {
      if (!exists("inputdata_spatial_discretization_planar_km", p)) {
        pAg$inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km
      } else {
        pAg$inputdata_spatial_discretization_planar_km = 1
      }
    }
    if (!exists("variabletomodel", pAg)) pAg$variabletomodel = "z"

    if (source_data_class=="rawdata") {
      B = bathymetry.db ( p=pAg, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
    } else if (source_data_class=="aggregated_rawdata") {
      B = bathymetry.db ( p=pAg, DS="aggregated_data" )
    } else if (source_data_class=="modelled_stmv") {
      pA = bathymetry_parameters(p=pAg, project_class="stmv")
      B = bathymetry.db(p=pAg, DS="baseline", varnames="all" )
    }

    B_varnames = names(B)
    if ( !is.null(vnames) ) {
      if ( !all(vnames %in% B_varnames)) stop("vnames not found:", vnames, "spatial_domain:", B_varnames)
    } else {
      vnames = setdiff( B_varnames, c( "plon", "plat", "lon", "lat") )
    }

    B = lonlat2planar( B, proj.type=pAg$aegis_proj4string_planar_km )
    B$plon = round(B$plon / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
    B$plat = round(B$plat / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
    B_map = paste(B$plon, B$plat, sep=".")

    locs = lonlat2planar( locs, proj.type=pAg$aegis_proj4string_planar_km )
    locs$plon = round(locs$plon / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
    locs$plat = round(locs$plat / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
    locs_map = paste(locs$plon, locs$plat, sep=".")

    locs_index = match( locs_map, B_map )
    out = B[locs_index, vnames]

    return(out)
  }


  if ( output_data_class=="areal_units") {

    pAg = bathymetry_carstm( p=p, DS="parameters", variabletomodel="z" )

    # locs must be sppoly with "AUID"
    require(raster)
    raster_template = raster(extent(locs)) # +1 to increase the area
    res(raster_template) = pAg$areal_units_resolution_km  # in meters
    crs(raster_template) = projection(locs) # transfer the coordinate system to the raster

    if (source_data_class=="rawdata") {
      B = bathymetry.db ( p=pAg, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
      B = sf::st_as_sf( B, coords=c("plon", "plat") )
    } else if (source_data_class=="aggregated_rawdata") {
      B = bathymetry.db ( p=pAg, DS="aggregated_data" )
      B = sf::st_as_sf( B, coords=c("plon", "plat") )
    } else if (source_data_class=="modelled_stmv") {
      pAg = bathymetry_parameters(p=p, project_class="stmv")
      B = bathymetry.db(p=pAg, DS="baseline", varnames="all" )
      B = sf::st_as_sf( B, coords=c("plon", "plat") )
      
    } else if (source_data_class=="modelled_carstm") {
      # copy of param list for global analysis in aegis.bathymetry/inst/scripts/02.bathymetry.carstm.R
      pAg = aegis.bathymetry::bathymetry_carstm(
        DS = "parameters",
        project_name = "bathymetry",
        spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
        variabletomodel ="z",
        carstm_model_label = "production",
        inputdata_spatial_discretization_planar_km = 1,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
        areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
        areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
        areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
        areal_units_overlay = "none"
      )

      B = carstm_summary( p=pAg ) # to load currently saved sppoly  ("locs")
      B_sppoly = areal_units( p=pAg )
      B_sppoly$z = NA
      B_sppoly$z = B$z.predicted[ match( B_sppoly$AUID, B$AUID )] ## <<< check order in case a proper merge is required
      B = as(B_sppoly, "sf")
    }

    B_varnames = names(B)
    if ( !is.null(vnames) ) {
      if ( !all(vnames %in% B_varnames)) stop("vnames not found:", vnames, "spatial_domain:", B_varnames)
    } else {
      vnames = setdiff( B_varnames, c( "plon", "plat", "lon", "lat") )
    }
    

    for (vn in vnames) {
      Bf = fasterize::fasterize( B), raster_template, field=vn )

      vn1 = paste(vn, "mean", sep="." )
      locs[vn1] = raster::extract( Bf, locs, fun=mean, na.rm=TRUE)

      vn2 = paste(vn, "sd", sep="." )
      locs[vn2] = raster::extract( Bf, locs, fun=sd, na.rm=TRUE)
    }

    out = locs[, setdiff( names(locs), c( "plon", "plat", "lon", "lat") ) ]

  }

  return(out)
}
