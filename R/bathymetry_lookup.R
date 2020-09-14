bathymetry_lookup = function( p, locs, vnames=NULL, data_class="rawdata" ){

  require(aegis.bathymetry)
  if ( data_class=="rawdata") {
    pB = spatial_parameters( spatial_domain=p$spatial_domain )
    if (!exists("inputdata_spatial_discretization_planar_km", pB)) {
      if (!exists("inputdata_spatial_discretization_planar_km", p)) {
        pB$inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km
      } else {
        pB$inputdata_spatial_discretization_planar_km = 1
      }
    }
    if (!exists("variabletomodel", pB)) pB$variabletomodel = "z"
    B = bathymetry.db ( p=pB, DS="aggregated_data" )
    B_varnames = names(B)
    if ( !is.null(vnames) ) {
      if ( !all(vnames %in% B_varnames)) stop("vnames not found:", vnames, "spatial_domain:", B_varnames)
    } else {
      vnames = setdiff( B_varnames, c( "plon", "plat") )
    }


    B = lonlat2planar( B, proj.type=pB$aegis_proj4string_planar_km )
    B$plon = round(B$plon / pB$inputdata_spatial_discretization_planar_km + 1 ) * pB$inputdata_spatial_discretization_planar_km
    B$plat = round(B$plat / pB$inputdata_spatial_discretization_planar_km + 1 ) * pB$inputdata_spatial_discretization_planar_km

    locs = lonlat2planar( locs, proj.type=pB$aegis_proj4string_planar_km )
    locs$plon = round(locs$plon / pB$inputdata_spatial_discretization_planar_km + 1 ) * pB$inputdata_spatial_discretization_planar_km
    locs$plat = round(locs$plat / pB$inputdata_spatial_discretization_planar_km + 1 ) * pB$inputdata_spatial_discretization_planar_km
    locs_map = paste(locs$plon, locs$plat, sep=".")

    B_map = paste(B$plon, B$plat, sep=".")
    locs_index = match( locs_map, B_map )
    out = B[locs_index, vnames]
    return(out)

  }

  if ( data_class=="stmv") {

    pB = bathymetry_parameters(p=p, project_class="stmv")
    B = bathymetry.db(p=pB, DS="baseline" )
    B_map = stmv::array_map( "xy->1", B[,c("plon","plat")], gridparams=pB$gridparams )
    B_varnames = names(B)
    if ( !is.null(vnames) ) {
      if ( !all(vnames %in% B_varnames)) stop("vnames not found:", vnames, "spatial_domain:", B_varnames)
    } else {
      vnames = setdiff( B_varnames, c( "plon", "plat") )
    }
    locs_map = stmv::array_map( "xy->1", locs, gridparams=pB$gridparams )
    locs_index = match( locs_map, B_map )
    out = B[locs_index, vnames]

  } else if ( data_class=="carstm_global") {

    pB = bathymetry_carstm(
      DS = "parameters",
      project_name = "bathymetry",
      spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
      variabletomodel ="z",
      carstm_model_label = "production",
      inputdata_spatial_discretization_planar_km = 1,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_overlay = "none"
    )
    res = carstm_summary( p=pB ) # to load currently saved sppoly

  } else if ( data_class=="carstm_local") {

    pB = bathymetry_carstm( p=p, DS="parameters", variabletomodel="z" )
    if (!(exists(pB$variabletomodel, M ))) M[,pB$variabletomodel] = NA
    kk =  which( !is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0) {
      AD = bathymetry.db ( p=pB, DS="aggregated_data"   )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( AD[, paste(pB$variabletomodel, "mean", sep="." )], AD$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      M[kk, pB$variabletomodel] = oo[jj ]
    }

  }

  return(out)
}
