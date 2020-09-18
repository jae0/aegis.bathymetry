
bathymetry_lookup = function( p, locs, vnames=NULL, output_data_class="points", source_data_class="aggregated_rawdata" ){

  # locs = output point locations or sp_polygon
  out = NULL

  require(aegis.bathymetry)

  # set up parameters for input data
  if ( source_data_class %in% c("rawdata", "aggregated_rawdata", "modelled_stmv" ) ) {
      pAg = spatial_parameters( spatial_domain=p$spatial_domain )

      if (source_data_class=="modelled_stmv") {

        pAg = bathymetry_parameters(p=p, project_class="stmv")

      } else {
        # minimal info required ...
        if (!exists("inputdata_spatial_discretization_planar_km", pAg)) {
          if (!exists("inputdata_spatial_discretization_planar_km", p)) {
            pAg$inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km
          } else {
            pAg$inputdata_spatial_discretization_planar_km = 1
          }
        }
        if (!exists("variabletomodel", pAg)) pAg$variabletomodel = "z"
      }

  } else if (source_data_class %in% "modelled_carstm" ) {
      # copy of param list for global analysis in aegis.bathymetry/inst/scripts/02.bathymetry.carstm.R
      pAg = aegis.bathymetry::bathymetry_carstm( DS = "parameters_production" )
  }


  # load input data or reformat it
   if (source_data_class=="rawdata") {

      B = bathymetry.db ( p=pAg, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
#      Bnames = c("lon", "lat", "z")
   } else if (source_data_class=="aggregated_rawdata") {

      B = bathymetry.db ( p=pAg, DS="aggregated_data" )
#       Bnames = c("z.mean", "z.sd",  "z.n", "plon", "plat", "lon", "lat")
      B$z = B$z.mean

   } else if (source_data_class=="modelled_stmv") {

      B = bathymetry.db(p=pAg, DS="baseline", varnames="all" )
    # Bnames = c( "plon", "plat", "z", "z.lb", "z.ub", "dZ", "ddZ", "b.sdTotal", "b.rsquared", "b.ndata", "b.sdSpatial", "b.sdObs", "b.phi", "b.nu", "b.localrange" )
      zname = "z"

   } else if (source_data_class=="modelled_carstm") {

      Bcarstm = carstm_summary( p=pAg ) # to load currently saved sppoly  ("locs")
      B = areal_units( p=pAg )
      bm = match( B$AUID, Bcarstm$AUID )
      B$z  = Bcarstm$z.predicted[ bm ]
      B$z.se = Bcarstm$z.predicted_se[ bm ]
      Bcarstm = NULL
      zname = "z"
  }

  require(raster)
  raster_template = raster(extent(locs)) # +1 to increase the area
  res(raster_template) = pAg$areal_units_resolution_km *1000 # in meters
  crs(raster_template) = projection(locs) # transfer the coordinate system to the raster

  if ( is.null(vnames) ) {
    vnames = zname
  } else {
    vnames = intersect( setdiff( names(B), c( "plon", "plat", "lon", "lat") ), vnames)
    if ( length(vnames) ==0 ) vnames = "z"
  }


  if ( output_data_class=="points" & source_data_class %in% c("rawdata", "aggregated_rawdata", "modelled_stmv" ) )  {

      B = lonlat2planar( B, proj.type=pAg$aegis_proj4string_planar_km )
      B$plon = round(B$plon / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
      B$plat = round(B$plat / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
      B_map = paste(B$plon, B$plat, sep=".")

      locs = lonlat2planar( locs, proj.type=pAg$aegis_proj4string_planar_km )
      locs$plon = round(locs$plon / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
      locs$plat = round(locs$plat / pAg$inputdata_spatial_discretization_planar_km + 1 ) * pAg$inputdata_spatial_discretization_planar_km
      locs_map = paste(locs$plon, locs$plat, sep=".")
      locs_index = match( locs_map, B_map )

      return( B[locs_index, vnames] )
  }

  if ( output_data_class=="points" & source_data_class=="modelled_carstm") {
      # convert to raster then match
      locs = as( locs, "sf" )
      for (vn in vnames) {
        Bf = fasterize::fasterize( as(B, "sf"), raster_template, field=vn )
        vn1 = paste(vn, "mean", sep="." )
        vn2 = paste(vn, "sd", sep="." )
        locs[, vn1] = raster::extract( Bf, locs, fun=mean, na.rm=TRUE)
        locs[, vn2] = raster::extract( Bf, locs, fun=sd, na.rm=TRUE)
      }
      return( as( locs, "Spatial" ) )
  }


  if ( output_data_class=="areal_units" & source_data_class %in% c("rawdata", "aggregated_rawdata", "modelled_stmv" ) ) {

      Bsf = sf::st_as_sf( B, coords=c("lon", "lat") )
      for (vn in vnames) {
        Bf = sp::over( as(B, "Spatial"), locs )
        vn1 = paste(vn, "mean", sep="." )
        locs[vn1] = raster::extract( Bf, locs, fun=mean, na.rm=TRUE)
        vn2 = paste(vn, "sd", sep="." )
        locs[vn2] = raster::extract( Bf, locs, fun=sd, na.rm=TRUE)
      }
      return(locs[, vnames] )
  }

  if (output_data_class=="areal_units" &  source_data_class=="modelled_carstm") {
      # convert to raster then match
      for (vn in vnames) {
        Bf = fasterize::fasterize( as(B, "sf"), raster_template, field=vn )
        Bf = sp::over( as(B, "Spatial"), locs )
        ii = which (is.finite( Bf ))

        vn1 = paste(vn, "mean", sep="." )
        locs[vn1] = raster::extract( Bf, locs, fun=mean, na.rm=TRUE)
        vn2 = paste(vn, "sd", sep="." )
        locs[vn2] = raster::extract( Bf, locs, fun=sd, na.rm=TRUE)
      }
      out = locs[, setdiff( names(locs), c( "plon", "plat", "lon", "lat") ) ]

      return(out)
  }

}


