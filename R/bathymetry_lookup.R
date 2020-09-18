
bathymetry_lookup = function( p, locs, vnames=NULL, output_data_class=NULL, source_data_class="aggregated_rawdata" ){

  # locs = output point locations or sp_polygon
  locs_class = class(locs)

  if (is.null(output_data_class)){
    if (any( grepl("sf", locs_class)) | any(grepl("SpatialPolygons", locs_class) )) {
      output_data_class = "areal_units"
    } else {
      output_data_class = "points"
    }
  }


  require(aegis.bathymetry)

  # set up parameters for input data
  if ( source_data_class %in% c("rawdata", "aggregated_rawdata", "modelled_stmv" ) ) {
      p_source = spatial_parameters( spatial_domain=p$spatial_domain )

      if (source_data_class=="modelled_stmv") {

        p_source = bathymetry_parameters(p=p, project_class="stmv")

      } else {
        # minimal info required ...
        if (!exists("inputdata_spatial_discretization_planar_km", p_source)) {
          if (!exists("inputdata_spatial_discretization_planar_km", p)) {
            p_source$inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km
          } else {
            p_source$inputdata_spatial_discretization_planar_km = 1
          }
        }
        if (!exists("variabletomodel", p_source)) p_source$variabletomodel = "z"
      }

  } else if (source_data_class %in% "modelled_carstm" ) {
      # copy of param list for global analysis in aegis.bathymetry/inst/scripts/02.bathymetry.carstm.R
      p_source = aegis.bathymetry::bathymetry_carstm( DS = "parameters_production" )
  }


  # load input data or reformat it
   if (source_data_class=="rawdata") {

      B = bathymetry.db ( p=p_source, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
#      Bnames = c("lon", "lat", "z")
   } else if (source_data_class=="aggregated_rawdata") {

      B = bathymetry.db ( p=p_source, DS="aggregated_data" )
#       Bnames = c("z.mean", "z.sd",  "z.n", "plon", "plat", "lon", "lat")
      B$z = B$z.mean
      B$z.mean = NULL

   } else if (source_data_class=="modelled_stmv") {

      B = bathymetry.db(p=p_source, DS="baseline", varnames="all" )
    # Bnames = c( "plon", "plat", "z", "z.lb", "z.ub", "dZ", "ddZ", "b.sdTotal", "b.rsquared", "b.ndata", "b.sdSpatial", "b.sdObs", "b.phi", "b.nu", "b.localrange" )
      zname = "z"

   } else if (source_data_class=="modelled_carstm") {

      Bcarstm = carstm_summary( p=p_source ) # to load currently saved sppoly  ("locs")
      B = areal_units( p=p_source )
      bm = match( B$AUID, Bcarstm$AUID )
      B$z  = Bcarstm$z.predicted[ bm ]
      B$z.se = Bcarstm$z.predicted_se[ bm ]
      Bcarstm = NULL
      zname = "z"
  }

  Bnames = setdiff( names(B), c("AUID", "uid", "layer", "plon", "plat", "lon", "lat", "au_sa_km2",
    "cfanorth_surfacearea", "cfasouth_surfacearea", "cfa23_surfacearea",  "cfa24_surfacearea", "cfa4x_surfacearea" ) )


  if ( output_data_class=="points" & source_data_class %in% c("rawdata", "aggregated_rawdata", "modelled_stmv" ) )  {
      B = lonlat2planar( B, proj.type=p_source$aegis_proj4string_planar_km )
      B$plon = round(B$plon / p_source$inputdata_spatial_discretization_planar_km + 1 ) * p_source$inputdata_spatial_discretization_planar_km
      B$plat = round(B$plat / p_source$inputdata_spatial_discretization_planar_km + 1 ) * p_source$inputdata_spatial_discretization_planar_km
      B_map = paste(B$plon, B$plat, sep=".")
      locs = lonlat2planar( locs, proj.type=p_source$aegis_proj4string_planar_km )
      locs$plon = round(locs$plon / p_source$inputdata_spatial_discretization_planar_km + 1 ) * p_source$inputdata_spatial_discretization_planar_km
      locs$plat = round(locs$plat / p_source$inputdata_spatial_discretization_planar_km + 1 ) * p_source$inputdata_spatial_discretization_planar_km
      locs_map = paste(locs$plon, locs$plat, sep=".")
      locs_index = match( locs_map, B_map )
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return( B[locs_index, vnames] )
  }


  if ( output_data_class=="points" & source_data_class=="modelled_carstm") {
      # convert to raster then match

      require(raster)
      raster_template = raster(extent(locs)) # +1 to increase the area
      res(raster_template) = p_source$areal_units_resolution_km  # crs usually in meters, but aegis's crs is in km
      crs(raster_template) = projection(locs) # transfer the coordinate system to the raster

      locs = sf::st_as_sf( as.data.frame(locs), coords=c(1, 2) )
      st_crs(locs) = crs(B)
      for (vn in Bnames) {
        Bf = fasterize::fasterize( as(B, "sf"), raster_template, field=vn )
        vn2 = paste(vn, "sd", sep="." )
        locs[, vn ] = raster::extract( Bf, locs, fun=mean, na.rm=TRUE)
        locs[, vn2] = raster::extract( Bf, locs, fun=sd, na.rm=TRUE)
      }
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return( as.matrix(locs[[vnames]]) )
  }


  if ( output_data_class=="areal_units" & source_data_class %in% c("rawdata", "aggregated_rawdata", "modelled_stmv" ) ) {
      Bsf = sf::st_as_sf( B, coords=c("lon", "lat") )
      st_crs(Bsf) = CRS( projection_proj4string("lonlat_wgs84") )
      Bsf = sf::st_transform( Bsf, crs=CRS(proj4string(locs)) )
      for (vn in Bnames) {
        vn2 = paste(vn, "sd", sep="." )
        slot(locs,"data")[,vn] = sp::over( locs, as(Bsf, "Spatial"), fn=mean, na.rm=TRUE )[,vn]
        slot(locs,"data")[,vn2] = sp::over( locs, as(Bsf, "Spatial"), fn=sd, na.rm=TRUE )[,vn]
      }
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return(locs[, vnames] )
  }


  if (output_data_class=="areal_units" &  source_data_class=="modelled_carstm") {
      # convert to raster then match
      require(raster)
      raster_template = raster(extent(locs)) # +1 to increase the area
      res(raster_template) = p_source$areal_units_resolution_km  # crs usually in meters, but aegis's crs is in km
      crs(raster_template) = projection(locs) # transfer the coordinate system to the raster

      for (vn in Bnames) {
        Bf = fasterize::fasterize( as(B, "sf"), raster_template, field=vn )
        vn2 = paste(vn, "sd", sep="." )
        slot(locs,"data")[,vn] = sp::over( locs, as(B, "SpatialPolygonsDataFrame"), fn=mean, na.rm=TRUE )[,vn]
        slot(locs,"data")[,vn2] = sp::over( locs, as(B, "SpatialPolygonsDataFrame"), fn=sd, na.rm=TRUE )[,vn]
      }
      vnames = intersect( names(B), vnames )
      if ( length(vnames) ==0 ) vnames=names(B) # no match returns all
      return(locs[,vnames])
  }

}


