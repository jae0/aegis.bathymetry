
  bathymetry_db = function( p=NULL, DS=NULL, varnames=NULL, redo=FALSE, sppoly=NULL, additional.data=c("sc_survey", "groundfish", "lobster", "sc_logbooks" ), ... ) {

    #\\ Note inverted convention: depths are positive valued
    #\\ i.e., negative valued for above sea level and positive valued for below sea level
    if ( is.null(p))  {
      p_add = list(...)
      if (length(p_add) > 0 ) {
        p = bathymetry_parameters(...)
      } else {
        p = bathymetry_parameters()
      }
    }


    if ( DS=="gebco") {
      #library(RNetCDF)
      # request at: https://www.bodc.ac.uk/data/online_delivery/gebco/ [ jae.choi@dfo ] ... / gate.gate
      # extent: (WSEN) = -+72,36,-45.,53
      # and saved as: bio.data/bathymetry/data/gebco.{xyz,nc}  # still waiting
      # and xz compressed
      fn = file.path( p$datadir, "bathymetry.gebco.rdq" )
      if (file.exists (fn) ) {
        return(read_write_fast(fn))
      }

      fn_local = file.path( p$datadir, "gebco.xyz.xz") # xz compressed file
      nc = open.nc(bathy_fname)
      read.nc(nc)
      array(tmp$z, dim=tmp$dim)
      gebco = read.table( xzfile( fn_local ) )
      names(gebco) = c("lon", "lat", "z")
      gebco$z = - gebco$z
      # levelplot( log(z+ min(gebco$z) ))~lon+lat, gebco, aspect="iso")
      read_write_fast( gebco, file=fn )
    }

    # --------------

    if ( DS=="etopo1") {
      # etopo1_bedrock.xyz ---> 1 min resolution
      # extent: (WSEN) = -72,36,-45.,53
      # download manually from:  http://maps.ngdc.noaa.gov/viewers/wcs-client/
      # and saved as: bio.data/bathymetry/data/etopo1_bedrock.xyz
      # and xz compressed
      fn = file.path( p$datadir, "bathymetry.etopo1.rdz" )
      if (file.exists (fn) ) {
        return(read_write_fast(fn))
      }
      fn_local = file.path( p$datadir, "etopo1_bedrock.xyz.xz") # xz compressed file
      etopo1 = read.table( xzfile( fn_local ) )
      names(etopo1) = c("lon", "lat", "z")
      etopo1$z = - etopo1$z
      # levelplot( log(z+ min(etopo1$z) ))~lon+lat, etopo1, aspect="iso")
      read_write_fast( etopo1, file=fn )
    }

    # --------------

    if ( DS =="Greenlaw_DEM") {
      # DEM created 2014
      # GCS_WGS_1984, UTM_Zone_20N; spheroid:: 6378137.0, 298.257223563
      # 322624071 "grid points
      # 50 m  horizontal resolution
      # depth range: -5053.6 to 71.48 m
      fn = file.path( p$datadir, "bathymetry.greenlaw.rdz" )
      if (file.exists (fn) ) {
        return(read_write_fast(fn))
      }

message("FIXE ME::: deprecated libs, use sf/stars")

#      require(rgdal)
      demfile.adf = file.path( p$datadir, "greenlaw_DEM", "mdem_50", "w001001.adf" )  # in ArcInfo adf format
      dem = new( "GDALReadOnlyDataset", demfile.adf )
      # gdem = asSGDF_GROD( dem, output.dim=dim(dem) ) # regrid to another dim
      # gdem = getRasterData(dem) # in matrix format
      gdem = getRasterTable(dem) # as a data frame
      names(gdem) = c("plon", "plat", "z")
      gdem = gdem[ is.finite( gdem$z ) , ]
      gdem$plon = gdem$plon / 1000
      gdem$plat = gdem$plat / 1000
      gdem = planar2lonlat( gdem, "utm20" )  # plon,plat in meters but crs for utm20 in km
      gdem = gdem[, c("lon", "lat", "z") ]
      read_write_fast( gdem, file=file.path( p$datadir, "bathymetry.greenlaw.rdz") )
    }




    # --------------

    if (  DS %in% c("z.lonlat.rawdata.redo", "z.lonlat.rawdata") ) {
			# raw data minimally modified all concatenated, dups removed
      fn = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.rdz" )

      if (DS =="z.lonlat.rawdata" ) {
        return( read_write_fast(fn) )
      }

			# this data was obtained from CHS via David Greenberg in 2004; range = -5467.020, 383.153; n=28,142,338
      fn_nwa = file.path( p$datadir, "nwa.chs15sec.xyz.xz") # xz compressed file
      chs15 = read.table( xzfile( fn_nwa ) )
      setDT(chs15)
      names(chs15) = c("lon", "lat", "z")
      chs15$z = - chs15$z
      chs15$source = "chs15sec"
 
      # pei = which( chs15$lon < -60.5 & chs15$lon > -64.5 & chs15$lat>45.5 & chs15$lat<48.5 )
      # levelplot( z~lon+lat, data=chs15[pei,] )

if (0) {
  
  # oversmooth?

      # Michelle Greenlaw's DEM from 2014
      # range -3000 to 71.5 m; n=155,241,029 .. but mostly interpolated
      gdem = bathymetry_db( DS="Greenlaw_DEM" )
      setDT(gdem)

      gdem$z = - gdem$z

      # pei = which( gdem$lon < -60.5 & gdem$lon > -65 & gdem$lat>45.5 & gdem$lat<49 )
      # levelplot( z~I(round(lon,3))+I(round(lat,3)), data=gdem[pei,] )

      # bad boundaries in Greenlaw's gdem:
      # southern Gulf of St lawrence has edge effects
      bd1 = rbind( c( -62, 46.5 ),
                   c( -61, 47.5 ),
                   c( -65, 47.5 ),
                   c( -65, 46.2 ),
                   c( -64, 46.2 ),
                   c( -62, 46.5 ) )

      a = which( point.in.polygon( gdem$lon, gdem$lat, bd1[,1], bd1[,2] ) != 0 )
      gdem = gdem[ -a,]

      # remove also the northern and eastern margins for edge effects
      gdem = gdem[ which( gdem$lat <  47.1) ,]
      gdem = gdem[ which( gdem$lon < -56.5) ,]
 
      # temporary break up of data to make it functional in smaller RAM systems
      gdem$source = "greenlaw50m"
    

}

      # chs and others above use chs depth convention: "-" is below sea level,
			# in snowcrab and groundfish convention "-" is above sea level
			# retain postive values at this stage to help contouring near coastlines
      

      etopo1 = bathymetry_db( DS="etopo1" )
      setDT(etopo1)
      etopo1$source = "etopo1min"

      bathy = rbind( chs15, etopo1 )
      rm(etopo1)
      rm(chs15) 
      gc()


			if ( "sc_survey" %in% additional.data ) {
        # range from 23.8 to 408 m below sea level ... these have dropped the "-" for below sea level; n=5925 (in 2014)
        # project.library( "bio.snowcrab")
        sc = bio.snowcrab::snowcrab.db( DS="set.clean")[, c("lon", "lat", "z") ]
				setDT(sc)

        sc = sc [ is.finite(lon) & is.finite(lat) &is.finite(z) ,]

				j = which(duplicated(sc))
        if (length (j) > 0 ) sc = sc[-j,]
        sc$source = "sc_survey"

        bathy = rbind( bathy, sc )
			  # p = p0
        rm (sc); gc()

        #sc$lon = round(sc$lon,1)
        #sc$lat = round(sc$lat,1)
        # contourplot( z~lon+lat, sc, cuts=10, labels=F )
      }


      # Too noisy:
			# if ( "sc_logbooks" %in% additional.data ) {
      #   # range from 23.8 to 408 m below sea level ... these have dropped the "-" for below sea level; n=5925 (in 2014)
      #   # project.library( "bio.snowcrab")
      #   scl = bio.snowcrab::logbook.db( DS="logbook")[, c("lon", "lat", "depth") ]
			# 	setDT(scl)
      #   setnames(scl, "depth", "z")

      #   scl = scl[ is.finite(lon) & is.finite(lat) &is.finite(z) ,]
      #   scl$depth[ scl$z > 10 & scl$z < 300 ]

			# 	j = which(duplicated(scl))
      #   if (length (j) > 0 ) scl = scl[-j,]
      #   scl$source = "sc_logbooks"

      #   bathy = rbind( bathy, scl )
			#   # p = p0
      #   rm (scl); gc()

      #   #scl$lon = round(scl$lon,1)
      #   #scl$lat = round(scl$lat,1)
      #   # contourplot( z~lon+lat, scl, cuts=10, labels=F )
      # }


      if ( "groundfish" %in% additional.data ) {
        # n=13031; range = 0 to 1054

        warning( "Should use bottom contact estimates as a priority ?" )
				gf = groundfish_survey_db( DS="set.base" )[, c("lon","lat", "sdepth") ]
				setDT(gf)
        names(gf) = c("lon", "lat", "z")

        gf = gf[ is.finite(lon) & is.finite(lat) &is.finite(z)  , ]
				j = which(duplicated(gf))
        if (length (j) > 0 ) gf = gf[-j,]
 				gf$source = "gf_survey"

        bathy = rbind( bathy, gf )
        rm (gf); gc()

        #gf$lon = round(gf$lon,1)
        #gf$lat = round(gf$lat,1)
        #contourplot( z~lon+lat, gf, cuts=10, labels=F )

			}

      if ( "lobster" %in% additional.data ) {
        current.year = lubridate::year(lubridate::now())
        p0temp = aegis.temperature::temperature_parameters( yrs=1900:current.year )
        lob = temperature_db( p=p0temp, DS="lobster", yr=1900:current.year ) # FSRS data ...  new additions have to be made at the rawdata level manually; yr must be passed to retrieve data ..
        lob = lob[, c("lon","lat", "z") ]
        setDT(lob)
        lob$source = "fsrs"

        lob = lob[ is.finite(lon) & is.finite(lat) &is.finite(z)   ,]
        j = which(duplicated(lob))
        if (length (j) > 0 ) lob = lob[-j,]
        bathy = rbind( bathy, lob )
        rm (lob); gc()

        #lob$lon = round(lob$lon,1)
        #lob$lat = round(lob$lat,1)
        #contourplot( z~lon+lat, lob, cuts=10, labels=F )

      }
  
      bid = paste( round(bathy$lon,4), round(bathy$lat,4), round(bathy$z, 1) ) # crude way to find dups
      bathy = bathy[!duplicated(bid),]
      
      read_write_fast( bathy, file=fn )
  
      return ( fn )
    }


    # ------------------------------

    if ( DS=="aggregated_data") {

      if (!exists("inputdata_spatial_discretization_planar_km", p) )  p$inputdata_spatial_discretization_planar_km = 1

      fn = file.path( p$datadir, paste( "bathymetry", "aggregated_data", p$spatial_domain,  round(p$inputdata_spatial_discretization_planar_km, 6) , "rdz", sep=".") )
      M = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          message("Using: ", fn)
          M=read_write_fast( fn)
          return( M )
        }
      }
      message("Making: ", fn)

      M = bathymetry_db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
   
      setDT(M)
      M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
      gc()

      M = lonlat2planar( M, proj.type=p$aegis_proj4string_planar_km)  # first ensure correct projection
      gc()

      M$lon = NULL
      M$lat = NULL
      setnames(M, p$variabletomodel, "z" )

      # thin data a bit ... remove potential duplicates and robustify

      M$plon = trunc(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      M$plat = trunc(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km

      gc()
       
      # filter_by_spatial_domain removes depths < 0 (above sea level)
      M = M[ filter_by_spatial_domain( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...   filter depths

      if (exists("quantile_bounds", p)) {
        TR = quantile(M$z, probs=p$quantile_bounds, na.rm=TRUE )
        keep = which( M$z >=  TR[1] & M$z <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
        keep = NULL
        gc()
      }

      M = M[, .( mean=mean(z, na.rm=TRUE), sd=sd(z, na.rm=TRUE), n=length(which(is.finite(z))) ), by=list(plon, plat) ]

      colnames(M) = c( "plon", "plat", paste( p$variabletomodel, c("mean", "sd", "n"), sep=".") )
      M = planar2lonlat( M, p$aegis_proj4string_planar_km )
 
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")
      setDF(M)
      read_write_fast(M, file=fn)

      return( M )
    }


    # ------------------------------

   
    if ( DS=="areal_units_input" ) {

      
      outdir = file.path( p$data_root, "modelled", p$carstm_model_label ) 
      fn = file.path( outdir, "areal_units_input.rdz"  )
      if ( !file.exists(outdir)) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )

      xydata = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          xydata = read_write_fast( fn)
          return( xydata )
        }
      }
      xydata = bathymetry_db( p=p, DS="aggregated_data"   )  #
      names(xydata)[which(names(xydata)=="z.mean" )] = "z"
      xydata = xydata[ , c("lon", "lat"  )]

      read_write_fast(xydata, file=fn )
      return( xydata )
    }

    # -----------------------


    if ( DS=="carstm_inputs") {

      # prediction surface
      #\\ Note inverted convention: depths are positive valued
      #\\ i.e., negative valued for above sea level and positive valued for below sea level
      crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
      if (is.null(sppoly)) sppoly = areal_units( p=p )  # will redo if not found
      sppoly = st_transform(sppoly, crs=crs_lonlat )
      areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

      fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs", areal_units_fn, sep="_") )
      if (p$carstm_inputs_prefilter =="rawdata") {
        fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs_rawdata", areal_units_fn, sep="_") )
      }
      if (p$carstm_inputs_prefilter =="sampled") {
        fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs_sampled", areal_units_fn, sep="_") )
      }

      # inputs are shared across various secneario using the same polys
      #.. store at the modeldir level as default
      outputdir = dirname( fn )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      if (!redo)  {
        if (file.exists(fn)) {
          M = read_write_fast( fn)
          return( M )
        }
      }

      # reduce size
      if (p$carstm_inputs_prefilter =="aggregated") {
        M = bathymetry_db ( p=p, DS="aggregated_data"   )  # 16 GB in RAM just to store!
        names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel

      } else if (p$carstm_inputs_prefilter =="sampled") {
        M = bathymetry_db ( p=p, DS="z.lonlat.rawdata"  )  # 16 GB in RAM just to store!

        require(data.table)
        setDT(M)
        names(M)[which(names(M)=="z") ] = p$variabletomodel

        M = M[ which( !duplicated(M)), ]
        M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")

    # thin data a bit ... remove potential duplicates and robustify
        M = lonlat2planar( M, proj.type=p$aegis_proj4string_planar_km )  # first ensure correct projection
        setDT(M)

        M$plon = trunc(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
        M$plat = trunc(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
    
        M = M[,.SD[sample(.N, min(.N, p$carstm_inputs_prefilter_n))], by =list(plon, plat) ]  # compact, might be slightly slower
        # M = M[ M[, sample(.N, min(.N, p$carstm_inputs_prefilter_n) ), by=list(plon, plat)], .SD[i.V1], on=list(plon, plat), by=.EACHI]  # faster .. just a bit
        setDF(M)

      } else {

        M = bathymetry_db ( p=p, DS="z.lonlat.rawdata"  )  # 16 GB in RAM just to store!
        names(M)[which(names(M)=="z") ] = p$variabletomodel
        M = M[ which( !duplicated(M)), ]
        M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
        # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")

      }
    

      if (p$carstm_inputs_prefilter != "aggregated") {
        if (exists("quantile_bounds", p)) {
          TR = quantile(M[[p$variabletomodel]], probs=p$quantile_bounds, na.rm=TRUE )
          keep = which( M[[p$variabletomodel]] >=  TR[1] & M[[p$variabletomodel]] <=  TR[2] )
          if (length(keep) > 0 ) M = M[ keep, ]
          keep = NULL
          gc()
        }
      }

        # if (exists("quantile_bounds", p)) {
        #   TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds, na.rm=TRUE )
        #   keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        #   if (length(keep) > 0 ) M = M[ keep, ]
        # }

      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")
          # p$quantile_bounds = c(0.0005, 0.9995)
 
      M = M[ which (M$z > 1), ]

      M$AUID = st_points_in_polygons(
        pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname="AUID"
      )
      M = M[ which(!is.na(M$AUID)),]
      M$AUID = as.character( M$AUID )  # match each datum to an area

      eps = .Machine$double.eps
      M$z = M$z + runif( nrow(M), min=-eps, max=eps )

      M$lon = NULL
      M$lat = NULL

      M$tag = "observations"

      APS = st_drop_geometry(sppoly)
      
      gc()

      APS[, p$variabletomodel] = NA
      APS$AUID = as.character( APS$AUID )
      APS$tag ="predictions"

      vn = c("z", "tag", "AUID" )

      M = rbind( M[, vn], APS[, vn] )
      APS = NULL

      #required for carstm formulae
      M$space = match( M$AUID, as.character(sppoly$AUID) ) ## must be index matching nb graphs
      
      read_write_fast( M, file=fn )
      return( M )
    }


    # ------------------------------


    if ( DS %in% c("bathymetry", "stmv_inputs", "stmv_inputs_redo" )) {

      fn = file.path( p$modeldir, paste( "bathymetry", "stmv_inputs", "rdz", sep=".") )
      if (DS %in% c("bathymetry", "stmv_inputs") ) {
        hm = read_write_fast( fn)
        return( hm )
      }

      B = bathymetry_db ( p=p, DS="aggregated_data", redo=TRUE )  # 16 GB in RAM just to store!
      B$lon = NULL
      B$lat = NULL
      names(B)[which(names(B) == paste(p$variabletomodel, "mean", sep="."))] = p$variabletomodel

      print( "Warning: this needs a lot of RAM .. ~60GB depending upon resolution of discretization .. a few hours " )

      hm = list( input=B, output=list( LOCS = spatial_grid(p) ) )
      B = NULL; gc()
      read_write_fast( hm, file=fn )
      hm = NULL
      gc()
      return(fn)

    }

    # ------------------------------


    if ( DS %in% c("stmv_inputs_highres", "stmv_inputs_highres_redo" )) {

      fn = file.path( p$modeldir, paste( "bathymetry", "stmv_inputs_highres", "rdz", sep=".") )
      if (DS %in% c("stmv_inputs_highres") ) {
        hm = read_write_fast( fn)
        return( hm )
      }

      B = bathymetry_db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!

      # p$quantile_bounds = c(0.0005, 0.9995)
      if (exists("quantile_bounds", p)) {
        TR = quantile(B[[p$variabletomodel]], probs=p$quantile_bounds, na.rm=TRUE )
        keep = which( B[[p$variabletomodel]] >=  TR[1] & B[[p$variabletomodel]] <=  TR[2] )
        if (length(keep) > 0 ) B = B[ keep, ]
      }

      # thin data a bit ... remove potential duplicates and robustify
      B = lonlat2planar( B, proj.type=p$aegis_proj4string_planar_km )  # first ensure correct projection

      B$lon = NULL
      B$lat = NULL

      attr( B, "proj4string_planar" ) =  p$aegis_proj4string_planar_km

      hm = list( input=B, output=list( LOCS = spatial_grid(p) ) )
      B = NULL; gc()
      read_write_fast( hm, file=fn )
      hm = NULL
      gc()
      return(fn)

    }

    # ----------------

    if ( DS == "landmasks.create" ) {

      ## NOTE::: This does not get used

      # on resolution of predictions
      V = SpatialPoints( planar2lonlat( spatial_grid(p), proj.type=p$aegis_proj4string_planar_km )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires", regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="predictions", crs=p$aegis_proj4string_planar_km)

      # on resolution of statistics
      sbox = list(
        plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
        plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
      V = as.matrix( expand.grid( sbox$plons, sbox$plats ))
      names(V) = c("plon", "plat")
      V = SpatialPoints( planar2lonlat( V, proj.type=p$aegis_proj4string_planar_km )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires",regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="statistics", crs=p$aegis_proj4string_planar_km)
    }

    #-------------------------

    if ( DS %in% c("complete", "complete.redo" )) {
      #// merge all stmv results and compute stats and warp to different grids
      outdir = file.path( p$modeldir, p$stmv_model_label, p$project_class, paste(  p$stmv_global_modelengine, p$stmv_local_modelengine, sep="_") )

      fn = file.path( outdir, paste( "bathymetry", "complete", p$spatial_domain, "rdz", sep=".") )

      if ( DS %in% c( "complete") ) {
        Z = NULL
        if ( file.exists ( fn) ) Z=read_write_fast( fn)
        return( Z )
      }

      nr = p$nplons
      nc = p$nplats

      # data prediction grid
      B = spatial_grid( p)
      Bmean = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
      Blb = stmv_db( p=p, DS="stmv.prediction", ret="lb" )
      Bub = stmv_db( p=p, DS="stmv.prediction", ret="ub" )
      Z = data.frame( cbind(B, Bmean, Blb, Bub) )
      names(Z) = c( "plon", "plat", "z", "z.lb", "z.ub") # really Z.mean but for historical compatibility "z"
      B = Bmean = Blb = Bub = NULL

      # # remove land
      # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" )
      # Z$z[oc] = NA
      # Z$z.sd[oc] = NA

      Zmn = matrix( Z[,p$variabletomodel], nrow=nr, ncol=nc )  # means

      # first order central differences but the central term drops out:
      # diffr = ( ( Zmn[ 1:(nr-2), ] - Zmn[ 2:(nr-1), ] ) + ( Zmn[ 2:(nr-1), ] - Zmn[ 3:nr, ] ) ) / 2
      # diffc = ( ( Zmn[ ,1:(nc-2) ] - Zmn[ ,2:(nc-1) ] ) + ( Zmn[ ,2:(nc-1) ] - Zmn[ ,3:nc ] ) ) / 2
      diffr =  Zmn[ 1:(nr-2), ] - Zmn[ 3:nr, ]
      diffc =  Zmn[ ,1:(nc-2) ] - Zmn[ ,3:nc ]
      rm (Zmn); gc()

      dZ = ( diffr[ ,2:(nc-1) ] + diffc[ 2:(nr-1), ] ) / 2
      dZ = rbind( dZ[1,], dZ, dZ[nrow(dZ)] )  # top and last rows are copies .. dummy value to keep dim correct
      dZ = cbind( dZ[,1], dZ, dZ[,ncol(dZ)] )

      Z$dZ =  abs(c(dZ))

      # gradients
      ddiffr =  dZ[ 1:(nr-2), ] - dZ[ 3:nr, ]
      ddiffc =  dZ[ ,1:(nc-2) ] - dZ[ ,3:nc ]
      dZ = diffc = diffr = NULL

      ddZ = ( ddiffr[ ,2:(nc-1) ] + ddiffc[ 2:(nr-1), ] ) / 2
      ddZ = rbind( ddZ[1,], ddZ, ddZ[nrow(ddZ)] )  # top and last rows are copies .. dummy value to keep dim correct
      ddZ = cbind( ddZ[,1], ddZ, ddZ[,ncol(ddZ)] )
      Z$ddZ = abs(c(ddZ))

      # merge into statistics
      BS = stmv_db( p=p, DS="stmv.stats" )
      colnames(BS) = paste("b", colnames(BS), sep=".")
      Z = cbind( Z, BS )

      read_write_fast( Z, file=fn)

      BS = ddZ = ddiffc = ddiffr = NULL
      gc()

      # now warp to the other grids
      p0 = p  # the originating parameters

      Z0 = Z  # rename as 'Z' will be overwritten below
      L0 = Z0[, c("plon", "plat")]
      L0i = array_map( "xy->2", L0, gridparams=p0$gridparams )

      voi_bathy = setdiff( names(Z0), c("plon", "plat", "lon", "lat") )
      grids = setdiff( unique( p0$spatial_domain_subareas ), p0$spatial_domain )

      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( spatial_domain=gr ) #target projection
        # warping
          L1 = spatial_grid( p=p1 )
          L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
          L1 = planar2lonlat( L1, proj.type=p1$aegis_proj4string_planar_km )
          Z = L1
          L1$plon_1 = L1$plon # store original coords
          L1$plat_1 = L1$plat
          L1 = lonlat2planar( L1, proj.type=p0$aegis_proj4string_planar_km )

          p1$wght = fields::setup.image.smooth(
            nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
            theta=p1$pres/3, xwidth=4*p1$pres, ywidth=4*p1$pres )
        # theta=p1$pres/3 assume at pres most of variance is accounted ... correct if dense pre-intepolated matrices .. if not can be noisy

          for (vn in voi_bathy) {
            Z[,vn] = spatial_warp( Z0[,vn], L0, L1, p0, p1, "fast", L0i, L1i )
          }
        Z = Z[ , names(Z0) ]

        fn = file.path( outdir, paste( "bathymetry", "complete", p1$spatial_domain, "rdz", sep=".") )
        read_write_fast (Z, file=fn)
      }

      return(fn)

      if (0) {
        aoi = which( Z$z > 10 & Z$z < 500 )
        datarange = log( quantile( Z[aoi,"z"], probs=c(0.001, 0.999), na.rm=TRUE ))
        dr = seq( datarange[1], datarange[2], length.out=100)

        levelplot( log(z) ~ plon + plat, Z[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), at=dr, col.regions=rev(color.code( "seis", dr)) )
        levelplot( log(phi) ~ plon + plat, Z[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(range) ~ plon + plat, Z[aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )

      }
    }


    # ------------

    if (DS %in% c("baseline", "baseline.redo") ) {
      # form prediction surface in planar coords over the ocean

      outdir = file.path( p$modeldir, p$stmv_model_label, p$project_class, paste(  p$stmv_global_modelengine, p$stmv_local_modelengine, sep="_") )

      if ( DS=="baseline" ) {
        #  used to obtain coordinates
        fn = paste( "bathymetry", "baseline", p$spatial_domain, "rdz" , sep=".")
        outfile =  file.path( outdir, fn )
        Z = read_write_fast( outfile )
        Znames = names(Z)
        if (is.null(varnames)) varnames =c("plon", "plat")  # default is to send locs only .. different reative to all other data streams
        varnames = intersect( Znames, varnames )  # send anything that results in no match causes everything to be sent
        if (length(varnames) == 0) varnames=Znames  # no match .. send all
        Z = Z[ , varnames]
        return (Z)
      }

      for (domain in unique( c(p$spatial_domain_subareas, p$spatial_domain ) ) ) {
        pn = spatial_parameters( p=p, spatial_domain=domain )
        # if ( pn$spatial_domain == "snowcrab" ) {
        #   # NOTE::: snowcrab baseline == SSE baseline, except it is a subset so begin with the SSE conditions
        #   pn = spatial_parameters( p=pn, spatial_domain="SSE" )
        # }
        Z = bathymetry_db( p=pn, DS="complete" )
        Z = Z[ filter_by_spatial_domain( spatial_domain=domain, Z=Z ), ]

        # range checks
        ii = which( Z$dZ < exp(-6))
        if (length(ii) > 0) Z$dZ[ii] = exp(-6)

        ii = which( Z$dZ > 50 )
        if (length(ii) > 0) Z$dZ[ii] = 50

        ii = which( Z$ddZ < exp(-6))
        if (length(ii) > 0) Z$ddZ[ii] = exp(-6)

        ii = which( Z$ddZ > 20 )
        if (length(ii) > 0) Z$ddZ[ii] = 20

        fn = paste( "bathymetry", "baseline", domain, "rdz" , sep="." )
        outfile =  file.path( outdir, fn )

        read_write_fast (Z, file=outfile )
        print( outfile )
      }
      # require (lattice); levelplot( z~plon+plat, data=Z, aspect="iso")
      return( "completed" )
    }

    # ------------


    if (DS %in% c("baseline_prediction_locations", "baseline_prediction_locations.redo") ) {
      # form prediction surface in planar coords over the ocean

      outdir = p$modeldir

      if ( DS=="baseline_prediction_locations" ) {
        #  used to obtain coordinates
        fn = paste( "bathymetry", "baseline_prediction_locations", p$spatial_domain, "rdz" , sep=".")
        outfile =  file.path( outdir, fn )
        Z = NULL
        if (file.exists(fn)) {
          Z = read_write_fast( outfile )
          Znames = names(Z)
          if (is.null(varnames)) varnames =c("plon", "plat")  # default is to send locs only .. different reative to all other data streams
          varnames = intersect( Znames, varnames )  # send anything that results in no match causes everything to be sent
          if (length(varnames) == 0) varnames=Znames  # no match .. send all
          Z = Z[ , varnames]
          return (Z)
        }
      }

      pn = spatial_parameters( p=p, spatial_domain=p$spatial_domain )
      Z = spatial_grid(p)
      Z = planar2lonlat( Z,  proj.type=p$aegis_proj4string_planar_km   )
      
      pL = aegis.bathymetry::bathymetry_parameters( project_class="stmv"  )
      LUT= aegis_survey_lookuptable( aegis_project="bathymetry", 
        project_class="core", DS="aggregated_data", pL=pL )

      Z$z = aegis_lookup( pL=pL, LUT=LUT, 
        LOCS=Z[, c("lon", "lat")], project_class="core", 
        output_format="points", space_resolution=pn$pres, variable.name="z.mean" )  # core == unmodelled
      
      Z = Z[ filter_by_spatial_domain( spatial_domain=p$spatial_domain, Z=Z ), ]

      fn = paste( "bathymetry", "baseline_prediction_locations", p$spatial_domain, "rdz" , sep=".")
      outfile =  file.path( outdir, fn )

      read_write_fast (Z, file=outfile )
      print( outfile )

      # require (lattice); levelplot( z~plon+plat, data=Z, aspect="iso")
      return( Z )
    }

    # ------------



  }  # end bathymetry_db
