
  bathymetry.db = function( p=NULL, DS=NULL, varnames=NULL, redo=FALSE, ... ) {

    #\\ Note inverted convention: depths are positive valued
    #\\ i.e., negative valued for above sea level and positive valued for below sea level
    if ( is.null(p)) p = bathymetry_parameters(...)

    if ( !exists("project.name", p)) p$project.name = "bathymetry"
    if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project.name )
    if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
    if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )


    if ( DS=="gebco") {
      #library(RNetCDF)
      # request at: https://www.bodc.ac.uk/data/online_delivery/gebco/ [ jae.choi@dfo ] ... / gate.gate
      # extent: (WSEN) = -+72,36,-45.,53
      # and saved as: bio.data/bathymetry/data/gebco.{xyz,nc}  # still waiting
      # and xz compressed
      fn = file.path( p$datadir, "bathymetry.gebco.rdata" )
      if (file.exists (fn) ) {
        load(fn)
        return(gebco)
      }

      fn_local = file.path( p$datadir, "gebco.xyz.xz") # xz compressed file
      nc = open.nc(bathy_fname)
      read.nc(nc)
      array(tmp$z, dim=tmp$dim)
      gebco = read.table( xzfile( fn_local ) )
      names(gebco) = c("lon", "lat", "z")
      gebco$z = - gebco$z
      # levelplot( log(z+ min(gebco$z) ))~lon+lat, gebco, aspect="iso")
      save( gebco, file=fn, compress=TRUE )
    }

    # --------------

    if ( DS=="etopo1") {
      # etopo1_bedrock.xyz ---> 1 min resolution
      # extent: (WSEN) = -72,36,-45.,53
      # download manually from:  http://maps.ngdc.noaa.gov/viewers/wcs-client/
      # and saved as: bio.data/bathymetry/data/etopo1_bedrock.xyz
      # and xz compressed
      fn = file.path( p$datadir, "bathymetry.etopo1.rdata" )
      if (file.exists (fn) ) {
        load(fn)
        return(etopo1)
      }
      fn_local = file.path( p$datadir, "etopo1_bedrock.xyz.xz") # xz compressed file
      etopo1 = read.table( xzfile( fn_local ) )
      names(etopo1) = c("lon", "lat", "z")
      etopo1$z = - etopo1$z
      # levelplot( log(z+ min(etopo1$z) ))~lon+lat, etopo1, aspect="iso")
      save( etopo1, file=fn, compress=TRUE )
    }

    # --------------

    if ( DS =="Greenlaw_DEM") {
      # DEM created 2014
      # GCS_WGS_1984, UTM_Zone_20N; spheroid:: 6378137.0, 298.257223563
      # 322624071 "grid points
      # 50 m  horizontal resolution
      # depth range: -5053.6 to 71.48 m
      fn = file.path( p$datadir, "bathymetry.greenlaw.rdata" )
      if (file.exists (fn) ) {
        load(fn)
        return(gdem)
      }

      require(rgdal)
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
      save( gdem, file=file.path( p$datadir, "bathymetry.greenlaw.rdata"), compress=TRUE )
    }




    # --------------

    if (  DS %in% c("z.lonlat.rawdata.redo", "z.lonlat.rawdata") ) {
			# raw data minimally modified all concatenated, dups removed
      fn = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.rdata" )

      if (DS =="z.lonlat.rawdata" ) {
        load( fn )
        return( bathy )
      }

      print( "This is going to take a lot of RAM!")

			# this data was obtained from CHS via David Greenberg in 2004; range = -5467.020, 383.153; n=28,142,338
      fn_nwa = file.path( p$datadir, "nwa.chs15sec.xyz.xz") # xz compressed file
      chs15 = read.table( xzfile( fn_nwa ) )
      names(chs15) = c("lon", "lat", "z")
      # chs15 = chs15[ which( chs15$z < 1000 ) , ]
      chs15$z = - chs15$z

      # temporary break up of data to make it functional in smaller RAM systems
      chs1000_5000 = chs15[ which( chs15$z > 1000 ), ]
      u =  which(duplicated( chs1000_5000))
      if (length(u)>0) chs1000_5000 = chs1000_5000[-u,]

      chs0_1000 = chs15[ which( chs15$z <= 1000 ), ]
      u =  which(duplicated( chs0_1000 ))
      if (length(u)>0) chs0_1000 = chs0_1000[-u,]

      rm ( chs15); gc()

      fn0 = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_0_1000.rdata" )
      fn1 = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_1000_5000.rdata" )

      save ( chs0_1000, file=fn0 )
      save ( chs1000_5000, file=fn1 )

      rm ( chs0_1000, chs1000_5000 )
      gc()

      # pei = which( chs15$lon < -60.5 & chs15$lon > -64.5 & chs15$lat>45.5 & chs15$lat<48.5 )
      # levelplot( z~lon+lat, data=chs15[pei,] )


      # Michelle Greenlaw's DEM from 2014
      # range -3000 to 71.5 m; n=155,241,029 .. but mostly interpolated
      gdem = bathymetry.db( DS="Greenlaw_DEM" )
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
      gdem = gdem[- a,]

      # remove also the northern and eastern margins for edge effects
      gdem = gdem[ which( gdem$lat <  47.1) ,]
      gdem = gdem[ which( gdem$lon < -56.5) ,]

      fn0g = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_0_1000_gdem.rdata" )
      fn1g = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_1000_5000_gdem.rdata" )

      # temporary break up of data to make it functional in smaller RAM systems
      gdem1000_5000 = gdem[ which( gdem$z > 1000 ), ]
      # u =  which(duplicated( gdem1000_5000))
      # if (length(u)>0) gdem1000_5000 = gdem1000_5000[-u,]
      save ( gdem1000_5000, file=fn1g )
      rm( gdem1000_5000 ); gc()

      gdem0_1000 = gdem[ which( gdem$z <= 1000 ), ]
      # u =  which(duplicated( gdem0_1000 ))
      # if (length(u)>0) gdem0_1000 = gdem0_1000[-u,]
      save ( gdem0_1000, file=fn0g )
      rm( gdem0_1000 )
      rm( gdem) ;gc()
      gc()


      # chs and others above use chs depth convention: "-" is below sea level,
			# in snowcrab and groundfish convention "-" is above sea level
			# retain postive values at this stage to help contouring near coastlines

      bathy = bathymetry.db( DS="etopo1" )

      additional.data=c("snowcrab", "groundfish", "lobster")

			if ( "snowcrab" %in% additional.data ) {
        # range from 23.8 to 408 m below sea level ... these have dropped the "-" for below sea level; n=5925 (in 2014)
        # project.library( "bio.snowcrab")
        sc = bio.snowcrab::snowcrab.db( DS="set.clean")[,c("lon", "lat", "z") ]
				sc = sc [ which (is.finite( rowSums( sc ) ) ) ,]
				j = which(duplicated(sc))
        if (length (j) > 0 ) sc = sc[-j,]
        bathy = rbind( bathy, sc )
			  # p = p0
        rm (sc); gc()

        #sc$lon = round(sc$lon,1)
        #sc$lat = round(sc$lat,1)
        # contourplot( z~lon+lat, sc, cuts=10, labels=F )
      }

      if ( "groundfish" %in% additional.data ) {
        # n=13031; range = 0 to 1054

        warning( "Should use bottom contact estimates as a priority ?" )
				gf = groundfish.db( DS="set.base" )[, c("lon","lat", "sdepth") ]
				gf = gf[ which( is.finite(rowSums(gf) ) ) ,]
        names(gf) = c("lon", "lat", "z")
				j = which(duplicated(gf))
        if (length (j) > 0 ) gf = gf[-j,]
 				bathy = rbind( bathy, gf )
        rm (gf); gc()

        #gf$lon = round(gf$lon,1)
        #gf$lat = round(gf$lat,1)
        #contourplot( z~lon+lat, gf, cuts=10, labels=F )

			}

      if ( "lobster" %in% additional.data ) {
        current.year = lubridate::year(lubridate::now())
        p0temp = aegis.temperature::temperature_parameters( yrs=1900:current.year )
        lob = temperature.db( p=p0temp, DS="lobster", yr=1900:current.year ) # FSRS data ...  new additions have to be made at the rawdata level manually; yr must be passed to retrieve data ..
        lob = lob[, c("lon","lat", "z") ]
        lob = lob[ which( is.finite(rowSums(lob) ) ) ,]
        j = which(duplicated(lob))
        if (length (j) > 0 ) lob = lob[-j,]
        bathy = rbind( bathy, lob )
        rm (lob); gc()

        #lob$lon = round(lob$lon,1)
        #lob$lat = round(lob$lat,1)
        #contourplot( z~lon+lat, lob, cuts=10, labels=F )

      }


      u =  which(duplicated( bathy ))
      if (length(u)>0) bathy = bathy[ -u, ]
      rm (u)

      bathy0 = bathy[ which(bathy$z <= 1000), ]
      bathy1 = bathy[ which(bathy$z  > 1000), ]
      rm(bathy)

      gc()

      load( fn0)
      bathy0 = rbind( bathy0, chs0_1000 )
      rm(chs0_1000) ;gc()

      load( fn0g )
      bathy0 = rbind( bathy0, gdem0_1000 )
      rm ( gdem0_1000 ) ;gc()

      u =  which(duplicated( bathy0 ))
      if (length(u)>0) bathy0 = bathy0[ -u, ]

      fn0b = file.path( p$datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_0_1000_bathy.rdata" )
      save ( bathy0, file=fn0b )
      rm (bathy0); gc()

    # ---

      load( fn1 )
      bathy1 = rbind( bathy1, chs1000_5000 )
      rm (  chs1000_5000 ) ; gc()

      load( fn1g )
      bathy1 = rbind( bathy1, gdem1000_5000 )
      rm ( gdem1000_5000 ) ; gc()

      u =  which(duplicated( bathy1 ))
      if (length(u)>0) bathy1 = bathy1[ -u, ]
      rm (u)

      load( fn0b )

      bathy = rbind( bathy0, bathy1 )
      rm( bathy1, bathy0 ) ; gc()

      bid = paste( bathy$lon, bathy$lat, round(bathy$z, 1) ) # crude way to find dups
      bathy = bathy[!duplicated(bid),]
      save( bathy, file=fn, compress=T )

      # save ascii in case someone needs it ...
      fn.bathymetry.xyz = file.path( p$datadir, "bathymetry.canada.east.xyz" )
      fn.xz = xzfile( paste( fn.bathymetry.xyz, ".xz", sep="" ) )
      write.table( bathy, file=fn.xz, col.names=F, quote=F, row.names=F)
      system( paste( "xz",  fn.bathymetry.xyz ))  # compress for space

      file.remove(fn.bathymetry.xyz)
      file.remove( fn0 )
      file.remove( fn1 )
      file.remove( fn0g )
      file.remove( fn1g )
      file.remove( fn0b )

      return ( fn )
    }

    # ------------------------------

    if ( DS %in% c("bathymetry", "stmv.inputs", "stmv.inputs.redo" )) {

      fn = file.path( p$modeldir, paste( "bathymetry", "stmv.inputs", "rdata", sep=".") )
      if (DS %in% c("bathymetry", "stmv.inputs") ) {
        print( "Warning: stmv.inputs is loading from a saved instance ... use stmv.inputs.redo if data need refresh" )
        load( fn)
        return( hm )
      }

      print( "Warning: this needs a lot of RAM .. ~60GB depending upon resolution of discretization .. a few hours " )

      B = bathymetry.db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
      # B = B[ which(B$z > -100),]  # take part of the land to define coastline
      B = lonlat2planar( B, proj.type=p$internal.crs )
      B$lon = NULL
      B$lat = NULL
      gc()

      # thin data a bit ... remove potential duplicates and robustify
      B$plon = floor(B$plon / p$pres_discretization_bathymetry) * p$pres_discretization_bathymetry
      B$plat = floor(B$plat / p$pres_discretization_bathymetry) * p$pres_discretization_bathymetry
      zrez = 1  # meter
      bid = paste( B$plon, B$plat  )
      B = B[!duplicated(paste( bid, floor(B$z / zrez) * zrez )),]
      gc()
      B = B[ order(bid),]
      rm(bid); gc()

      # reduce size of B in steps ... this is too large to run in one go

      nr = nrow(B)
      ns = 20 # number of segments
      nl = floor( nr / ns )

      pos = 0 # index start position
      tmpfn = list()

      for ( nn in 1:ns ) {
        print( nn)
        ii = pos + (1:nl)
        out = B[ii,]
        tmpfn[[nn]] =  file.path( p$data_root, paste("tmp", nn, "rdata", sep=".") )
        save(out, file=tmpfn[[nn]])
        out = NULL; gc()
        pos = ii[nl] # update last index
      }

      if (pos < nr) {
        # get stragglers
        ii = (pos + 1) : nr
        out = B[ii,]
        tmpfn[[ns+1]] =  file.path( p$data_root, paste("tmp", (ns+1), "rdata", sep=".") )
        save(out, file=tmpfn[[ns+1]])
        out = NULL; gc()
      }
      B = NULL;
      gc()

      for ( nn in 1:ns ) {
        print(nn)
        load( tmpfn[[nn]] )
        bb = tapply( X=out$z, INDEX=list(paste( out$plon, out$plat) ), FUN = function(w) {median(w, na.rm=TRUE)}, simplify=TRUE )
        locs = sapply( strsplit(rownames(bb), " ", fixed=TRUE ), unlist)
        rownames(bb) = NULL
        B = rbind( B, cbind(bb, as.numeric(locs[1,]), as.numeric(locs[2,]) ) )
        bb = NULL
        locs = NULL
        out = NULL
        gc()
      }
      if (pos < nr) {
        # get stragglers
        load( tmpfn[[ns+1]] )
        bb = tapply( X=out$z, INDEX=as.list(out[, c("plon", "plat")]), FUN = function(w) {median(w, na.rm=TRUE)}, simplify=TRUE )
        locs = sapply( strsplit(rownames(bb), " ", fixed=TRUE ), unlist)
        rownames(bb) = NULL
        B = rbind( B, cbind(bb, as.numeric(locs[1,]), as.numeric(locs[2,]) ) )
        bb = NULL
        locs = NULL
        out = NULL
        gc()
      }
      bb = NULL
      gc()
      B = as.data.frame( B )
      names(B) = c("z", "plon", "plat")
      B = B[ which( is.finite( B$z )) ,]

      hm = list( input=B, output=list( LOCS = spatial_grid(p) ) )
      B = NULL; gc()
      save( hm, file=fn, compress=TRUE)

      for ( fnt in tmpfn ) {
        if (file.exists (fnt)) file.remove(fnt)
      }


    }

    # ----------------

    if ( DS == "landmasks.create" ) {

      ## NOTE::: This does not get used

      # on resolution of predictions
      V = SpatialPoints( planar2lonlat( spatial_grid(p), proj.type=p$internal.crs )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires", regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="predictions" ,internal.crs=p$internal.crs)

      # on resolution of statistics
      sbox = list(
        plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
        plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
      V = as.matrix( expand.grid( sbox$plons, sbox$plats ))
      names(V) = c("plon", "plat")
      V = SpatialPoints( planar2lonlat( V, proj.type=p$internal.crs )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires",regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="statistics" ,internal.crs=p$internal.crs)
    }

    #-------------------------

    if ( DS %in% c("complete", "complete.redo" )) {
      #// merge all stmv results and compute stats and warp to different grids

      fn = file.path( p$modeldir, paste( "bathymetry", "complete", p$spatial.domain, "rdata", sep=".") )

      if ( DS %in% c( "complete") ) {
        Z = NULL
        if ( file.exists ( fn) ) load( fn)
        return( Z )
      }

      nr = p$nplons
      nc = p$nplats

      # data prediction grid
      B = spatial_grid( p)
      Bmean = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
      Blb = stmv_db( p=p, DS="stmv.prediction", ret="lb" )
      Bub = stmv_db( p=p, DS="stmv.prediction", ret="ub" )
      Z = cbind(B, Bmean, Blb, Bub)
      names(Z) = c( "plon", "plat", "z", "z.lb", "z.ub") # really Z.mean but for historical compatibility "z"
      B = Bmean = Blb = Bub = NULL

      # # remove land
      # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" )
      # Z$z[oc] = NA
      # Z$z.sd[oc] = NA

      Zmn = matrix( Z$z, nrow=nr, ncol=nc )  # means

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

      save( Z, file=fn, compress=TRUE)

      BS = ddZ = ddiffc = ddiffr = NULL
      gc()

      # now warp to the other grids
      p0 = p  # the originating parameters

      Z0 = Z  # rename as 'Z' will be overwritten below
      L0 = Z0[, c("plon", "plat")]
      L0i = array_map( "xy->2", L0, gridparams=p0$gridparams )

      voi_bathy = setdiff( names(Z0), c("plon", "plat", "lon", "lat") )
      grids = setdiff( unique( p0$spatial.domain.subareas ), p0$spatial.domain )

      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( spatial.domain=gr ) #target projection
          L1 = spatial_grid( p1 )
          L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
          L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
          Z = L1
          L1$plon_1 = L1$plon # store original coords
          L1$plat_1 = L1$plat
          L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
          p1$wght = fields::setup.image.smooth(
            nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
            theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
          for (vn in voi_bathy) {
            Z[,vn] = spatial_warp( Z0[,vn], L0, L1, p0, p1, "fast", L0i, L1i )
          }
        Z = Z[ , names(Z0) ]
        fn = file.path( p$modeldir, paste( "bathymetry", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (Z, file=fn, compress=TRUE)
      }

      return(fn)

      if (0) {
        aoi = which( B$z > 10 & B$z < 500 )
        datarange = log( quantile( B[aoi,"z"], probs=c(0.001, 0.999), na.rm=TRUE ))
        dr = seq( datarange[1], datarange[2], length.out=100)

        levelplot( log(z) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), at=dr, col.regions=rev(color.code( "seis", dr)) )
        levelplot( log(phi) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(range) ~ plon + plat, B[aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )

      }
    }


    # ------------

    if (DS %in% c("baseline", "baseline.redo") ) {
      # form prediction surface in planar coords over the ocean

      if ( DS=="baseline" ) {

#        print(p$spatial.domain)
        outfile =  file.path( p$modeldir, paste( "bathymetry", "baseline", p$spatial.domain, "rdata" , sep=".") )
        Z = NULL
        load( outfile )
        Znames = names(Z)
        if (is.null(varnames)) varnames =c("plon", "plat")  # default is to send locs only .. different reative to all other data streams
        varnames = intersect( Znames, varnames )  # send anything that results in no match causes everything to be sent
        if (length(varnames) == 0) varnames=Znames  # no match .. send all
        Z = Z[ , varnames]
        return (Z)
      }

      for (domain in unique( c(p$spatial.domain.subareas, p$spatial.domain ) ) ) {
        pn = spatial_parameters( spatial.domain=domain )
        # if ( pn$spatial.domain == "snowcrab" ) {
        #   # NOTE::: snowcrab baseline == SSE baseline, except it is a subset so begin with the SSE conditions
        #   pn = spatial_parameters( p=pn, spatial.domain="SSE" )
        # }
        Z = bathymetry.db( p=pn, DS="complete" )
        Z = geo_subset( spatial.domain=domain, Z=Z )

        # range checks
        ii = which( Z$dZ < exp(-6))
        if (length(ii) > 0) Z$dZ[ii] = exp(-6)

        ii = which( Z$dZ > 50 )
        if (length(ii) > 0) Z$dZ[ii] = 50

        ii = which( Z$ddZ < exp(-6))
        if (length(ii) > 0) Z$ddZ[ii] = exp(-6)

        ii = which( Z$ddZ > 20 )
        if (length(ii) > 0) Z$ddZ[ii] = 20

        outfile =  file.path( p$modeldir, paste( "bathymetry", "baseline", domain, "rdata" , sep=".") )

        save (Z, file=outfile, compress=T )
        print( outfile )
      }
      # require (lattice); levelplot( z~plon+plat, data=Z, aspect="iso")
      return( "completed" )
    }




    if (DS=="areal_units") {
      fn = file.path( p$datadir, paste( "areal_units", p$spatial.domain, p$resolution, "rdata", sep=".") )
      sppoly = NULL
      if (!redo) {
        if (file.exists(fn)) load(fn)
        return(sppoly)
      }

      sppoly = areal_units(
        strata_type="aegis_lattice",
        resolution=p$resolution,
        spatial.domain=p$spatial.domain,
        proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
        overlay="none",
        redo=TRUE
      )

      W.nb = poly2nb(sppoly, row.names=sppoly$StrataID, queen=TRUE)  # slow .. ~1hr?
      W.remove = which(card(W.nb) == 0)

      if ( length(W.remove) > 0 ) {
        # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
        W.keep = which(card(W.nb) > 0)
        W.nb = nb_remove( W.nb, W.remove )
        sppoly = sppoly[W.keep,]
        row.names(sppoly) = as.character(sppoly$StrataID)
        sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
        sppoly$StrataID = factor( as.character(sppoly$StrataID) )
        sppoly$strata = as.numeric( sppoly$StrataID )
        sppoly = sppoly[order(sppoly$strata),]
      }

      attr(sppoly, "nb") = W.nb
      save(sppoly, file=fn, compress=TRUE)
      return( sppoly )
    }





    if (DS=="carstm.inputs") {

    }


  }  # end bathymetry.db
