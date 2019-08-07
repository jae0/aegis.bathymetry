
bathymetry.figures = function( p=NULL, varnames="z", logyvar=FALSE, isodepths = c( 100, 300, 500 ), savetofile="", width=1365, height=1024, pointsize=12, res=96, quality=80 ) {

  b = bathymetry.db( p=p, DS="complete" )
  b = b[b$z > 0 , ]

  for (vn in varnames) {
    print(vn)

    if (logyvar) {
      b = b[ b[,vn]>0, ]
      b[,vn] = log(b[,vn])
    }

    # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions",internal.crs=p$internal.crs )
    # if (length(oc) > 0) {
    #   b = b[oc,]
    # }

    datarange = quantile( b[,vn], probs=c(0.05, 0.95), na.rm=TRUE )
    dr = seq( datarange[1], datarange[2], length.out=100)

    levplt = levelplot( b[,vn] ~ plon + plat, data=b[,], aspect="iso", main=NULL,
      at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
      panel = function(x, y, subscripts, ...) {
        panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
        sp.lines( isobath.db( p=p, DS="isobath", depths=isodepths, crs=p$internal.crs ), col = "gray80", cex=0.1 )
        sp.lines( coastline.db( p=p, crs=p$internal.crs ), col = "steelblue", cex=0.1 )
      }
    )

    if ( savetofile != "" ) {
      outdir = file.path(p$data_root, "bathymetry", "maps", p$spatial.domain )
      for (i in 1:length(savetofile)){
        devtype = savetofile[i]
        if (devtype =="jpeg") devtype="jpg"
        fn = file.path( outdir, paste( "bathymetry", vn, p$spatial.domain, devtype, sep=".") )
        print(fn)
        if (devtype == "pdf" ) {
          pdf(file=fn, width=5, height=4, bg='white')
        } else if (devtype == "png" ) {
          png(filename=fn, width=width, height=height, pointsize=pointsize, res=res, bg='white' )
        } else if (devtype == "jpg" ) {
          jpeg(filename=fn, width=width, height=height, pointsize=pointsize, res=res, bg='white', quality=quality )
        } else {
          stop( "device not supported ... add it here")
        }
        print( levplt )
        dev.off()
      }
    } else {
      print( levplt )
    }

  }
  return(varnames)
}
