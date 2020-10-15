
map_area_unit_problem = function( inp = NULL, logz=FALSE, just_return_results=FALSE ) {

  require(aegis)
  require(aegis.bathymetry)

  if (is.null(inp)) {
    p = list()
    p$data_root = project.datadirectory( "aegis", "bathymetry" )
    p$datadir = file.path( p$data_root, "data" )
    fn = file.path(p$datadir, "maup.rdata")
    fn_res = file.path(p$datadir, "maup_summary.rdata")
    if (logz) fn_res = file.path(p$datadir, "maup_log_summary.rdata")
    if (just_return_results) {
      load(fn_res)
      return(maup)
    }
    if (file.exists(fn) ) load( fn )
    if (is.null(inp)) {
      inp = bathymetry_db( p=p, "z.lonlat.rawdata")
      box = list( lon=c(-55,-45), lat=c(40, 50) )
      inp = inp[ which(
        inp$lon > box$lon[1] & inp$lon < box$lon[2] &
        inp$lat > box$lat[1] & inp$lat < box$lat[2]
      ), ]
      inp = lonlat2planar( inp, proj.type=projection_proj4string("utm20"))
      save (inp, file=file.path(p$datadir, "maup.rdata"), compress=TRUE )
    }
  }

  xvalues = inp[,"plon"]
  yvalues = inp[,"plat"]
  zvalues = inp[,3]

  xrange = range( xvalues, na.rm=TRUE )
  yrange = range( yvalues, na.rm=TRUE )
  dmin = min( diff(xrange),  diff(yrange) )

  if (logz) zvalues = log(zvalues+6000)

  results=list()
  results$global = list(
    mean = mean( zvalues, na.rm=TRUE),
    sd = sd(zvalues, na.rm=TRUE),
    n = length(zvalues),
    min = min(zvalues, na.rm=TRUE),
    max = max(zvalues, na.rm=TRUE),
    median = median(zvalues, na.rm=TRUE),
    resolution = sqrt(mean( diff( sort(unique(xvalues))))^2 + mean(diff(sort(unique(yvalues)) ))^2)
  )

  for (discretized_n in c(2, 3, 4, 5, 7, 10, 15, 20, 40, 60, 80, 100, 150, 200, 300, 400, 800, 1000, 1100, 1500, 2000, 5000, 8000, 10000) ) {
    message( discretized_n
    )

    minresolution = rep( dmin / discretized_n, 2)

    # basic aggregation in 2D
    res = tapply(
      X=zvalues,
      INDEX=list(
        aegis_floor( xvalues / minresolution[1] + 1) * minresolution[1],
        aegis_floor( yvalues / minresolution[2] + 1) * minresolution[2]),
        FUN = function(w) {mean(w, na.rm=TRUE)},
        simplify=TRUE
    )
    res = as.data.frame( as.table (res) )
    res[,1] = as.numeric(as.character( res[,1] ))
    res[,2] = as.numeric(as.character( res[,2] ))
    res = res[ which( is.finite( res[,3] )) ,]
    names(res) =c("x", "y", "z")

    results[[as.character(discretized_n)]] = list(
      mean = mean( res[,3], na.rm=TRUE),
      sd = sd(res[,3], na.rm=TRUE),
      n = length(res[,3]),
      min = min(res[,3], na.rm=TRUE),
      max = max(res[,3], na.rm=TRUE),
      median = median(res[,3], na.rm=TRUE),
      resolution = diff(xrange)/discretized_n
    )

  }

  maup = do.call(rbind.data.frame, results)  # list to data frame
  maup = maup[ order(maup$resolution), ]

  require(stmv)
  gr = stmv_variogram( cbind(xvalues, yvalues), zvalues, methods="fft", plotdata=FALSE )

  attr( maup, "variogram" ) = gr

  save(maup, file=fn_res, compress=TRUE)

  if (0) {

    x = maup$resolution
    # x = log(  maup$n / maup$resolution^2  )  # data density
    yrange = range( maup$min, maup$max )
    plot(mean ~ x, maup, pch=20, ylim=yrange)
    lines( median ~ x, maup, col="green", lwd=1)
    lines(min ~ x, maup, col="red", lwd=4)
    lines(max ~ x, maup, col="blue", lwd=4)
    lines( I(mean+sd) ~ x, maup, col="gray", lwd=2)
    lines( I(mean-sd) ~ x, maup, col="gray", lwd=2)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 378 km

    plot( sd~ x, maup)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

  }

  return(maup)
}
