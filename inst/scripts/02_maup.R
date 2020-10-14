

### demonstration of the map area unit problem example:

maup = map_area_unit_problem( just_return_results=TRUE )  #default is bathymetry data

    x = maup$resolution
    # x = log(  maup$n / maup$resolution^2  )  # data density
    yrange = range( maup$min, maup$max )
    plot(mean ~ x, maup, pch=20, ylim=yrange)
    lines( median ~ x, maup, col="green", lwd=1)
    lines(min ~ x, maup, col="red", lwd=4)
    lines(max ~ x, maup, col="blue", lwd=4)
    lines( I(mean+sd) ~ x, maup, col="gray", lwd=2)
    lines( I(mean-sd) ~ x, maup, col="gray", lwd=2)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

    plot( sd~ x, maup)
    abline( v=attr(maup, "variogram")$fft$localrange ) # 380 km

