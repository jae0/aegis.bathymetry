
# ------------------------------------------------
# Areal unit modelling of bathymetry


reset_input_data = FALSE
if (0) reset_input_data = TRUE # choose this if we are redoing input data "views"


# --------------------------------
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p = aegis.bathymetry::bathymetry_parameters(
  project.mode="carstm",
  id = "bathymetry",
  resolution = 20, # km
  internal.crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
  internal.crs_planar = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
  boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
  libs = RLibrary ( "sp", "spdep", "rgeos", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry")
)



p$variables$Y = "z"  # name to give variable in extraction and model



# set up default map projection
p = c(p, coastline_layout( p=p, redo=reset_input_data ) )
p$mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")


# --------------------------------
# Get the data


if (exists(tmpdir)) {
  setwd( tmpdir )  # temp files saved here
} else {
  tmpdir = getwd()
}

# --------------------------------
# ensure if polys exist and create if required
# for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygons_managementarea( species="snowcrab", au))

sppoly = areal_units(
  strata_type="aegis_lattice",
  resolution=p$resolution,
  spatial.domain=p$spatial.domain,
  proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  overlay="none",
  redo=reset_input_data
)

W.nb = neighbourhood_structure( sppoly=sppoly, strata_type="aegis_lattice" )



## ----------------------------------
# covariate estimates for prediction in strata and year
# collapse PS vars with time into APS (and regrid via raster)
  APS = aegis_db_extract(
    vars=p$lookupvars,
    yrs=p$yrs,
    spatial.domain=p$spatial.domain,
    dyear=p$prediction.dyear,
    resolution=p$resolution,
    proj4string=sp::CRS(p$internal.crs),
    returntype="data.frame",
    redo = reset_input_data
  )

  APS$yr = as.numeric( APS$year)
  APS$Y = NA
  APS$data_offset = 1  # force to be density n/km^2
  APS$tag = "predictions"


  # StrataID reset to be consistent in both data and prediction areal units
  o = over( SpatialPoints( set[,c("plon", "plat")], sp::CRS(p$internal.crs) ), spTransform(sppoly, sp::CRS(p$internal.crs) ) ) # match each datum to an area
  set$StrataID = o$StrataID


  o = over( SpatialPoints( APS[,c("plon", "plat")], sp::CRS(p$internal.crs) ), spTransform(sppoly, sp::CRS(p$internal.crs) ) ) # match each datum to an area
  APS$StrataID = o$StrataID

  o = NULL

  #  good data
  ok = which(
    is.finite(set[,p$variables$Y]) &   # INLA can impute Y-data
    is.finite(set$data_offset) &
    is.finite(set$StrataID)
  )


# construct meanweights matrix
weight_year = meanweights_by_strata( set=set, StrataID=as.character( sppoly$StrataID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
# weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
# weight_year = weight_year[ match(as.character(sppoly$StrataID), rownames(weight_year) )]


varstokeep = unique( c( "Y", "StrataID", "yr", "data_offset", "tag", p$lookupvars) )

M = rbind( set[ok, varstokeep], APS[,varstokeep] )

M = M[ which(
      is.finite(M$data_offset) &
      is.finite(M$StrataID)
    ) , ]

M$yr_factor = factor( as.character(M$yr) )
M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) )
M$strata  = as.numeric( M$StrataID)
M$year  = as.numeric( M$yr_factor)
M$iid_error = 1:nrow(M) # for inla indexing for set level variation


M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$tsd[!is.finite(M$tsd)] = median(M$tsd, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$tmin[!is.finite(M$tmin)] = median(M$tmin, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$tmax[!is.finite(M$tmax)] = median(M$tmax, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$degreedays[!is.finite(M$degreedays)] = median(M$degreedays, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$dZ[!is.finite(M$dZ)] = median(M$dZ, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$ddZ[!is.finite(M$ddZ)] = median(M$ddZ, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$substrate.grainsize[!is.finite(M$substrate.grainsize)] = median(M$substrate.grainsize, na.rm=TRUE )  # missing data .. quick fix .. do something better


dd_t = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
dd_tsd = seq( min(M$tsd, na.rm=TRUE), max(M$tsd, na.rm=TRUE), length.out=15 )
dd_tmin = seq( -2, 10, by=1 ) # seq( min(M$tmin, na.rm=TRUE), max(M$tmin, na.rm=TRUE), length.out=15 )
dd_tmax = seq( 2, 14, by=1 ) # seq( min(M$tmax, na.rm=TRUE), max(M$tmax, na.rm=TRUE), length.out=15 )
dd_dday = c(10, 100, 200, 400, 800, 1000, 2000, 3000, 4000, 5000)
dd_z = c(2.5, 5, 10, 20, 40, 80, 120, 160, 220, 320, 480, 640, 1000 )
dd_dz = c(0.01, 0.1, 1, 2, 4, 6, 8, 10, 12)
dd_ddz = c(0.01, 0.1, 0.2, 0.4, 0.8, 1, 2, 3, 4)
dd_sstr = c( 0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 4.0, 8.0 )


M$ti = as.numeric( as.character( cut( M$t, breaks=dd_t, labels=diff(dd_t)/2 + dd_t[-length(dd_t)], include.lowest=TRUE ) ))
M$tisd = as.numeric( as.character( cut( M$tsd, breaks=dd_tsd, labels=diff(dd_tsd)/2 + dd_tsd[-length(dd_tsd)], include.lowest=TRUE ) ))
M$timin = as.numeric( as.character( cut( M$tmin, breaks=dd_tmin, labels=diff(dd_tmin)/2 + dd_tmin[-length(dd_tmin)], include.lowest=TRUE ) ))
M$timax = as.numeric( as.character( cut( M$tmax, breaks=dd_tmax, labels=diff(dd_tmax)/2 + dd_tmax[-length(dd_tmax)], include.lowest=TRUE ) ))
M$tidday = as.numeric( as.character( cut( M$degreedays, breaks=dd_dday, labels=diff(dd_dday)/2 + dd_dday[-length(dd_dday)], include.lowest=TRUE ) ))
M$zi = as.numeric( as.character( cut( M$z, breaks=dd_z, labels=diff(dd_z)/2 + dd_z[-length(dd_z)], include.lowest=TRUE ) ))
M$zid = as.numeric( as.character( cut( M$dZ, breaks=dd_dz, labels=diff(dd_dz)/2 + dd_dz[-length(dd_dz)], include.lowest=TRUE ) ))
M$zidd = as.numeric( as.character( cut( M$ddZ, breaks=dd_ddz, labels=diff(dd_ddz)/2 + dd_ddz[-length(dd_ddz)], include.lowest=TRUE ) ))
M$si = as.numeric( as.character( cut( M$substrate.grainsize, breaks=dd_sstr, labels=diff(dd_sstr)/2 + dd_sstr[-length(dd_sstr)], include.lowest=TRUE ) ))


totest = setdiff(1:ncol(M), which(names(M) %in% c("Y", "StrataID", "tag", "yr_factor") ))
ii = which(is.finite(rowSums(M[,totest])))

M = M[ii,]

# ---------------------
# generic PC priors
m = log( {set$Y / set$data_offset}[ok] )
m[!is.finite(m)] = min(m[is.finite(m)])




carstm_hyperparameters = function( reference_sd, alpha=0.5, reference_mean=0 ) {
  # some generic PC priors, scaled by sd of data
  # pc.prior to median .. minimally info. scale

  hyper = list(

    iid = list(
      prec = list(
        prior = "pc.prec",  # exponential decay
        param = c(reference_sd, alpha)
      )
    ),

    # means informative, sd marginally diffuse
    # see: inla.set.control.fixed.default() for defaults
    fixed = list(
        mean.intercept = reference_mean,
        prec.intercept = 1e-3,
        mean=0,
        prec=1e-2
    ),


    # param=c(u, alpha); u=sigma; alpha=prob;
    # see inla.doc("pc.rw2") inla.doc("pc.prec")  ..prior sd attributable to rw2
    rw2 = list(
      prec = list(
        prior = "pc.prec",  # exponential decay
        param = c(reference_sd, alpha)
      )
    ),

    # see inla.doc("ar1") ; theta0, theta1 are expected
    # param=c(u, alpha); u=sigma; alpha=prob;
    # see inla.doc("pc.prec")  ..prior sd attributable to autocor rho
    # param=c(u, alpha); rho = 0.5; u=sqrt(1-rho); alpha=prob; see inla.doc("pc.cor1")
    ar1 = list(
      prec = list(
        prior = "pc.prec",  # exponential decay
        param = c(reference_sd, alpha)
      ),
      rho = list(
        prior = "pc.cor0", # inla.doc("pc.cor0") ..base model: rho = 0  --- expoential; will tend to 0 unless there is info
        param = c(sqrt(1-0.5), 0.1)  # rho=0.5; u=sqrt(1-rho)  ... 100-10% of probablity weight to rho 0.5 or less .. forces smooth and only goes high if really high
      )
    ),

    # param=c(u, alpha); u=phi (proportion spatial); alpha=prob
    bym2 = list(
      prec = list(
        prior = "pc.prec",
        param = c(reference_sd, alpha)
      ),
      phi = list(
        prior="pc",  # see bottom of inla.doc("bym2")
        param=c(0.5, 0.5) # c(phi=0.5, alpha=0.5)
      )
    )
  )

  return(hyper)
}


H = carstm_hyperparameters( sd(m), alpha=0.5, median(m) )
# H$prec$prec.intercept = 1e-9




# -------------------------------------
# simple glm
fit = glm(
  formula = Y ~ 1 + offset( log( data_offset) ) + StrataID + yr_factor,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M[ which(M$tag=="observations"), ]
)

s = summary(fit)
AIC(fit)  # 77326

# reformat predictions into matrix form
ii = which(
  M$tag=="predictions" &
  M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] &
  M$yr_factor %in% M[ which(M$tag=="observations"), "yr_factor"]
)

preds = predict( fit, newdata=M[ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

out = reformat_to_matrix(
  input = preds$fit,
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
# convert numbers/km to biomass/strata (kg)..
RES$poisson_glm = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_glm_cfanorth = colSums( {out * weight_year * sppoly$cfanorth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_glm_cfasouth = colSums( {out * weight_year * sppoly$cfasouth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_glm_cfa4x = colSums( {out * weight_year * sppoly$cfa4x_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km

plot( poisson_glm ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")
plot( poisson_glm_cfanorth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
plot( poisson_glm_cfasouth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
plot( poisson_glm_cfa4x ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")

# map it
vn = "pred"
yr = "2018"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# -------------------------------------
# simple gam
require(mgcv)
fit = gam(
  formula = Y ~ 1 + offset( log( data_offset) ) + StrataID + yr_factor + s(strata, year, bs="ts"),
  family = "poisson", # "zeroinflatedpoisson0",
  data= M[ which(M$tag=="observations"), ]
)

s = summary(fit)
AIC(fit)  # 76752

# reformat predictions into matrix form
ii = which(
  M$tag=="predictions" &
  M$StrataID %in% M[ which(M$tag=="observations"), "StrataID"] &
  M$yr_factor %in% M[ which(M$tag=="observations"), "yr_factor"]
)

preds = predict( fit, newdata=M[ii,], type="response", na.action=na.omit, se.fit=TRUE )  # no/km2

out = reformat_to_matrix(
  input = preds$fit,
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
# convert numbers/km to biomass/strata (kg)..
RES$poisson_gam = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_gam_cfanorth = colSums( {out * weight_year * sppoly$cfanorth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_gam_cfasouth = colSums( {out * weight_year * sppoly$cfasouth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_gam_cfa4x = colSums( {out * weight_year * sppoly$cfa4x_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km

plot( poisson_gam ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")
plot( poisson_gam_cfanorth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
plot( poisson_gam_cfasouth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
plot( poisson_gam_cfa4x ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")

# map it
vn = "pred"
yr = "2018"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# - -----------------------------
# simple with default priors
fit = inla(
  formula = Y ~ 1 + offset( log( data_offset) ) + StrataID + yr_factor + f(iid_error, model="iid", hyper=H$iid) ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
 # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  # num.threads=4,
  # blas.num.threads=4,
  verbose=TRUE
)

fn_test ="~/tmp/snowcrab_20km_bym.Rdata"
# save( fit, file=fn_test)
# load(fn_test)



s = summary(fit)
s$dic$dic  # 32072
s$dic$p.eff # 4865


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_basic = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE )  / 10^6 # kt
RES$poisson_basic_cfanorth = colSums( {out * weight_year * sppoly$cfanorth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_basic_cfasouth = colSums( {out * weight_year * sppoly$cfasouth_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km
RES$poisson_basic_cfa4x = colSums( {out * weight_year * sppoly$cfa4x_surfacearea}, na.rm=TRUE ) / 10^6  # 10^6 kg -> kt # kg/km * km

lines( poisson_basic ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")
lines( poisson_basic_cfanorth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
lines( poisson_basic_cfasouth ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")
lines( poisson_basic_cfa4x ~ yr, data=RES, lty=1, lwd=2.5, col="green", type="b")


# --------------------------------------
# simple with priors
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(strata, model="iid", hyper=H$iid)
      + f(year, model="iid", hyper=H$iid )
      + f(iid_error, model="iid", hyper=H$iid)
   ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
  #  # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  num.threads=4,
  blas.num.threads=4,
  verbose=TRUE
)

fn3 ="~/tmp/snowcrab_30km_iid_space_time.Rdata"
# save( fit, file=fn3)
# load(fn3)

# Fixed effects:
#               mean   sd 0.025quant 0.5quant 0.975quant   mode kld
# (Intercept) -6.269 1.12     -8.468   -6.269     -4.071 -6.269   0

# Random effects:
#   Name	  Model
#     strata BYM2 model
#    year IID model
#    iid_error IID model
#    ti RW2 model
#    tidday RW2 model
#    zi RW2 model
#    zid RW2 model
#    zidd RW2 model
#    si RW2 model

# Model hyperparameters:
#                           mean    sd 0.025quant 0.5quant 0.975quant   mode
# Precision for strata     0.738 0.017      0.707    0.738      0.772  0.736
# Phi for strata           0.328 0.007      0.314    0.328      0.343  0.328
# Precision for year      18.798 0.425     17.974   18.794     19.645 18.788
# Precision for iid_error  1.140 0.021      1.099    1.140      1.182  1.140
# Precision for ti         3.352 0.070      3.216    3.351      3.493  3.349
# Precision for tidday     3.524 0.102      3.329    3.523      3.728  3.520
# Precision for zi         3.465 0.076      3.322    3.462      3.621  3.453
# Precision for zid        4.847 0.111      4.622    4.850      5.058  4.864
# Precision for zidd      54.613 1.207     52.279   54.599     57.025 54.572
# Precision for si         0.054 0.002      0.051    0.054      0.057  0.054

# Expected number of effective parameters(stdev): 5008.73(2.07)
# Number of equivalent replicates : 1.52

# Deviance Information Criterion (DIC) ...............: 32072.22
# Deviance Information Criterion (DIC, saturated) ....: 13129.59
# Effective number of parameters .....................: 4865.00

# Marginal log-Likelihood:  -23435.97


s = summary(fit)
s$dic$dic  # 32072
s$dic$p.eff # 4865


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_basic2 = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) /10^6  # km
lines( poisson_basic2 ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")


# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# -----------------------------------------
# car simple each posterior config takes @30km:: 10sec .. x 25 configs = 4 min  //  @20km :: 40 sec x 25 config = 20 min // @ 10 km 123.90s tot 45 min
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(strata, model="bym2", graph=W.nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      + f(iid_error, model="iid", hyper=H$iid)
    ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  num.threads=2,
  blas.num.threads=2,
  verbose=TRUE
)


fn80 ="~/tmp/snowcrab_30km_bym_iid.Rdata"
# save( fit, file=fn80)
# load(fn80)

s = summary(fit)
s$dic$dic  # 31225
s$dic$p.eff # 5200

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car_simple = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) / 1e6

lines( poisson_car_simple ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# -----------------------------------------
# simple car grouped by year  27 configs x 100s each = 60 min
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(strata, model="bym2", graph=W.nb, group=year, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      + f(iid_error, model="iid", hyper=H$iid)
      # + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tisd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timin, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timax, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tidday, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zid, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zidd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(si, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2),
      ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=FALSE
)

fn90 ="~/tmp/snowcrab_30km_bym|yr.Rdata"
# save( fit, file=fn90)
# load(fn90)

s = summary(fit)
s$dic$dic  # 31107
s$dic$p.eff # 5124

# Fixed effects:
#              mean     sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 5.564 0.1096      5.348    5.564       5.78 5.565   0

# Random effects:
# Name	  Model
#  strata   BYM2 model
# year   IID model
# iid_error   IID model

# Model hyperparameters:
#                           mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for strata    0.3016 0.0408     0.2274   0.2996     0.3879 0.2963
# Phi for strata          0.9970 0.0051     0.9827   0.9991     1.0000 0.9999
# GroupRho for strata     0.9176 0.0146     0.8864   0.9185     0.9434 0.9201
# Precision for year      6.3591 2.4395     2.9280   5.9114    12.3501 5.1296
# Precision for iid_error 0.7329 0.0203     0.6955   0.7320     0.7751 0.7291

# Expected number of effective parameters(std dev): 5278.80(16.21)
# Number of equivalent replicates : 1.399

# Deviance Information Criterion (DIC) ...............: 31107.11
# Deviance Information Criterion (DIC, saturated) ....: 12918.07
# Effective number of parameters .....................: 5124.12

# Marginal log-Likelihood:  -18060.75

# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car.year_iid_yr = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) / 10^6
plot( poisson_car.year_iid_yr ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )








# -----------------------------------------
# simple car grouped by year  45 configs x 22s each = 60 min // 25 min @ 20 km
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(strata, model="bym2", graph=W.nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      + f(iid_error, model="iid", hyper=H$iid)
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tisd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timin, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(timax, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(tidday, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
       + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zid, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(zidd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      # + f(si, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") , # to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=TRUE
)

fn100 ="~/tmp/snowcrab_30km_bym_envir_temp.Rdata"
# save( fit, file=fn100)
# load(fn100)

s = summary(fit)
s$dic$dic  # 30288
s$dic$p.eff # 4942

# Fixed effects:
#              mean     sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 1.324 0.0679      1.191    1.324      1.458 1.324   0

# Random effects:
# Name	  Model
#  strata   BYM2 model
# year   IID model
# iid_error   IID model
# ti   RW2 model
# zi   RW2 model

# Model hyperparameters:
#                            mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for strata     0.3665 0.0582     0.2686   0.3606     0.4970 0.3479
# Phi for strata           0.4514 0.1155     0.2228   0.4554     0.6655 0.4796
# Precision for year       0.5738 0.1993     0.2793   0.5421     1.0536 0.4841
# Precision for iid_error  0.6356 0.0168     0.6032   0.6354     0.6693 0.6349
# Precision for ti        12.5146 9.8082     2.0243   9.9645    38.1867 5.4997
# Precision for zi         0.1644 0.1040     0.0375   0.1412     0.4284 0.0958

# Expected number of effective parameters(std dev): 5348.15(21.85)
# Number of equivalent replicates : 1.381

# Deviance Information Criterion (DIC) ...............: 30288.36
# Deviance Information Criterion (DIC, saturated) ....: 12099.32
# Effective number of parameters .....................: 4941.85

# Marginal log-Likelihood:  -23563.46
# Posterior marginals for linear predictor and fitted values computed


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car.year_iid_yr = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) / 10^6
plot( poisson_car.year_iid_yr ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )






# -----------------------------------------
# simple car grouped by year nn configs x nn s each = nn min
fit = inla(
  formula =
    Y ~ 1 + offset( log( data_offset) )
      + f(strata, model="bym2", graph=W.nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
      + f(year, model="iid", hyper=H$iid )
      + f(iid_error, model="iid", hyper=H$iid)
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(tisd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(timin, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(timax, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(tidday, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zid, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zidd, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(si, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      ,
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list(int.strategy="eb") , # to get empirical Bayes results much faster.
   # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=TRUE
)

fn101 ="~/tmp/snowcrab_30km_bym_envir_all.Rdata"
# save( fit, file=fn101)
# load(fn101)

s = summary(fit)
s$dic$dic  # 30288
s$dic$p.eff # 4942

# Fixed effects:
#              mean     sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 1.324 0.0679      1.191    1.324      1.458 1.324   0

# Random effects:
# Name	  Model
#  strata   BYM2 model
# year   IID model
# iid_error   IID model
# ti   RW2 model
# zi   RW2 model

# Model hyperparameters:
#                            mean     sd 0.025quant 0.5quant 0.975quant   mode
# Precision for strata     0.3665 0.0582     0.2686   0.3606     0.4970 0.3479
# Phi for strata           0.4514 0.1155     0.2228   0.4554     0.6655 0.4796
# Precision for year       0.5738 0.1993     0.2793   0.5421     1.0536 0.4841
# Precision for iid_error  0.6356 0.0168     0.6032   0.6354     0.6693 0.6349
# Precision for ti        12.5146 9.8082     2.0243   9.9645    38.1867 5.4997
# Precision for zi         0.1644 0.1040     0.0375   0.1412     0.4284 0.0958

# Expected number of effective parameters(std dev): 5348.15(21.85)
# Number of equivalent replicates : 1.381

# Deviance Information Criterion (DIC) ...............: 30288.36
# Deviance Information Criterion (DIC, saturated) ....: 12099.32
# Effective number of parameters .....................: 4941.85

# Marginal log-Likelihood:  -23563.46
# Posterior marginals for linear predictor and fitted values computed


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_matrix(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( StrataID=M$StrataID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( StrataID=sppoly$StrataID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$poisson_car.year_iid_yr = colSums( {out * weight_year * sppoly$sa_strata_km2}, na.rm=TRUE ) / 10^6
plot( poisson_car.year_iid_yr ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly@data[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )






# ------------------------------------------------
# Model 1: a binomial (presence-absence) aks, habitat probability model with linear covariate effects

set$Y = trunc(set$Y)

fit = glm(
  formula = Y ~  offset(data_offset) + 1 + z + dZ + ddZ + t + tsd, #+ tmin + tmax + degreedays + log(substrate.grainsize),
  family=poisson(link="log"),
  data=set[ok,]
)


fit = glm(
  formula = Y ~ 1 + z + dZ + ddZ + t + tsd + tmin + tmax + degreedays + log(substrate.grainsize),
  family=binomial(link="logit"),
  data=set[ok,]
)

s = summary(fit)
AIC(fit)  # 20686688

toplot = expand.grid( t=seq(-1,20), z=c(5,10,20,40,80,160,320,640),degreedays=seq(0, 5000, by=100) )
toplot$predictions = predict(fit, newdata=toplot, type="response", se.fit=FALSE  )

plot( predictions ~ z, toplot[ which( {toplot$t==min(toplot$t)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ t, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ degreedays, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$t==min(toplot$t)} ), ], type="b")


# predicted probabilities of observing cod, given covariates (temperature, depth, etc)

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)
APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$StrataID, sppoly$StrataID )

predictions = predict(fit, APS, type="response", se.fit=TRUE  )
APS$predictions = predictions$fit
APS$predictions.se = predictions$se.fit

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$StrataID), ncol=length(p$yrs), dimnames=list( sppoly$StrataID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$predictions
RES$habitat_glm = colSums( {out * sppoly$sa_strata_km2 }, na.rm=TRUE ) /sum(sppoly$sa_strata_km2) # sa weighted average prob habitat


# map it onto strata means of temperature and depth
aps = APS[ APS$year==2017,  ]
iy = match( as.character(sppoly$StrataID), aps$StrataID )
vn = "pred"
sppoly@data[,vn] = APS$predictions[iy]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# ------------------------------------------------
# Model 5b: habitat model with a smoothed covariate effect
fit = mgcv::gam(
  formula = pa ~  1 + s(t, k=3, bs="ts")  +  s(z, k=3, bs="ts")  +  s(degreedays, k=3, bs="ts"),
  family=binomial(link="logit"),
  data=set[ok,]
)

fit = mgcv::gam(
  formula = Y ~ offset(log(data_offset)) + 1 + StrataID:yr_factor + StrataID + yr_factor + s(t, bs="tp", k=3)  + s(z, bs="tp", k=3),
  family=poisson(link="log"),
  data=set[ok, ]
)

s = summa

inverse.logit = function( x ) {
  # x should be the log odds ratio
  oddsratio = exp(x)
  prob = oddsratio / (1 + oddsratio )
  return (prob)
}


plot(fit, all.terms=TRUE, trans=inverse.logit, seWithMean=TRUE, jit=TRUE, rug=TRUE )
s = summary(fit)
AIC(fit)  #  10579 .. a little better than the glm

toplot = expand.grid( t=seq(-1,20), z=c(5,10,20,40,80,160,320,640),degreedays=seq(0, 5000, by=100) )
toplot$predictions = predict(fit, newdata=toplot, type="response", se.fit=FALSE  )

plot( predictions ~ z, toplot[ which( {toplot$t==min(toplot$t)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ t, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ degreedays, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$t==min(toplot$t)} ), ], type="b")

# predicted probabilities of observing cod, given temperature and depth

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)
APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$StrataID, sppoly$StrataID )

predictions = predict(fit, APS, type="response", se.fit=TRUE  )
APS$predictions = predictions$fit
APS$predictions.se = predictions$se.fit

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$StrataID), ncol=length(p$yrs), dimnames=list( sppoly$StrataID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$predictions

RES$habitat_gam = colSums( {out * sppoly$sa_strata_km2 }, na.rm=TRUE ) /sum(sppoly$sa_strata_km2) # sa weighted average prob habitat


# map it
iy = which( APS$year=="2017")
it = match( as.character(sppoly$StrataID), APS$StrataID[iy] )
vn = "pred"
sppoly@data[,vn] = APS$predictions[iy][it]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


# ------------------------------------------------
# Model 5c: as above but via INLA  ... very slow
## similar to GAM

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "StrataID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$yr_factor = factor( M$yr, levels=p$yrs )
M$StrataID  = factor( M$StrataID, levels=levels(sppoly$StrataID ))

dtemps = seq( min(M$t, na.rm=TRUE), max(M$t, na.rm=TRUE), length.out=13 )
ddepths = c(2.5, 5, 10, 20, 40, 80, 160, 320, 640 )
ddegdays = c(10, 100, 200, 400, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000)

M$ti = as.numeric( as.character( cut( M$t, breaks=dtemps, labels=diff(dtemps)/2 + dtemps[-length(dtemps)], include.lowest=TRUE ) ))
M$zi = as.numeric( as.character( cut( M$z, breaks=ddepths, labels=diff(ddepths)/2 + ddepths[-length(ddepths)], include.lowest=TRUE ) ))
M$di = as.numeric( as.character( cut( M$degreedays, breaks=ddegdays, labels=diff(ddegdays)/2 + ddegdays[-length(ddegdays)], include.lowest=TRUE ) ))


fit = inla(
  formula = pa ~ 1
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec)
      + f(di, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)
plot(fit )
s = summary(fit)
s$dic$dic  # 10585   .. not sure why ..
s$dic$p.eff # 17.73

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )
APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$StrataID, sppoly$StrataID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$StrataID), ncol=length(p$yrs), dimnames=list( sppoly$StrataID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$mean
RES$habitat_inla = colSums( {out * sppoly$sa_strata_km2 }, na.rm=TRUE ) /sum(sppoly$sa_strata_km2) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly@data[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# NOTE most of the decline occured when "habitat" was stable !
# NOTE habitat has become more variable since 1995
# NOTE habitat has declined in 2017
# to do : compute SE's and add to the graph



## -------------------------------------------------------
# bym, iid_year
# Model XX:


APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "StrataID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$iid_error = 1:nrow(M)
M$yr_factor = factor( as.character(M$yr) )
M$StrataID  = factor( M$StrataID, levels=levels(sppoly$StrataID ))
M$strata  = as.numeric( M$StrataID)
M$year  = as.numeric( M$yr_factor)

dtemps = seq( min(M$t, na.rm=TRUE), max(M$t, na.rm=TRUE), length.out=13 )
ddepths = c(2.5, 5, 10, 20, 40, 80, 160, 320, 640 )
ddegdays = c(10, 100, 200, 400, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000)

M$ti = as.numeric( as.character( cut( M$t, breaks=dtemps, labels=diff(dtemps)/2 + dtemps[-length(dtemps)], include.lowest=TRUE ) ))
M$zi = as.numeric( as.character( cut( M$z, breaks=ddepths, labels=diff(ddepths)/2 + ddepths[-length(ddepths)], include.lowest=TRUE ) ))
M$di = as.numeric( as.character( cut( M$degreedays, breaks=ddegdays, labels=diff(ddegdays)/2 + ddegdays[-length(ddegdays)], include.lowest=TRUE ) ))

fit = inla(
  formula = pa ~ 1
  + f(iid_error, model="iid", hyper=H$iid)
  + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(di, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(year, model="iid", hyper=H$iid)
  + f(strata, model="bym2", graph=W.nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)
plot(fit )
s = summary(fit)
s$dic$dic  #8885
s$dic$p.eff #163.8

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )

APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$StrataID, sppoly$StrataID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$StrataID), ncol=length(p$yrs), dimnames=list( sppoly$StrataID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$mean
RES$habitat_bym_yriid = colSums( {out * sppoly$sa_strata_km2 }, na.rm=TRUE ) /sum(sppoly$sa_strata_km2) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly@data[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



if (0) {
 fn = file.path( "~", "tmp", "RES.rdata" )
 # save(RES, file=fn)
 # load(fn)
}



dev.new(width=11, height=7)
# col = c("slategray", "turquoise", "darkorange", "green", "blue", "darkred", "cyan", "darkgreen", "purple" )
col = c( "darkorange", "green", "blue",  "cyan" )
pch = c(20, 21, 22, 23)#, 24, 25, 26, 27, 20)
lty = c(1, 3, 4, 5)# , 6, 7, 1, 3, 4 )
lwd = c(4, 4, 4, 4)# , 4, 4, 4, 4, 4 )
type =c("l", "l", "l", "l")#, "l", "l", "l", "l", "l")
legend=c("GLM", "GAM", "INLA", "INLA CAR")#, "INLA Envir AR1 CAR", "INLA Envir AR1 CAR|year", "INLA Envir AR1|strata CAR", "INLA Envir AR1|strata CAR|year", "INLA Envir CAR|year")

plot( habitat_glm  ~ yr, data=RES, lty=lty[1], lwd=lwd[1], col=col[1], pch=pch[1], type=type[1], ylim=c(0.375,0.825), xlab="Year", ylab="kg")
lines( habitat_gam ~ yr, data=RES, lty=lty[2], lwd=lwd[2], col=col[2], pch=pch[2], type=type[2])
lines( habitat_inla ~ yr, data=RES, lty=lty[3], lwd=lwd[3], col=col[3], pch=pch[3], type=type[3])
lines( habitat_bym_yriid ~ yr, data=RES, lty=lty[4], lwd=lwd[4], col=col[4], pch=pch[4], type=type[4])  # yr_iid
# lines( INLA.Envir.AR1.CAR ~ yr, data=RES, lty=lty[5], lwd=lwd[5], col=col[5], pch=pch[5], type=type[5])
# lines( INLA.Envir.AR1.CAR_year ~ yr, data=RES, lty=lty[6], lwd=lwd[6], col=col[6], pch=pch[6], type=type[6])
# lines( INLA.Envir.AR1_strata.CAR ~ yr, data=RES, lty=lty[7], lwd=lwd[7], col=col[7], pch=pch[7], type=type[7])
# lines( INLA.Envir.AR1_strata.CAR_year ~ yr, data=RES, lty=lty[8], lwd=lwd[8], col=col[8], pch=pch[8], type=type[8])
# lines( INLA.Envir.yr_iid.CAR_year ~ yr, data=RES, lty=lty[9], lwd=lwd[9], col=col[9], pch=pch[9], type=type[9])



legend("topright", legend=legend, lty=lty, col=col, lwd=lwd )


# end
