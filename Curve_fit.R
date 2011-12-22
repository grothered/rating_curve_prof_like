## R code to fit a velocity curve to the montalban data.

# Read in the data, which was pasted into data_montalban.txt from Excel
measurements=read.delim(file='data_montalban.txt',header=F)

# Optionally, we add names to the columns in measurements
names(measurements) = c('Date', 'Site', 'Time', 'Stage', 'Distance', 'Duration', 'Speed', 'Mean Speed', 'Cross-sectional Area', 'Discharge')


# NOTE: On July 23, 2004, there are 2 measurements -- one which an anomolously
# high drouge time. We will have to drop this measurement.
measurements = measurements[measurements$Duration<100,]

# Plot measurements, with the date overlayed
plot(measurements[,4], measurements[,8], xlab='Stage', ylab='Mean Speed', col='steelblue', pch=19)
text(measurements[,4], measurements[,8], labels=measurements[,1], cex=0.6)


# Try fitting a relation:
# Speed ~ a*depth^b.
# using nonlinear least squares. 
# Assume depth = Stage - 20.5 (Seems a reasonable guess).
speed=measurements[,8]
stage=measurements[,4]
curve1=nls(speed ~ a*(stage-20.5)^b, start=list(a=0.5, b=0.5)) # This is the curve fitting routine
summary(curve1)
confint(curve1) # Compute confidence intervals for a and b

pdf(file='Montalban_speed_velocity_fit.pdf', width=8.0, height=5.0)
plot(stage, speed, xlab='Stage', ylab='Mean Speed',pch=19,col='steelblue', xlim=c(20.5,30), ylim=c(0,4))
points(seq(20.5, 30,len=20), 0.8162*( (seq(20.5,30,len=20)-20.5)^0.6919),t='l',col=2)
legend('topleft', c('Data', 'v = 0.82*(Stage - 20.5)^0.69'), lty=c(0,1), lwd=c(0,1), col=c('steelblue', 'red'), pch=c(19,NA))
dev.off()
# Note that b ~= 2/3
# This is consistent with a Manning relation: Speed = (n/sqrt(Sf))*(depth)^(2/3)
#
# We may approximate Sf as the water slope, ~= 5/8700 based on comparing water
# levels at Sto Nino and Nangka (which are separated by about 8700m) 
# Substituting this into the equation suggests n~=0.02 if a=0.8
# This seems physically reasonable, i.e. within the typical range of
# manning's coefficient 

curve2=nls(speed ~ a*(stage-20.5)^0.67, start=list(a=0.5)) # This is the curve fitting routine

# We can combine this with a Stage- Area relation to get a rating curve.

# To get a stage-area relation that is also valid for high flows, we need data extending up to that range.

#xsect=read.delim('montalban_xsection_2006.txt')
xsect=read.table('Composite_cross_section.txt')
xsect[,2]=xsect[,2]+10.47 # Add MSL to make the script consistent with previous versions

cs_area<-function(Stage){

    # Approximate the cross-section using 1000 equally spaced points
    f<-approx(xsect[,1], xsect[,2], n=1000)

    output=Stage*NA
    # Numerical integration to get cross-sectional area.
    for(i in 1:length(Stage)){
        s=Stage[i]
        output[i] = sum(pmax(s-f$y,0.0))*(f$x[2]-f$x[1])    
    }

    output
}

discharge_2006<-function(Stage){
    cs_area(Stage)*0.82*(Stage-20.5)^0.69
}


# Read in the montalban data
mont_WL = read.csv(file='montWL.csv', header=F, na.strings='no_data')
in2006 = which(mont_WL[,2]==2006)
in2006_2011=which((mont_WL[,2]==2006)|(mont_WL[,2]==2007)|(mont_WL[,2]==2008)|(mont_WL[,2]==2009)|(mont_WL[,2]==2010)|(mont_WL[,2]==2011))

# Compute discharge in 2006
mont_Q_2006 = discharge_2006(mont_WL[in2006, 6])
mont_Q_2006_2011 = discharge_2006(mont_WL[in2006_2011, 6])

# Write out the stage and discharge in 2006
write.table(cbind(mont_WL[in2006,1:6], mont_Q_2006), file='montalban_2006.csv',row.names=F,col.names=F, quote=F, sep=",") 
write.table(cbind(mont_WL[in2006_2011,1:6], mont_Q_2006_2011), file='montalban_2006_2011.csv',row.names=F,col.names=F, quote=F, sep=",") 



### Let's consider confidence intervals for the regression.
### Variables are 'speed' and 'depth' = stage -20.5
depth=stage-20.5
# If speed=a*depth^(b)*epsilon then log(speed) = log(a) + b*log(depth) + log(epsilon)
# Although this assumes a particular error structure
logspeed=log(speed)
logdepth=log(depth)

mod_log=lm(logspeed~logdepth)

# Compute prediction confidence intervals
drange=log(seq(0.1,15,len=15))
log_prediction=predict(mod_log, data.frame(logdepth=drange), interval='prediction')

# Plot
plot(depth, speed, ylim=c(0,10))
a=exp(coef(mod_log)[1])
b=coef(mod_log)[2]

points(depth, a*depth^b,t='l')
points(exp(drange), exp(log_prediction[,2]),t='l',col=2)
points(exp(drange), exp(log_prediction[,3]),t='l',col=2)

#########################################################################################
#
# Profile likelihood methods
#
#
# Best to fit the non-linear model, compute the confidence limits by
# exhaustively searching over the parameter bounds
#
#

# NOTE: Measurements 1:3 on August 4, 2005 are inconsistent with the EFCOS water levels. So I will cut them out of this.
include_pts=7:length(measurements[,1])

# Select points with < 10% error in the area -- see below for justification
# include_pts=which( abs(measurements[,9]/cs_area(measurements[,4]) -1.0) < 0.1)
# Note -- these points are not independent. 
# WITH THE ABOVE CHOICE, WE CANNOT REASONABLY CONSTRAIN THE PARAMETERS. THE ONLY THING WE COULD DO IS DO A FIT OF speed = a*(stage-constant)^b

power_loglik<-function(pars){
    a=pars[1]
    b=pars[2]
    sigma=pars[3]
    # Log likelihood
    sum(dnorm(speed[include_pts], mean= a*depth[include_pts]^b ,sd= sigma, log=T))
    }

# Compute parameter values
myfit=optim(par=c(1,1,1),fn=power_loglik, control=list(fnscale=-1))

# Plot
drange=seq(0,9,len=20)
plot(depth,speed,pch=19,col='steelblue', ylim=c(0,4))
points(drange, myfit$par[1]*drange^(myfit$par[2]),t='l')

# Compute confidence intervals by exhaustively computing the liklihood region
n_increm=200
arange=seq(0.5,1.5,len=n_increm)
brange=seq(0.2,1.,len=n_increm)
sigrange=seq(0.05,0.6,len=n_increm)
loglik_out=array(NA,dim=c(n_increm, n_increm, n_increm))
for(i in 1:n_increm){
    print(i)
    for(j in 1:n_increm){
        for(k in 1:n_increm){
        loglik_out[i,j,k] = power_loglik(c(arange[i], brange[j], sigrange[k]))
        }
    }
}
# According to Bolker, approx confidence intervals come from searching within
# the zone where the log-liklihood is within 2 units of the max
# Alternatively, the likelihood ratio test that Bolker reports gives the following: 
#CI_zone=loglik_out>=(max(loglik_out)- qchisq(0.95,1)/2) # EQN BOLKER

# The above result disagrees with R's confint function. 
# Other texts cite the following formula for the liklihood region FOR SINGLE PARAMETERS
df=length(include_pts)-2
CI_zone_onevar= (loglik_out>=(max(loglik_out)- qf(0.95, 1, df)/2) ) # EQN Donaldson & Schnabel (1985)
# Note that in the limit as the dataset size --> inf, this approaches the
# formula given in Bolker. See the useful facts below

# Donaldson and Schnabel (1985) report that the CONFIDENCE
# REGION (where all parameters are allowed to vary), is
# defined by: 
#
# S(theta) - S(theta_hat) <= (S(theta_hat)/(n-p))*p*qf(1-alpha, p, n-p) 
#
# where S(theta) is the residual sum of squares evaluated at parameter values 'theta',
# 'theta_hat' is the vector of parameter values that minimise S,
# p is  the number of parameters in the distribution, 
# n is the number of data points,
# and alpha is the confidence level (e.g. 5%).
# Note also that S(theta_hat)/(n-p) is an estimate of the residual variance of the model fit 
# The above result can be written as:
# -2*(loglik(theta)-loglik(theta_hat)) <= p*qf(1-alpha, p, n-p)

# On the other hand, the confidence interval about theta_j, a particular variable in theta, is
# the range of values of theta_j which:
#
# maximise( (theta_j - theta_j_hat)^2 ) subject to
# S(theta) - S(theta_hat) <=  S(theta_hat)/(n-p) * qf(1-alpha, 1, n-p)
#
# Note that this formula leads to decent agreement with R.
# It can be written as:
# -2*(loglik(theta) - loglik(theta_hat)) <= qf(1-alpha, 1, n-p)
 
# Compute single variable confidence limits
CI_inds=which(CI_zone_onevar==TRUE,arr.ind=T)
a_ci=arange[range(CI_inds[,1])]
b_ci=brange[range(CI_inds[,2])]
sig_ci=sigrange[range(CI_inds[,3])]

# Compare with nls
v1=speed[include_pts]
v2=depth[include_pts]
mod1=nls(v1~a*v2^b, start=list(a=1,b=1))
confint(mod1)

# Compute joint confidence limits
# EQN Donaldson & Schnabel (1985), and Bates and Watts (1988)
CI_zone_allvar=(loglik_out>=(max(loglik_out)- 2*qf(0.95, 2, df-2)/2) )
CI_inds_all=which(CI_zone_allvar==TRUE,arr.ind=T)
a_ci_all=arange[range(CI_inds_all[,1])]
b_ci_all=brange[range(CI_inds_all[,2])]
sig_ci_all=sigrange[range(CI_inds_all[,3])]

# Useful Facts
# Note: qchisq(prob, 1) = lim(df1-->inf) F(prob,1,df1)
#       df1*lim(df2--> inf) F(prob,df1, df2) = qchisq(prob, df1) 
#        sum( rnorm(k, mean=0, sd=1)^2) has a chisq(x, k) distribution 
#       ( t(n-p, 1-a/2))^2 has an F(1, n-p, 1-a) distribution


drange=seq(0,9,len=80)
# It is expensive to plot all the liklihood points -- sample them instead
speed_upper=drange*NA
speed_lower=drange*NA
for(i in 1:length(drange) ){
    print(i)
    speed_range=range( arange[CI_inds_all[,1]]*drange[i]^brange[CI_inds_all[,2]] ) 
    #speed_range=range( arange[CI_inds[,1]]*drange[i]^brange[CI_inds[,2]] ) 
    speed_lower[i]=speed_range[1]
    speed_upper[i]=speed_range[2]
}

#pdf(file='Montalban_stage_speed_proflik.pdf', width=5.4,height=5.4)
plot(depth, speed, ylim=c(0.1,5),main='Montalban Rating Curve data, UNREALISTIC')
points(drange, speed_upper,t='l',col=2)
points(drange,speed_lower, t='l',col=2)
points(drange, myfit$par[1]*drange^(myfit$par[2]),t='l', col='blue')
points(depth[include_pts], speed[include_pts], pch=19)
text(depth, speed, labels=measurements[,1], cex=0.6)
#dev.off()
# For prediction intervals, see Vecchia and Cooley, 1987. 
# Obviously the data shows evidence of time changes in the rating curve.
# However, how about we assume that the observations characterise the
# uncertainties in the curve that we would expect to experience in applying it.

# Question: Is it plausible that the difference in the July 2006 and August 2004 observations could be due to changes in cross-sectional area?
pdf('Area_fitted_and_observed.pdf', width=5.4, height=5.4)
plot(measurements[,9], cs_area(measurements[,4]),log='xy')
points(1:500,1:500,t='l',col=2)
text(measurements[,9], cs_area(measurements[,4]), measurements[,1], cex=0.6)
points(measurements[include_pts,9], cs_area(measurements[include_pts,4]),col=2,pch=19)
dev.off()
# So either 1) The data is very wrong, or 2) There have been changes in the cross-section. However, there do not seem to be many changes since 2006, in that part of the cross-section captured in the LiDAR. In either case, perhaps it is best to only fit our model to the data that agrees with the modern measurements? [NOTE: If the area estimates were wrong, but the velocity measurements were still correct, then it would be okay to use the old data. However, if the cross-section has just changed, then it would not be okay to use the old data. So the 'safe' thing is to just use the more recent data].

# Select points with < 10% error in the area
keepers=which(abs(measurements[,9]/cs_area(measurements[,4]) -1.0) < 0.1)
# This is so few points that we cannot constrain the error in our curve. The only thing we can do is fit the following.
v_1=speed[keepers]
v_2=measurements$Stage[keepers]
mod2=nls(v_1~a*(v_2-bed)^0.67, start=list(a=1, bed=20))
# This is a crap fit. Suggests that the roughness is really high at low Stages, which is plausible (vegetation etc in the river?). OR, that we have a looped rating curve, which is significantly unsteady. Check water slope / stage relation at this site.
# NOTE: Analysis suggests that although the flow is significantly unsteady, the variation in sqrt(surface slope) (max-min)/min should be only around 20%. We expect the discharge to roughly scale with this. 

# Alternative: Interpolation?
#vel_stage=aggregate(c(0, v_1), list(c(20.1, v_2)), mean)
#velfun=approxfun(vel_stage[,1], vel_stage[,2])
# This is unlikely to be much good.

# Based on analysis of the water levels at montalban and Nangka, we expect the water slope during a flood to varying around (5-7)/5000 or so, due to unsteadyness. This might change the discharge by 20%, assuming velocity and slope^0.5 are linearly related (which won't be true during unsteady flow!)


# Try working directly with discharge (so assuming their xsectional area is correct)
pdf('Rating_curve_performance.pdf', width=10,height=7)
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,3,0))
# Plot of 'measured' with native cross-sectional area, and fitted using 2006 cross-sectional area
Q_tmp = measurements[,8]*measurements[,9]
plot(measurements[include_pts,4], Q_tmp[include_pts], xlab='Stage (m)', ylab='Discharge (m^3/s)', log='y', xlim=c(20,28), pch=19)
points(seq(20.5,27,len=40), discharge_2006(seq(20.5,27,len=40)),t='l',col=2)
text(measurements[include_pts,4], Q_tmp[include_pts],measurements[include_pts,1])
legend('bottomright', c('Measured', 'Rating_curve'), pch=c(19,NA), lty=c(NA,1), col=c(1,2))

# Relative error
plot( measurements[include_pts,4], Q_tmp[include_pts]/discharge_2006(measurements[include_pts,4]), xlab='Stage (m)', ylab='(Measured Q)/(Rating_curve_Q)', log='y', xlim=c(20,28), pch=19)
abline(h=1)
text(measurements[include_pts,4], Q_tmp[include_pts]/discharge_2006(measurements[include_pts,4]), measurements[include_pts,1])
abline(h=0.5,lty='dashed')
abline(h=2,lty='dashed')

title('Rating curve, measurements and errors @ Montalban', cex.main=2, outer=T)
dev.off()
