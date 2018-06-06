## a function that builds on gasExchangeFlex.R to allow for Sensitivity testing. Assumptions
##    are set to follow those used for the routines CO2 calculations, and the function now
##    indepdendently calculates the pco2atm. 
## NOTE!!! This version uses k600 values from Vachon and Prairie  dx.doi.org/10.1139/cjfas-2013-0241

## required variables: Temperature Celsius ("temp"), pH ("ph"), Conductivity [uS/cm] ("cond"), 
## Atmospheric pCO2 ppm ("pco2atm"), Wind m/s ("wind"), Salinity ("salt"), lake Area [km2] ("larea")
## "either or" etc variables: Barometric pressure kP ("kpa") & Altitude (m) ("alt");   
##                        DIC uM ("dic") may be calculated from cond for kerri

gasExchangeK <- function(df) { 
  
  temp <- df[,1]
  cond <- df[,2]
  ph <- df[,3]
  wind <- df[,4]
  salt <- df[,5]
  dic <- df[,6]
  kpa <- df[,7]
  pco2atm <- df[,8]
  larea <- df[,9]
  
  r1a <- -2225.22 # all r1 are Khco3...
  r1b <- -0.049
  r1d <- 8.91
  r1e <- 336.6
  
  r2a <- 1346.24 # all r2 are kd
  r2b <- -0.126
  r2d <- -6.44
  r2e <- -196.4
  localc <- 1.056
  
  ## other data for scenario choices
  
  # =IF(G2<>" ",G2*7*0.0000025,AF$2*7*0.0000025)
  ##   i.e. if cond missing, grab mean cond value - but for now we are omitting missing
  ##   values so no mean calcs in this script
  ionstrength <- cond*7*0.0000025 
  
  # =(0.000011*F2^2-0.012*F2+6.58)-(0.5*SQRT(K2))/(1+1.4*SQRT(K2))
  pk1 <- (0.000011*temp^2 - 0.012*temp + 6.58) - (0.5*sqrt(ionstrength))/(1 + 1.4*sqrt(ionstrength))
  
  # =(0.00009*F2^2-0.0137*F2+10.62)-(2*SQRT(K2))/(1+1.4*SQRT(K2))
  pk2 <- (0.00009*temp^2 - 0.0137*temp + 10.62) - (2*sqrt(ionstrength))/(1 + 1.4*sqrt(ionstrength))
  
  # =1/(1+10^(-L2+D2)+10^(-(L2+M2)+2*D2))
  ao <- 1/(1 + 10^(-pk1 + ph) + 10^(-(pk1 + pk2) + 2*ph))
  
  # =1/(10^(-D2+L2)+1+10^(-M2+D2))
  a1 <- 1/(10^(-ph+pk1) + 1 + 10^(-pk2 + ph))
  
  # =1/(10^(-2*D2+(L2+M2))+10^(-D2+M2)+1)
  a2 <- 1/(10^(-2*ph + (pk1 + pk2)) + 10^(-ph + pk2) + 1)
  
  # =1.11+0.016*F2-0.00007*F2^2
  pkh <- 1.11 + 0.016*temp - 0.00007*temp^2
  
  ## =1000000*(E2*0.000001-10^(-14+D2)+10^(-D2))/(O2+2*P2) 
  ## dic <- 1000000*(alk*0.000001 - 10^(-14 + ph) + 10^(-ph))/(a1 + 2*a2) 
  ##   used in case alk exists, not dic. NB equation not tested to be correct, nor was it used by Kerri. 
  ##   She did a different regression-based calculation for missing (or ALL!!!) DIC values using 
  ##   conductivity-DIC relationship:
  
  # =((O4+2*P4)*(R4)-1000000*10^(-D4)+1000000*10^(-14)+D4)
  alk <- ((a1 + 2*a2)*(dic) - 1000000*10^(-ph) + 1000000*10^(-14) + ph) 
  
  # =R2*N2 
  co2 <- dic*ao
  
  # =S2/(10^-Q2) .... mcAtm
  pco2 <- co2/(10^-pkh) 
  
  ## =IF(J2="",+(I2/101.325)*((+H2*0.000001)*10^(-Q2+6)),
  ##      EXP(-0.12806*(+J2/1000))*((+H2*0.000001)*10^(-Q2+6)))
  ## if no alt, use kpa
  ## when both available, use alt. But kpa actually better to default to (discussions with kerri)
  co2eq <- (kpa/101.325) * ((pco2atm * 0.000001) * 10^(-pkh + 6))
  
  
  # log(S2)
  co2log <- log10(co2) 
  
  # =10^(-D2)
  ah <- 10^(-ph)
  
  # =2.07+0.215*F3^1.7 k cm/h ..  2.07 + 0.215*wind^1.7 is now replaced by this:
  # NOTE: unreliable for lakes < 0.1kmsq....approximately. see vachon and prairie
  kspeed1 <- 2.51 + 1.48*wind + 0.39*wind*log10(larea)
  
  # =G2/60/60 k cm/s ..
  kspeed2 <- kspeed1/60/60
  
  # =10^-(0.000152*(C2^2)-0.0129*C2+6.5784)
  k1 <- 10^-(0.000152*(temp^2) - 0.0129*temp + 6.5784)
  
  # =10^-(0.000121*(C2^2)-0.0149*C2+10.626)
  k2 <- 10^-(0.000121*(temp^2) - 0.0149*temp + 10.626)
  
  # =10^-(0.0002*C2^2-0.0444*C2+14.953)
  kw <- 10^-(0.0002*temp^2 - 0.0444*temp + 14.953)
  
  # =(0.61563+0.05316*C2)*0.00001 D cm2/s .. no idea what D is...
  dspeed <- (0.61563+0.05316*temp)*0.00001
  
  # =EXP($X$2+$Y$2*$T2^0.5+(10^4)*$Z$2/(C2+273)+$AA$2*LN(C2+273)) 
  r1 <- exp(r1a + r1b*salt^0.5 + (10^4)*r1d/(temp + 273) + r1e*log(temp + 273)) 
  
  # =EXP($X$3+$Y$3*$T2^0.5+(10^4)*$Z$3/(C2+273)+$AA$3*LN(C2+273))
  r2 <- exp(r2a + r2b*salt^0.5 + (10^4)*r2d/(temp + 273) + r2e*log(temp + 273)) 
  
  # =M2+N2*K2/E2
  r <- r1 + r2*kw/ah 
  
  # =1+E2^2*(I2*J2+I2*E2)^-1
  t <- 1 + ah^2*(k1*k2 + k1*ah)^-1
  
  # =(O2*R2/L2)^0.5
  q <- (r*t/dspeed)^0.5
  
  # =L2/H2
  z <- dspeed/kspeed2
  
  # =R2/((R2-1)+TANH(P2*Q2)/(P2*Q2))
  alpha <- t/((t - 1) + tanh(q*z)/(q*z))
  
  # =(S2-U2)*0.8
  flux <- (co2 - co2eq)*localc # latter requires calculation for local conditions; 
  # one stored by EW is what Kerri used
  
  # =Y2*X2
  fluxenh <- flux*alpha
  #data.frame(fluxenh, pco2)
}

