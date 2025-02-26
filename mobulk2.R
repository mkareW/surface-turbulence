# rm(list=ls())
 Zv=2
 Zt=Zv
 Zq=Zv
 Zrh=Zv
z0<-0.001
z0t<-0.001
z0q<-0.0001

V<-5
Ta<-283.15
RH<- 30
p<-530
  Ts<-273.15 #
# # # 
#   Ts<-on_gla$Tsurf[i]+273.15
#   V<-on_gla$WSPD1[i]
#   p<-on_gla$PRES[i]
#   RH<-on_gla$RH[i]
#   Ta<-on_gla$TAIR[i]

#V windspeed
#Ta Tair in Kelvin
#Ts temperature surface in Kelvin
#RH air humidity %
#p barometric pressure mbrs
#Zv,Zt,Zrh, height of measurements of the different V,T, RH (m)
#z0,z0t,z0q: surface roughness lengths

MObulk <- function(V,Ta,RH,p,Zv,Zt,Zrh,Ts,z0,z0t,z0q) { 
#V windspeed in ms-1
#Ta air temperature in K
#RH air relative humidity in %
#p air barometric pressure in mbars
#Zv,Zt,Zrh : height of the measurements (V,T,RH) in m 
#Ts: temperature of the surface in K
#z0,z0t,z0q: surface roughness lengths in m

#choice of stability functions set  1=monin-obukhov, 2=brutsaert
func<-2

#stability functions definitions 
# momentum unstable
psimun<-function(monin){
x<-(1-16*monin)^(1/4)
psimunstabl<-2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2
return(psimunstabl)}

#heat unstable
psihun <- function(monin){
x<-(1-16*monin)^(1/4)
psihunstabl<-2*log((1+x^2)/2)
return(psihunstabl)} 

#neutral
neut<- function(monin){
neutra<-0*monin
return(neutra)}

#brutsaert heat stable
psihst2 <- function(monin){
psihstabl2<- -5*monin
return(psihstabl2)}

#momentum unstable
psimst2 <- function(monin){
psimstabl2<- -5*monin
return(psimstabl2)}

#heat super stable
psihsst2 <- function(monin){ 
psihsstabl2<- -5*(log(monin)+1)
return(psihsstabl2)}

#momentum super stable
psimsst2 <- function(monin){ 
psimsstabl2<- -5*(log(monin)+1)
return(psimsstabl2)}


g<-9.81 #gravity acceleration ms-2
k<-0.4  #constant, von-karmann
LL1<-2830000   #sublimation heat J/kg;
LL2<-2500000   #vaporization heat J/Kg;

Tapot=Ta*(1024/p)^0.286 #potential air Temperature
Tspot=Ts*(1024/p)^0.286 #potential surface temperature

#Latent heat of - sublimation or evaporation
LL<-LL1
LL[which(Ts>=273.15)]<-LL2


#Air density
rho<-1.29*p/1024;

#Specific humidity calculation

#RH extrapolated to the height of V measurements, assuming log profile (neutral) 
alpha<-(RH-100)/(log(Zrh/z0q))
RHz<-100+alpha*log(Zt/z0q)

Zrh<-Zt 
#partial pressure of water in the air mbars
e<-(RHz/100)*6.1078*exp(17.08085*(Ta-273.15)/(234.175+Ta-273.15))
#partial pressure of water at the surface (considered saturated- hypothesis to check) mbars
es<-6.1078*exp(17.08085*(Ts-273.15)/(234.175+Ts-273.15))
      
q<-e*0.622/p  # spécific humidity of the air kg/kg
qs<-es*0.622/p #spécific humidity of the surface kg/kg

Ta_v=Ta*(1+0.61*q) # virtual T air 
Ts_v=Ts*(1+0.61*qs)# virtual T surf

Tav_pot=Ta_v*(1024/p)^0.286 # virtual potential T air
Tsv_pot=Ts_v*(1024/p)^0.286 # virtual potential T surf 
  
Tapot_v=Tapot*(1+0.61*q) #potential virtual Tair
Tspot_v=Tspot*(1+0.61*qs) #potential virtual Ts

  
#air specific heat
Cp=1005*(1+0.84*q)


#boucle flux
zsurL<-NA*V
if (Zv>z0 & Zt>z0t & Zrh>z0q & V>0 & is.na(V)==FALSE & Ta!=Ts)
{flag<-1} else {flag<-0}


if (flag==1) {       

ustar_i<-k*V/log(Zv/z0)
tstar_i<-k*(Ta-Ts)/log(Zt/z0t)
qstar_i<-k*(q-qs)/log(Zrh/z0q)

}

if (flag==1) {    
 Lstar<-(Ta_v*ustar_i^2)/(k*g*(tstar_i+0.61*qstar_i*Ta))
 zsurL<-Zv/Lstar}


#debut de la boucle
if (flag==1&is.na(zsurL)!=TRUE) {ite<-1
             delta<-1000
while (ite<51 & delta>0.01) {

if (zsurL<0&is.na(zsurL)!=TRUE) {                                                    
  stab_funct_m<-psimun
  stab_funct_h<-psihun 
}
  
  if (zsurL==0&is.na(zsurL)!=TRUE) {                                                    
    stab_funct_m<-neut
    stab_funct_h<-neut
  }

if (zsurL>0&is.na(zsurL)!=TRUE) {
stab_funct_m<-psimst2
stab_funct_h<-psihst2 
}

if (zsurL>=1&is.na(zsurL)!=TRUE)  {
stab_funct_m<-psimsst2
stab_funct_h<-psihsst2 
}


ustar<-k*(V)/(log(Zv/z0)-stab_funct_m(Zv/Lstar)+stab_funct_m(z0/Lstar))  
tstar<-k*(Ta-Ts)/(log(Zt/z0t)-stab_funct_h(Zt/Lstar)+stab_funct_h(z0t/Lstar))
qstar<-k*(q-qs)/(log(Zrh/z0q)-stab_funct_h(Zrh/Lstar)+stab_funct_h(z0q/Lstar))

Lstarold<-Lstar
Lstar<-(Ta_v*ustar^2)/(k*g*(tstar+0.61*qstar*Ta))
delta=abs((Lstar-Lstarold))

# show(i)
# show(V)
# show(Ta)
# show(Ts)
# show(RH)
# show(zsurL)
# show(Lstar)
# points(ite+1,Lstar)
ite<-ite+1

}

#flux calculation
H<-rho*Cp*ustar*tstar;#*(p/1024)^0.28;
LE<-LL*rho*ustar*qstar;

#calcul de l'erreur//erreurs de mesures;  //paper AMTD
dlnz0<-2.24#2.6//1.5;
dlnz0t<-2.51#2.6//1.5;
dV<-0.1
dT<-0.1
dTs<-0.35
dZ<-0.1
dRH<-0.03

erreurV<-dV/V
erreurT<-dTs/abs(Ta-Ts)
erreurQ<-dRH/RH+dT*234.175*17.08085/(234.17+Ta)^2
erreurQs<-dRH/RH+dTs*234.175*17.08085/(234.17+Ts)^2
erreurdQ<-sqrt(erreurQ^2+erreurQs);
erreurz0T<-(dlnz0t)*1/(log(Zt/z0t)-stab_funct_h(z0t/Lstar))
erreurz0<-(dlnz0)*1/(log(Zv/z0)-stab_funct_m(z0t/Lstar))
erreurZ<- -(dZ/Zv)*(stab_funct_m(z0t/Lstar)-log(Zv/z0))^(-1)*(stab_funct_h(z0t/Lstar)-log(Zt/z0t))^(-1)
erreurtot_sqrH=sqrt(erreurz0*erreurz0+erreurz0T*erreurz0T+erreurV*erreurV+erreurZ*erreurZ+erreurT*erreurT); 
erreurtot_sqrLE=sqrt(erreurz0*erreurz0+erreurz0T*erreurz0T+erreurV*erreurV+erreurZ*erreurZ+erreurdQ*erreurdQ); 

} else { H<- NA
LE<- NA
erreurV<- NA
erreurT<- NA
erreurQ<- NA
erreurQs<- NA
erreurdQ<- NA
erreurz0T<- NA
erreurz0<- NA
erreurZ<- NA
erreurtot_sqrH<- NA 
erreurtot_sqrLE<- NA
ustar<- NA
tstar<- NA
qstar<- NA
zsurL<- NA
}


 return(c(H,LE,ustar,tstar,qstar,zsurL,erreurtot_sqrH,erreurtot_sqrLE,erreurV,erreurT,erreurQ,erreurQs,erreurdQ,erreurz0T,erreurz0,erreurZ))
}

