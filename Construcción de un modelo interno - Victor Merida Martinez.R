library(goftest)

frecuencia = read.csv2("https://raw.githubusercontent.com/UC3M-student/Modelo-Interno/main/Frecuencia.csv", header = TRUE)$Frecuencia
severidad = read.csv2("https://raw.githubusercontent.com/UC3M-student/Modelo-Interno/main/Severidad.csv", header = TRUE)$Severidad

set.seed(1234)
########### FRECUENCIA ################

lambda_frecuencia = mean(frecuencia) # Se valida asi segun IOPA


########### SEVERIDAD  ################

muestra = sort(severidad)

severidad_depurado = muestra[muestra < muestra[2232*0.995]] #Nos quedamos con el 99.5 % de los valores y el 0.5 % no los estudiamos por atipicos


#### Kernel Density Estimation####

kernel_density_estimation = function(muestra, min, max) {
  puntos = seq(min,max,1)
  n = length(muestra)
  h1=1.06*sd(muestra)/(n^(1/5))
  f<-c()
  for (j in 1:length(puntos)){
    f[j]=(1/(h1*n))*sum(dnorm((puntos[j]-muestra)/h1,0,1))
  }
  df_kernel <<-data.frame(cbind(puntos,f))
  plot(df_kernel)
  aux<<-cbind(puntos,cumsum(f)) 
}

kernel_density_estimation(severidad_depurado,0,50000)

#### Estimacion Parametrica ####

# Segun Rene Doff, la cuantia de un siniestro puede seguir una distribucion Gamma, Normal o Log-Normal.

#### 1º Validacion - ¿ Cual es el mejor procedimiento en terminos de Sesgo , MSE y Consistencia ? #####
set.seed(1234)
## Distribucion Gamma
parametro = 300 
poblacion_gamma<-rgamma(10000,parametro,300)

M_sim<-c()
PM_sim<-c()
ML_sim<-c()

for (i in 1:100){ 
  muestra99<-sample(poblacion_gamma,193,replace=TRUE)
  
  #Metodo de los momentos
  param<-c()
  dif<-function(param){ 
    r1<-(param[1]/param[2]-mean(muestra99))^2
    r2<-(param[1]/(param[2]^2)-var(muestra99))^2
    return(r1+r2)
  }
  MM<-optim(c(2,1),dif,method="L-BFGS-B")$par
  
  M_sim[i] = MM[1]
  
  #Metodo Percentile Matching
  
  param<-c()
  dif<-function(param){
    r1<-(qgamma(0.1,param[1],param[2])-quantile(muestra99,0.1))^2
    r2<-(qgamma(0.8,param[1],param[2])-quantile(muestra99,0.8))^2
    return(r1+r2)  
  }
  PM<-optim(c(2,1),dif,method="L-BFGS-B",lower=c(0,0))$par
  PM_sim[i] = PM[1]
  
  #Metodo Maximo Verosimilitud
  
  param<-c()
  LL <- function(param) {
    -sum(dgamma(muestra99, param[1], param[2], log=TRUE)) 
  }
  ML <- optim(c(2,1),LL,method="L-BFGS-B")$par
  ML_sim[i]<-ML[1]
}

# SESGO - Gamma
Bias_M=parametro-mean(M_sim) #Sesgo Metodo Momentos
Bias_PM=parametro-mean(PM_sim) #Sesgo Percentil Matching
Bias_ML=parametro-mean(ML_sim) #Sesgo Maximum Likelihood

# MSE - MEAN SQUARED ERROR - ERROR CUADRATICO MEDIO - Gamma
MSE_M=Bias_M^2+var(M_sim) #MSE Metodo Momentos
MSE_PM=Bias_PM^2+var(PM_sim) #MSE Percentil Matching
MSE_ML=Bias_ML^2+var(ML_sim) #MSE Maximum Likelihood

#CUADRO DE RESULTADOS - Distribucion Gamma
Resultados_Gamma<-matrix(c(Bias_M,Bias_PM,Bias_ML,MSE_M,MSE_PM,MSE_ML),3,2)
rownames(Resultados_Gamma)<-c("Momentos","PM","ML") 
colnames(Resultados_Gamma)<-c("sesgo","MSE")
Resultados_Gamma  

# Consistencia
set.seed(1234)
parametro=130
p=130
N=10000
limite=300
min=20
consistencia_Gamma<-rgamma(N,parametro,p) 
M_sim<-matrix(0,limite-min,1)
PM_sim<-matrix(0,limite-min,1)
ML_sim<-matrix(0,limite-min,1)
bar <- txtProgressBar(0,limite,style=3)
i=1

for (n in min:limite){ 
  muestra99<-sample(consistencia_Gamma,n,replace=FALSE)
  
  #Metodo de los momentos
  dif<-function(param){
    r1<-(param[1]/param[2]-mean(muestra99))^2
    r2<-(param[1]/(param[2]^2)-var(muestra99))^2
    return(r1+r2)
  }
  MM<-optim(c(2,1),dif,method="L-BFGS-B")
  M_sim[i]<-MM$par[1]
  
  #Metodo Percentile Matching
  param<-c()
  
  dif<-function(param){
    r1<-(qgamma(0.1,param[1],param[2])-quantile(muestra99,0.1))^2
    r2<-(qgamma(0.8,param[1],param[2])-quantile(muestra99,0.8))^2
    return(r1+r2)  
  }
  
  PM<-optim(c(2,1),dif,method="L-BFGS-B")
  PM_sim[i]<-PM$par[1]
  
  #Metodo Maximo Verosimilitud
  param<-c()
  LL <- function(param) {
    -sum(dgamma(muestra99, param[1], param[2], log=TRUE))
  }
  ML <- optim(c(2,1),LL,method="L-BFGS-B")
  ML_sim[i]<-ML$par[1]
  setTxtProgressBar(bar, n)
  i=i+1
}

ts.plot(M_sim)
abline(parametro,0,col="red") #Metodo de los momentos
ts.plot(PM_sim)
abline(parametro,0,col="red") #Percentil Matching
ts.plot(ML_sim)
abline(parametro,0,col="red") #Maximum Likelihood

# Resultado: tras evaluar el mejor procedimiento en terminos de Sesgo, MSE y Consistencia -> Metodo Metodo Maximum Likelihood

## Distribucion Normal
set.seed(1234) #Ver si tengo que cambiar el set.seed
parametro = 300 
poblacion_normal<-rnorm(10000,parametro,300)

M_sim<-c()
PM_sim<-c()
ML_sim<-c()

for (i in 1:100){ 
  muestra99<-sample(poblacion_normal,193,replace=TRUE)
  
  #Metodo de los momentos
  M_sim[i] = mean(muestra99)
  
  #Metodo Percentile Matching
  
  param<-c()
  dif<-function(param){
    
    (qnorm(0.5, param[1], param[2]) - quantile(muestra99, 0.5))^2
    
  }
  
  PM<-optim(c(2,1),dif,method="L-BFGS-B",lower=c(0,0))$par
  PM_sim[i] = PM[1]
  
  #Metodo Maximo Verosimilitud
  
  param<-c()
  LL <- function(param) {
    -sum(dnorm(muestra99, param[1], param[2], log=TRUE)) 
  }
  ML <- optim(c(2,1),LL,method="L-BFGS-B")$par
  ML_sim[i]<-ML[1]
  
}

# SESGO - Normal
Bias_M=parametro-mean(M_sim) #Sesgo Metodo Momentos
Bias_PM=parametro-mean(PM_sim) #Sesgo Percentil Matching
Bias_ML=parametro-mean(ML_sim) #Sesgo Maximum Likelihood

# MSE - MEAN SQUARED ERROR - ERROR CUADRATICO MEDIO - Noraml
MSE_M=Bias_M^2+var(M_sim) #MSE Metodo Momentos
MSE_PM=Bias_PM^2+var(PM_sim) #MSE Percentil Matching
MSE_ML=Bias_ML^2+var(ML_sim) #MSE Maximum Likelihood

#CUADRO DE RESULTADOS - Normal
Resultados_Normal<-matrix(c(Bias_M,Bias_PM,Bias_ML,MSE_M,MSE_PM,MSE_ML),3,2)
rownames(Resultados_Normal)<-c("Momentos","PM","ML") 
colnames(Resultados_Normal)<-c("sesgo","MSE")
Resultados_Normal 

#Consistencia
set.seed(1234)
parametro=130
p=130
N=10000
limite=300
min=20
consistencia_Normal<-rnorm(N,parametro,p) 
M_sim<-matrix(0,limite-min,1)
PM_sim<-matrix(0,limite-min,1)
ML_sim<-matrix(0,limite-min,1)
bar <- txtProgressBar(0,limite,style=3)
i=1

for (n in min:limite){ 
  
  muestra99<-sample(consistencia_Normal,n,replace=FALSE)
  
  #Metodo de los momentos
  
  M_sim[i]<- mean(muestra99)
  
  #Metodo Percentile Matching
  param<-c()
  
  dif<-function(param){
    (qnorm(0.5, param[1], param[2]) - quantile(muestra99, 0.5))^2
  }
  
  PM<-optim(c(2,1),dif,method="L-BFGS-B")
  PM_sim[i]<-PM$par[1]
  
  #Metodo Maximo Verosimilitud
  param<-c()
  LL <- function(param) {
    -sum(dnorm(muestra99, param[1], param[2], log=TRUE))
  }
  ML <- optim(c(2,1),LL,method="L-BFGS-B")
  ML_sim[i]<-ML$par[1]
  setTxtProgressBar(bar, n)
  i=i+1
}

ts.plot(M_sim)
abline(parametro,0,col="red") #Metodo de los momentos
ts.plot(PM_sim)
abline(parametro,0,col="red") #Percentil Matching
ts.plot(ML_sim)
abline(parametro,0,col="red") #Maximum Likelihood

# Resultado: tras evaluar el mejor procedimiento en terminos de Sesgo, MSE y Consistencia -> Metodo Maximum Likelihood

## Log - normal 
set.seed(123)
parametro = 1
poblacion_lognormal<-rlnorm(10000,parametro,1)

M_sim<-c()
PM_sim<-c()
ML_sim<-c()

for (i in 1:100){ 
  muestra99<-sample(poblacion_lognormal,193,replace=TRUE)
  
  #Metodo de los momentos
  M_sim[i] = mean(muestra99)
  
  #Metodo Percentile Matching
  param<-c()
  dif<-function(param){
    
    (qlnorm(0.5, param[1], param[2]) - quantile(muestra99, 0.5))^2
    
  }
  
  PM<-optim(c(1,1),dif,method="L-BFGS-B",lower=c(0,0))$par
  PM_sim[i] = PM[1]
  
  #Metodo Maximo Verosimilitud
  
  param<-c()
  LL <- function(param) {
    -sum(dlnorm(muestra99, param[1], param[2], log=TRUE)) 
  }
  ML <- optim(c(2,1),LL,method="Nelder-Mead")$par
  ML_sim[i]<-ML[1]
  
}

# SESGO - Log-Normal
Bias_M=parametro-mean(M_sim) #Sesgo Metodo Momentos
Bias_PM=parametro-mean(PM_sim) #Sesgo Percentil Matching
Bias_ML=parametro-mean(ML_sim) #Sesgo Maximum Likelihood

# MSE - MEAN SQUARED ERROR - ERROR CUADRATICO MEDIO - Log-Normal
MSE_M=Bias_M^2+var(M_sim) #MSE Metodo Momentos
MSE_PM=Bias_PM^2+var(PM_sim) #MSE Percentil Matching
MSE_ML=Bias_ML^2+var(ML_sim) #MSE Maximum Likelihood

#CUADRO DE RESULTADOS - Log-Normal
Resultados_LogNormal<-matrix(c(Bias_M,Bias_PM,Bias_ML,MSE_M,MSE_PM,MSE_ML),3,2)
rownames(Resultados_LogNormal)<-c("Momentos","PM","ML") 
colnames(Resultados_LogNormal)<-c("sesgo","MSE")
Resultados_LogNormal  

# Consistencia
set.seed(1234)
parametro=1
p=1
N=10000
limite=300
min=20

consistencia_LogNormal<-rlnorm(N,parametro,p) 
M_sim<-matrix(0,limite-min,1)
PM_sim<-matrix(0,limite-min,1)
ML_sim<-matrix(0,limite-min,1)
bar <- txtProgressBar(0,limite,style=3)
i=1

for (n in min:limite){ 
  
  muestra99<-sample(consistencia_LogNormal,n,replace=FALSE)
  
  #Metodo de los momentos
  
  M_sim[i]<- mean(log(muestra99))
  
  #Metodo Percentile Matching
  param<-c()
  
  dif<-function(param){
    (qlnorm(0.5, param[1], param[2]) - quantile(muestra99, 0.5))^2
  }
  
  PM<-optim(c(2,1),dif,method="L-BFGS-B")
  PM_sim[i]<-PM$par[1]
  
  #Metodo Maximo Verosimilitud
  param<-c()
  LL <- function(param) {
    -sum(dlnorm(muestra99, param[1], param[2], log=TRUE))
  }
  ML <- optim(c(2,1),LL,method="Nelder-Mead")
  ML_sim[i]<-ML$par[1]
  setTxtProgressBar(bar, n)
  i=i+1
}

ts.plot(M_sim)
abline(parametro,0,col="red") #Metodo de los momentos
ts.plot(PM_sim)
abline(parametro,0,col="red") #Percentil Matching
ts.plot(ML_sim)
abline(parametro,0,col="red") #Maximum Likelihood

# Resultado: tras evaluar el mejor procedimiento en terminos de Sesgo, MSE y Consistencia -> Metodo Maximum Likelihood


### Una vez elegido el metodo Maximum Likelihood para estimar los parametros, ahora identificaremos cual se ajusta mas a los datos


## Distribucion Gamma 
set.seed(1234)
df_gamma = data.frame()

param<-c()
LL_Gamma <- function(param) {
  -sum(dgamma(severidad_depurado, param[1], param[2], log=TRUE)) 
}

ML_Gamma <- optim(c(1,1),LL_Gamma,method="Nelder-Mead")$par
ML_Gamma_1 = ML_Gamma

df_gamma = rbind(df_gamma, c(ML_Gamma[1], ML_Gamma[2], cvm.test(severidad_depurado,"pgamma", ML_Gamma[1], ML_Gamma[2], estimated = TRUE)$statistic,"Gamma"))
colnames(df_gamma) = c("Param_1", "Param_2","CvM","Distribution")

## Distribucion Normal 
set.seed(1234)
df_normal = data.frame()

param<-c()
LL_Normal <- function(param) {
  -sum(dnorm(severidad_depurado, param[1], param[2], log=TRUE)) 
}

ML_Normal <- optim(c(1,1),LL_Normal,method="Nelder-Mead")$par

df_normal = rbind(df_normal, c(ML_Normal[1], ML_Normal[2], cvm.test(severidad_depurado,"pnorm", ML_Normal[1], ML_Normal[2], estimated = TRUE)$statistic, "Normal"))

colnames(df_normal) = c("Param_1", "Param_2","CvM","Distribution")


## Distribucion Log-normal 
set.seed(1234)
df_lnormal = data.frame()

param<-c()
  
LL_LogNormal <- function(param) {
  -sum(dlnorm(severidad_depurado, param[1], param[2], log=TRUE)) 
}

ML_Log_Normal <- optim(c(1,1),LL_LogNormal,method="Nelder-Mead")$par

df_lnormal = rbind(df_lnormal, c(ML_Log_Normal[1], ML_Log_Normal[2], cvm.test(severidad_depurado,"plnorm", ML_Log_Normal[1], ML_Log_Normal[2], estimated = TRUE)$statistic, "Log-Normal"))

colnames(df_lnormal) = c("Param_1", "Param_2","CvM","Distribution")


##Tabla de resultados

df_resultado_cvm = data.frame()

df_resultado_cvm = rbind(df_resultado_cvm, df_gamma, df_lnormal, df_normal)
df_resultado_cvm$Param_1 = round(as.numeric(df_resultado_cvm$Param_1),4)
df_resultado_cvm$Param_2 = round(as.numeric(df_resultado_cvm$Param_2),7)
df_resultado_cvm$CvM = round(as.numeric(df_resultado_cvm$CvM),4)
df_resultado_cvm$Distribution = c("Gamma","Log-Normal", "Normal")

df_resultado_cvm

### Resultado: por el criterio de valoracion de la estimacion de densidades de Cramer y Von Mises, elegimos la distribucion Gamma

## Para apoyar nuestra postura,usamos una 2º validacion: Cross Section Validation


#### 2º Validacion: Cross Section Validation  ####
set.seed(1234)
n<-length(severidad_depurado)
sim=100

Dif_Gamma = c()
Dif_Norm = c()
Dif_Lnorm = c()

for (rep in 1:sim){
  prop=0.9
  selec <- sample(length(severidad_depurado),size = n*prop) 
  train<-severidad_depurado[selec] 
  val<-severidad_depurado[-selec] 
  
  LL <- function(param) {
    -sum(dgamma(train, param[1], param[2], log=TRUE))
  }
  
  ML <- optim(c(1,1),LL,method="Nelder-Mead")
  
  prediccion_Gamma<-sort(rgamma(length(val),ML$par[1], ML$par[2]))
  
  
  LL <- function(param) {
    -sum(dnorm(train, param[1], param[2], log=TRUE)) 
  }
  
  ML <- optim(c(1,1),LL,method="Nelder-Mead")
  prediccion_Normal<-sort(rnorm(length(val),ML$par[1], ML$par[2]))
  
  
  LL <- function(param) {
    -sum(dlnorm(train, param[1], param[2], log=TRUE)) 
  }
  
  ML <- optim(c(1,1),LL,method="Nelder-Mead")
  prediccion_Lnormal<-sort(rlnorm(length(val),ML$par[1], ML$par[2]))
  
  
  val_ordenados<-sort(val) 
  Dif_Gamma[rep] = mean((prediccion_Gamma-val_ordenados)^2)
  Dif_Norm[rep] = mean((prediccion_Normal-val_ordenados)^2)
  Dif_Lnorm[rep] = mean((prediccion_Lnormal-val_ordenados)^2)
  
}

Resultados_Cross_Validation = matrix(c(mean(Dif_Gamma),mean(Dif_Norm),mean(Dif_Lnorm)),1,3)
rownames(Resultados_Cross_Validation)<-c("MSE") 
colnames(Resultados_Cross_Validation)<-c("Gamma","Normal","LogNormal")

Resultados_Cross_Validation

# Resultado: La distribucion Gamma tiene el minimo MSE

#### DISTRIBUCION DEL ESTIMADOR ####

##Lo calculamos por bootstrap

set.seed(1234)

Shape_bootstrap = c()
Rate_bootstrap = c()

for (i in 1:100) {
  
  muestra_bootstrap = sample(severidad_depurado,length(severidad_depurado), replace = TRUE)
  
  param<-c()
  LL_Gamma <- function(param) {
    -sum(dgamma(muestra_bootstrap, param[1], param[2], log=TRUE)) 
  }
  
  ML_Gamma <- optim(c(1,1),LL_Gamma,method="Nelder-Mead")$par
  
  Shape_bootstrap[i] = ML_Gamma[1]
  Rate_bootstrap[i] = ML_Gamma[2]
  
}

c(quantile(Shape_bootstrap,0.025), quantile(Shape_bootstrap,0.975)) #Bandas de confianza de Shape
c(quantile(Rate_bootstrap,0.025), quantile(Rate_bootstrap,0.975)) #Bandas de confianza de Rate

#### Bandas de confianza - CvM (Ajuste de la densidad) ###### 


CvM_Gamma_distribution = c() 

for (i in 1:100) {  
  
  CvM_Gamma_distribution[i] = cvm.test(severidad_depurado,"pgamma",Shape_bootstrap[i],Rate_bootstrap[i], estimated = TRUE)$statistic
  
}

c(quantile(CvM_Gamma_distribution,0.025),quantile(CvM_Gamma_distribution,0.975)) #Intervalo de Confianza CvM

#### MODELO AGREGADO ####
set.seed(123456)
simulaciones=100

coste_cartera<-c()
for (sim in 1:simulaciones){
  coste_poliza<-c()
  for (i in 1:23236){ 
    frecuencia<-rpois(1,lambda_frecuencia)
    severidad<-rgamma(frecuencia,ML_Gamma_1[1],ML_Gamma_1[2])
    coste_poliza[i]<-sum(severidad) 
  }
  
  coste_cartera[sim]<-sum(coste_poliza) 
}

hist(coste_poliza)
hist(coste_cartera)


#### VaR 99.5 % ####
var_cartera_modelo_agregado<-quantile(coste_cartera,0.995)

#### TVaR ####
TVaR_cartera_modelo_agregado = mean(coste_cartera[coste_cartera>var_cartera_modelo_agregado])

#### Best-Estimate ####
Best_Estimate_cartera_modelo_agregado = mean(coste_cartera) 

#### CAPITAL ECONOMICO ####
SCR_cartera_modelo_agregado = var_cartera_modelo_agregado - mean(coste_cartera)


########### TABLAS DE MORTALIDAD DE JAPON - 19 ###########

mortalidad = read.csv2("https://raw.githubusercontent.com/UC3M-student/Modelo-Interno/main/Mortalidad.csv", header = TRUE)
mortalidad$qx = as.numeric(mortalidad$qx)

# La tabla de mortalidad es estimada por el metodo de la distribucion empirica y, por tanto, no hacen falta nuevos calculos de estimacion

set.seed(1234444)

suma_asegurada = c()

for (w in 1:100) {
  
  total_fallecidos = 0
  
  for (i in 1:99) {
    
    Prob_morir = mortalidad$qx[i] 
    
    for (z in 1:mortalidad$nx[i]) {
      
      Monte_Carlo_number = runif(1)
      
      if (Monte_Carlo_number <= Prob_morir ) {
        
        total_fallecidos = total_fallecidos + 1
        
      }
    }
  }
  
  suma_asegurada[w] = total_fallecidos*1200000
  
}

hist(suma_asegurada)

#### VaR 99.5 % ####
VaR_995_mortalidad = quantile(suma_asegurada,0.995)

#### TVaR ####
TVaR_mortalidad = mean(suma_asegurada[suma_asegurada > VaR_995_mortalidad])

#### Best-Estimate ####
Best_Estimate_mortalidad = mean(suma_asegurada) #REPASAR SI ESTA BIEN ; Se puede dejar asi porque sigue una distribucion Normal

#### CAPITAL ECONOMICO ####
SCR_mortalidad = VaR_995_mortalidad - Best_Estimate_mortalidad


############# RIESGOS TOTALES DE LA COMPAÑIA

# Usamos el Teorema de Fisher -> suma de distribuciones Normales sigue una Normal


Riesgo_Compañia = sort(suma_asegurada) + sort(coste_cartera)

#### VaR 99.5 % ####
VaR_995_compañia = quantile(Riesgo_Compañia,0.995)

#### TVaR ####
TVaR_compañia = mean(Riesgo_Compañia[Riesgo_Compañia > VaR_995_compañia])

#### Best-Estimate ####
Best_Estimate_Compañia = mean(Riesgo_Compañia)

#### CAPITAL ECONOMICO ####
SCR_compañia = VaR_995_compañia - Best_Estimate_Compañia


