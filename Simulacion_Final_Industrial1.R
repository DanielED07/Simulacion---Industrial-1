#---------------------------------------------------------------#
# LAS SIGUIENTE FUNCION AYUDA A RECALCULAR LIMITES Y LA MUESTRA #
#---------------------------------------------------------------#
data_recalc <- function(datos,idx,n,C4){
  
  datosCombinados <- datos[,-idx] # Se elimina la muestra fuera de control
  # Se recalculan los estadisticos para las muestras
  if(is.null(dim(datosCombinados))){
    xbars <- mean(datosCombinados) # Medias en las  muestras
    x_barbar <- mean(xbars) # X barra barra
    
    ss <- sd(datosCombinados) # Desviaciones estandar de las  muestras
    sbar <- mean(ss) # S barra
  }else{
    xbars <- apply(datosCombinados, 2, mean) # Medias en las  muestras
    x_barbar <- mean(xbars) # X barra barra
    
    ss <- apply(datosCombinados, 2, sd) # Desviaciones estandar de las  muestras
    sbar <- mean(ss) # S barra
  }
  
  # Sigmas para x y s
  sigma_x <- sbar/(sqrt(n)*C4)
  sigma_s <-(sbar/C4)*sqrt(1-C4^2)
  limiteControlInferiorXbar <- LCI(cl=x_barbar,sigma = sigma_x,L=L)
  limiteControlSuperiorXbar <- LCS(cl=x_barbar,sigma = sigma_x,L=L)
  limiteControlInferiorSbar <- LCI(cl=sbar,sigma = sigma_s,L=L)
  limiteControlSuperiorSbar <- LCS(cl=sbar,sigma = sigma_s,L=L)
  
  return(list("datosCombinados"=datosCombinados,"xbars"=xbars,"ss"=ss,"sbar"=sbar,"x_barbar"=x_barbar,
              "limiteControlInferiorXbar"=limiteControlInferiorXbar,"limiteControlSuperiorXbar"=limiteControlSuperiorXbar,
              "limiteControlInferiorSbar"=limiteControlInferiorSbar,"limiteControlSuperiorSbar"=limiteControlSuperiorSbar))
}
#--------------------------------------#
# FUNCION PARA REPLICAR EL EXPERIMENTO #
#--------------------------------------#
F1_per <- function(L,m,n,c,mu1){
require(zeallot)
set.seed(123)
L=1;m=25;n=5;mu1=3;c=0.04

numeroSubgruposContaminados <- ceiling(m*c)
numeroSubgruposLimpios <- m-numeroSubgruposContaminados

xLimpios <- rnorm(n*numeroSubgruposLimpios,mean=0,sd=1)
xContaminados <- rnorm(n*numeroSubgruposContaminados,mean=mu1,sd=1)

nombresContaminados <- paste0("C",1:numeroSubgruposContaminados)
nombresLimpios <- paste0("L",1:numeroSubgruposLimpios)

datosContaminados <- data.frame(matrix(rep(0,n*numeroSubgruposContaminados),ncol=numeroSubgruposContaminados)) # Datos fuera de control
datosLimpios <- data.frame(matrix(rep(0,n*numeroSubgruposLimpios),ncol=numeroSubgruposLimpios)) # Datos en control

colnames(datosContaminados) <- nombresContaminados
colnames(datosLimpios) <- nombresLimpios

for (col in 1:numeroSubgruposContaminados) {
  datos <- sample(xContaminados, n, replace = FALSE) # Muestra aleatoria de datos generados
  datosContaminados[,col] <- datos # Se agrega la muestra aleatoria al data.frame
  idx <- match(datos,xContaminados) # Indices de la muestra aleatoria
  xContaminados <- xContaminados[-idx] # Se eliminan de los datos generados
  
}

for (col in 1:numeroSubgruposLimpios) {
  datos <- sample(xLimpios, n, replace = FALSE) # Muestra aleatoria de datos generados
  datosLimpios[,col] <- datos # Se agrega la muestra aleatoria al data.frame
  idx <- match(datos,xLimpios) # Indices de la muestra aleatoria
  xLimpios <- xLimpios[-idx] # Se eliminan de los datos generados
  
}

datosCombinados <- cbind(datosLimpios,datosContaminados)
#apply(datosCombinados, 2,mean)

xbars <- apply(datosCombinados, 2, mean) # Medias en las m muestras
x_barbar <- mean(xbars) # X barra barra
 
ss <-apply(datosCombinados, 2, sd) 
sbar <- mean(ss) # S barra

C4 <- c4(n) # c4
sigma_s <-(sbar/cte)*sqrt(1-C4^2)
sigma_x <- sbar/(sqrt(n)*C4)


#L <- 2
limiteControlInferiorXbar <- LCI(cl=x_barbar,sigma = sigma_x,L=L)
limiteControlSuperiorXbar <- LCS(cl=x_barbar,sigma = sigma_x,L=L)
limiteControlInferiorSbar <- LCI(cl=sbar,sigma = sigma_s,L=L)
limiteControlSuperiorSbar <- LCS(cl=sbar,sigma = sigma_s,L=L)

par(mfrow=c(1,2))

plot(x=1:ncol(datosCombinados),y=ss,type="line",ylim=c(0,2))
abline(h=c(limiteControlInferiorSbar,limiteControlSuperiorSbar),col="red")
abline(h=sbar,col="blue")

plot(x=1:ncol(datosCombinados),y=xbars,type="line",ylim=c(-3,3))
abline(h=c(limiteControlInferiorXbar,limiteControlSuperiorXbar),col="red")
abline(h=x_barbar,col="blue")

criterioParada=TRUE
numeroIteraciones=0
Ts <- 0
Fs <- 0
while (criterioParada) {
  numeroIteraciones=numeroIteraciones+1
  # FUERA DE CONTROL EN SS
  idxFueraDeControlS_s <-as.integer(which(c(ss>limiteControlSuperiorSbar)==TRUE)) # Indice fuera de control superior para grafica s
  idxFueraDeControlI_s <-as.integer(which(c(ss<limiteControlInferiorSbar)==TRUE)) # Indice fuera de control inferior para grafica s
  idx_s <- c(idxFueraDeControlS_s,idxFueraDeControlI_s) # indices fuera de control para s
  nombresMuestra_s <- names(datosCombinados[,idx_s]) # Muestras fuera de control en grafico de s
  nReales_s <- sum(nombresMuestra_s %in% "C1") # Muestras que disparan señal cuando realmente provienen de una causa asignable
  nFalsos_s <- sum(nombresMuestra_s %in% c(nombresContaminados,nombresLimpios)) # Muestras que disparan señal cuando realmente provienen de una en control
  Ts <- Ts+nReales_s
  Fs <- Fs+nFalsos_s
  #---------------------------------------------#
  # SI NO HAY MUESTRAS QUE SE SALGAN DE CONTROL #
  # EN EL GRAFICO DE S, SE PROCEDE A VERIFICAR  #
  # EL GRAFICO DE X                             #
  #---------------------------------------------#
  if(length(idx_s)==0){
    idxFueraDeControlS_x <-as.integer(which(c(xbars>limiteControlSuperiorXbar)==TRUE)) # Indice fuera de control superior para grafica x
    idxFueraDeControlI_x <-as.integer(which(c(xbars<limiteControlInferiorXbar)==TRUE)) # Indice fuera de control inferior para grafica x
    idx_x <- c(idxFueraDeControlS_x,idxFueraDeControlI_x) # indices fuera de control para x
    nombresMuestra_x <- names(datosCombinados[,idx_x]) # Muestras fuera de control en grafico de x
    nReales_x <- sum(nombresMuestra_x %in% "C1") # Muestras que disparan señal cuando realmente provienen de una causa asignable
    nFalsos_x <- sum(nombresMuestra_x %in% c(nombresContaminados,nombresLimpios)) # Muestras que disparan señal cuando realmente provienen de una en control
    Ts <- Ts+nReales_x
    Fs <- Fs+nFalsos_x
    
    if(length(idx_x)==0){break}
    c(datosCombinados,xbars,ss,sbar,x_barbar,
      limiteControlInferiorXbar,limiteControlSuperiorXbar,
      limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                            idx=idx_x,
                                                                            n=n,
                                                                            C4=C4)
  }else{
    c(datosCombinados,xbars,ss,sbar,x_barbar,
      limiteControlInferiorXbar,limiteControlSuperiorXbar,
      limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                            idx=idx_s,
                                                                            n=n,
                                                                            C4=C4)
   
    
    }
  }
return(list("ANI"=numeroIteraciones,"T"=(Ts/(m*c))*100,"F"=(Fs/(m-m*c))*100))
}
# Replicas
v_F1_per <- Vectorize(F1_per) # Se vectoriza la funcion

set.seed(12)
Ls <- seq(1,5,0.1)
simulacion_1 <-  replicate(100,expr = v_F1_per(L=Ls,m=25,n=5,mu1=1,c=0.04))
simulacion_2 <-  replicate(100,expr = v_F1_per(L=Ls,m=25,n=5,mu1=2,c=0.04))
simulacion_3 <-  replicate(100,expr = v_F1_per(L=Ls,m=25,n=5,mu1=3,c=0.04))

ANI_1 <- apply(simulacion_1,2,function(x){return(unlist(x["ANI",]))})
ANI_2 <- apply(simulacion_2,2,function(x){return(unlist(x["ANI",]))})
ANI_3 <- apply(simulacion_3,2,function(x){return(unlist(x["ANI",]))})

ANIs_1 <- apply(ANI_1, 2, mean)
ANIs_2 <- apply(ANI_2, 2, mean)
ANIs_3 <- apply(ANI_3, 2, mean)
ANIs <- c(ANIs_1,ANIs_2,ANIs_3)

T_1 <- apply(simulacion_1,2,function(x){return(unlist(x["T",]))})
T_2 <- apply(simulacion_2,2,function(x){return(unlist(x["T",]))})
T_3 <- apply(simulacion_3,2,function(x){return(unlist(x["T",]))})

Ts_1 <- apply(T_1, 2, mean)
Ts_2 <- apply(T_2, 2, mean)
Ts_3 <- apply(T_3, 2, mean)
Ts <- c(Ts_1,Ts_2,Ts_3)

F_1 <- apply(simulacion_1,2,function(x){return(unlist(x["F",]))})
F_2 <- apply(simulacion_2,2,function(x){return(unlist(x["F",]))})
F_3 <- apply(simulacion_3,2,function(x){return(unlist(x["F",]))})

Fs_1 <- apply(F_1, 2, mean)
Fs_2 <- apply(F_2, 2, mean)
Fs_3 <- apply(F_3, 2, mean)
Fs <- c(Fs_1,Fs_2,Fs_3)

datos_plot <- data.frame("L"=rep(Ls,3),"ANI"=ANIs,"TAP"=Ts,"FAP"=Fs,"Shift"=rep(c("Shift:1","Shift:2","Shift:3"),each=length(Ls)))

library(plotly)
plot_ly(data=datos_plot,x=~L,y=~ANI,  color = ~Shift) %>% 
  add_lines()

plot_ly(data=datos_plot,x=~L,y=~FAP,  color = ~Shift) %>% 
  add_lines()

plot_ly(data=datos_plot,x=~L,y=~TAP,  color = ~Shift) %>% 
  add_lines()

#-----------------------
#-----------------------
#-----------------------
set.seed(12)
simulacion <- replicate(100,expr = v_F1_per(L=Ls,m=25,n=5,mu1=1,c=0.04))
ANI <- apply(simulacion,2,function(x){return(unlist(x["ANI",]))})
ANIss <- apply(ANI, 2, mean)
Ts <- apply(simulacion,2,function(x){return(unlist(x["T",]))})
Tss <- apply(Ts, 2, mean)
Fs <- apply(simulacion,2,function(x){return(unlist(x["F",]))})
Fss <- apply(Fs, 2, mean)
datos <- data.frame("L"=Ls,"ANI"=ANIss,"TAP"=Tss,"FAP"=Fss)

plot_ly(data=datos,x=~L,y=~ANI) %>% 
  add_lines()

plot_ly(data=datos,x=~L,y=~FAP) %>% 
  add_lines()

plot_ly(data=datos,x=~L,y=~TAP) %>% 
  add_lines()
