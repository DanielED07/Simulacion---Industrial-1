#-------------#
# MAS NUCLEOS # 
#-------------#
library(parallel)
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl,library(zeallot))
clusterSetRNGStream(cl)

set.seed(12)
#--------------#
# SIMULACIONES #
#--------------#
simulacion_1 <- parSapply(cl,Ls, function(i) {
  #--------------------------------------------------#
  # FUNCION PARA CALCULAR LIMITE DE CONTROL INFERIOR #
  #--------------------------------------------------#
  LCI <- function(cl,sigma,L){
    return(cl-L*sigma)
  }
  #--------------------------------------------------#
  # FUNCION PARA CALCULAR LIMITE DE CONTROL SUPERIOR #
  #--------------------------------------------------#
  LCS <- function(cl,sigma,L){
    return(cl+L*sigma)
  }
  #--------------------------#
  # FUNCION PARA CALCULAR c4 #
  #--------------------------#
  c4 <- function(n){
    return(sqrt(2/(n-1))*(gamma(n/2)/gamma((n-1)/2)))
  } 
  #---------------------------------------------------------------#
  # LAS SIGUIENTE FUNCION AYUDA A RECALCULAR LIMITES Y LA MUESTRA #
  #---------------------------------------------------------------#
  data_recalc <- function(datos,idx,n,C4,L){
    
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
    # set.seed(12)
    # L=3;m=25;n=5;mu1=3;c=0.04
    
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
    sigma_s <-(sbar/C4)*sqrt(1-C4^2)
    sigma_x <- sbar/(sqrt(n)*C4)
    
    limiteControlInferiorXbar <- LCI(cl=x_barbar,sigma = sigma_x,L=L)
    limiteControlSuperiorXbar <- LCS(cl=x_barbar,sigma = sigma_x,L=L)
    limiteControlInferiorSbar <- LCI(cl=sbar,sigma = sigma_s,L=L)
    limiteControlSuperiorSbar <- LCS(cl=sbar,sigma = sigma_s,L=L)
    
    # par(mfrow=c(1,2))
    # 
    # plot(x=1:ncol(datosCombinados),y=ss,type="line",ylim=c(0,2))
    # abline(h=c(limiteControlInferiorSbar,limiteControlSuperiorSbar),col="red")
    # abline(h=sbar,col="blue")
    # 
    # plot(x=1:ncol(datosCombinados),y=xbars,type="line",ylim=c(-3,3))
    # abline(h=c(limiteControlInferiorXbar,limiteControlSuperiorXbar),col="red")
    # abline(h=x_barbar,col="blue")
    
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
      nombresMuestra_s <- names(datosCombinados)[idx_s] # Muestras fuera de control en grafico de s
      nReales_s <- sum(nombresMuestra_s %in% nombresContaminados) # Muestras que disparan señal cuando realmente provienen de una causa asignable
      nFalsos_s <- sum(nombresMuestra_s %in% nombresLimpios) # Muestras que disparan señal cuando realmente provienen de una en control
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
        nombresMuestra_x <- names(datosCombinados)[idx_x] # Muestras fuera de control en grafico de x
        nReales_x <- sum(nombresMuestra_x %in% nombresContaminados) # Muestras que disparan señal cuando realmente provienen de una causa asignable
        nFalsos_x <- sum(nombresMuestra_x %in% nombresLimpios) # Muestras que disparan señal cuando realmente provienen de una en control
        Ts <- Ts+nReales_x
        Fs <- Fs+nFalsos_x
        
        if(length(idx_x)==0){break}
        c(datosCombinados,xbars,ss,sbar,x_barbar,
          limiteControlInferiorXbar,limiteControlSuperiorXbar,
          limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                                idx=idx_x,
                                                                                n=n,
                                                                                C4=C4,L=L)
      }else{
        c(datosCombinados,xbars,ss,sbar,x_barbar,
          limiteControlInferiorXbar,limiteControlSuperiorXbar,
          limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                                idx=idx_s,
                                                                                n=n,
                                                                                C4=C4,L=L)
        
        
      }
    }
    return(list("ANI"=numeroIteraciones,"T"=(Ts/(m*c))*100,"F"=(Fs/(m-m*c))*100,"Xbar"=(x_barbar-0)^2,"S_c4"=((sbar/C4)-1)^2))
  }
  # Se vectoriza la anterior función
  v_F1_per <- Vectorize(F1_per) # Se vectoriza la funcion
  
  #------------#
  # SIMULACION #
  #------------#
  N_sim <- 100000 # Numero de replicas
  replicate(N_sim,v_F1_per(L=i,m=25,n=5,mu1=1,c=0.04))
} )
#---------------------------------------------------------------------#
#---------------------------------------------------------------------#
#---------------------------------------------------------------------#
simulacion_2 <- parSapply(cl,Ls, function(i) {
  #--------------------------------------------------#
  # FUNCION PARA CALCULAR LIMITE DE CONTROL INFERIOR #
  #--------------------------------------------------#
  LCI <- function(cl,sigma,L){
    return(cl-L*sigma)
  }
  #--------------------------------------------------#
  # FUNCION PARA CALCULAR LIMITE DE CONTROL SUPERIOR #
  #--------------------------------------------------#
  LCS <- function(cl,sigma,L){
    return(cl+L*sigma)
  }
  #--------------------------#
  # FUNCION PARA CALCULAR c4 #
  #--------------------------#
  c4 <- function(n){
    return(sqrt(2/(n-1))*(gamma(n/2)/gamma((n-1)/2)))
  } 
  #---------------------------------------------------------------#
  # LAS SIGUIENTE FUNCION AYUDA A RECALCULAR LIMITES Y LA MUESTRA #
  #---------------------------------------------------------------#
  data_recalc <- function(datos,idx,n,C4,L){
    
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
    # set.seed(12)
    # L=3;m=25;n=5;mu1=3;c=0.04
    
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
    sigma_s <-(sbar/C4)*sqrt(1-C4^2)
    sigma_x <- sbar/(sqrt(n)*C4)
    
    limiteControlInferiorXbar <- LCI(cl=x_barbar,sigma = sigma_x,L=L)
    limiteControlSuperiorXbar <- LCS(cl=x_barbar,sigma = sigma_x,L=L)
    limiteControlInferiorSbar <- LCI(cl=sbar,sigma = sigma_s,L=L)
    limiteControlSuperiorSbar <- LCS(cl=sbar,sigma = sigma_s,L=L)
    
    # par(mfrow=c(1,2))
    # 
    # plot(x=1:ncol(datosCombinados),y=ss,type="line",ylim=c(0,2))
    # abline(h=c(limiteControlInferiorSbar,limiteControlSuperiorSbar),col="red")
    # abline(h=sbar,col="blue")
    # 
    # plot(x=1:ncol(datosCombinados),y=xbars,type="line",ylim=c(-3,3))
    # abline(h=c(limiteControlInferiorXbar,limiteControlSuperiorXbar),col="red")
    # abline(h=x_barbar,col="blue")
    
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
      nombresMuestra_s <- names(datosCombinados)[idx_s] # Muestras fuera de control en grafico de s
      nReales_s <- sum(nombresMuestra_s %in% nombresContaminados) # Muestras que disparan señal cuando realmente provienen de una causa asignable
      nFalsos_s <- sum(nombresMuestra_s %in% nombresLimpios) # Muestras que disparan señal cuando realmente provienen de una en control
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
        nombresMuestra_x <- names(datosCombinados)[idx_x] # Muestras fuera de control en grafico de x
        nReales_x <- sum(nombresMuestra_x %in% nombresContaminados) # Muestras que disparan señal cuando realmente provienen de una causa asignable
        nFalsos_x <- sum(nombresMuestra_x %in% nombresLimpios) # Muestras que disparan señal cuando realmente provienen de una en control
        Ts <- Ts+nReales_x
        Fs <- Fs+nFalsos_x
        
        if(length(idx_x)==0){break}
        c(datosCombinados,xbars,ss,sbar,x_barbar,
          limiteControlInferiorXbar,limiteControlSuperiorXbar,
          limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                                idx=idx_x,
                                                                                n=n,
                                                                                C4=C4,L=L)
      }else{
        c(datosCombinados,xbars,ss,sbar,x_barbar,
          limiteControlInferiorXbar,limiteControlSuperiorXbar,
          limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                                idx=idx_s,
                                                                                n=n,
                                                                                C4=C4,L=L)
        
        
      }
    }
    return(list("ANI"=numeroIteraciones,"T"=(Ts/(m*c))*100,"F"=(Fs/(m-m*c))*100,"Xbar"=(x_barbar-0)^2,"S_c4"=((sbar/C4)-1)^2))
  }
  # Se vectoriza la anterior función
  v_F1_per <- Vectorize(F1_per) # Se vectoriza la funcion
  
  #------------#
  # SIMULACION #
  #------------#
  N_sim <- 100000 # Numero de replicas
  replicate(N_sim,v_F1_per(L=i,m=25,n=5,mu1=2,c=0.04))
} )
#---------------------------------------------------------------------#
#---------------------------------------------------------------------#
#---------------------------------------------------------------------#
simulacion_3 <- parSapply(cl,Ls, function(i) {
  #--------------------------------------------------#
  # FUNCION PARA CALCULAR LIMITE DE CONTROL INFERIOR #
  #--------------------------------------------------#
  LCI <- function(cl,sigma,L){
    return(cl-L*sigma)
  }
  #--------------------------------------------------#
  # FUNCION PARA CALCULAR LIMITE DE CONTROL SUPERIOR #
  #--------------------------------------------------#
  LCS <- function(cl,sigma,L){
    return(cl+L*sigma)
  }
  #--------------------------#
  # FUNCION PARA CALCULAR c4 #
  #--------------------------#
  c4 <- function(n){
    return(sqrt(2/(n-1))*(gamma(n/2)/gamma((n-1)/2)))
  } 
  #---------------------------------------------------------------#
  # LAS SIGUIENTE FUNCION AYUDA A RECALCULAR LIMITES Y LA MUESTRA #
  #---------------------------------------------------------------#
  data_recalc <- function(datos,idx,n,C4,L){
    
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
    # set.seed(12)
    # L=3;m=25;n=5;mu1=3;c=0.04
    
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
    sigma_s <-(sbar/C4)*sqrt(1-C4^2)
    sigma_x <- sbar/(sqrt(n)*C4)
    
    limiteControlInferiorXbar <- LCI(cl=x_barbar,sigma = sigma_x,L=L)
    limiteControlSuperiorXbar <- LCS(cl=x_barbar,sigma = sigma_x,L=L)
    limiteControlInferiorSbar <- LCI(cl=sbar,sigma = sigma_s,L=L)
    limiteControlSuperiorSbar <- LCS(cl=sbar,sigma = sigma_s,L=L)
    
    # par(mfrow=c(1,2))
    # 
    # plot(x=1:ncol(datosCombinados),y=ss,type="line",ylim=c(0,2))
    # abline(h=c(limiteControlInferiorSbar,limiteControlSuperiorSbar),col="red")
    # abline(h=sbar,col="blue")
    # 
    # plot(x=1:ncol(datosCombinados),y=xbars,type="line",ylim=c(-3,3))
    # abline(h=c(limiteControlInferiorXbar,limiteControlSuperiorXbar),col="red")
    # abline(h=x_barbar,col="blue")
    
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
      nombresMuestra_s <- names(datosCombinados)[idx_s] # Muestras fuera de control en grafico de s
      nReales_s <- sum(nombresMuestra_s %in% nombresContaminados) # Muestras que disparan señal cuando realmente provienen de una causa asignable
      nFalsos_s <- sum(nombresMuestra_s %in% nombresLimpios) # Muestras que disparan señal cuando realmente provienen de una en control
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
        nombresMuestra_x <- names(datosCombinados)[idx_x] # Muestras fuera de control en grafico de x
        nReales_x <- sum(nombresMuestra_x %in% nombresContaminados) # Muestras que disparan señal cuando realmente provienen de una causa asignable
        nFalsos_x <- sum(nombresMuestra_x %in% nombresLimpios) # Muestras que disparan señal cuando realmente provienen de una en control
        Ts <- Ts+nReales_x
        Fs <- Fs+nFalsos_x
        
        if(length(idx_x)==0){break}
        c(datosCombinados,xbars,ss,sbar,x_barbar,
          limiteControlInferiorXbar,limiteControlSuperiorXbar,
          limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                                idx=idx_x,
                                                                                n=n,
                                                                                C4=C4,L=L)
      }else{
        c(datosCombinados,xbars,ss,sbar,x_barbar,
          limiteControlInferiorXbar,limiteControlSuperiorXbar,
          limiteControlInferiorSbar,limiteControlSuperiorSbar) %<-% data_recalc(datos=datosCombinados,
                                                                                idx=idx_s,
                                                                                n=n,
                                                                                C4=C4,L=L)
        
        
      }
    }
    return(list("ANI"=numeroIteraciones,"T"=(Ts/(m*c))*100,"F"=(Fs/(m-m*c))*100,"Xbar"=(x_barbar-0)^2,"S_c4"=((sbar/C4)-1)^2))
  }
  # Se vectoriza la anterior función
  v_F1_per <- Vectorize(F1_per) # Se vectoriza la funcion
  
  #------------#
  # SIMULACION #
  #------------#
  N_sim <- 100000 # Numero de replicas
  replicate(N_sim,v_F1_per(L=i,m=25,n=5,mu1=3,c=0.04))
} )
stopCluster(cl)

simulaciones <- grep("simulacion",ls(),value = T)
#----------------------#
#Organizemos los datos #
#----------------------#

datos_plot <- data.frame("L"=rep(Ls,3),"ANI"=0,"TAP"=0,"FAP"=0,"MSE_X"=0,"MSE_Sc4"=0,"Shift"=rep(c("Shift:1","Shift:2","Shift:3"),each=length(Ls)))
a <- c(1:41) # Indices para la simulacion 1 (Shift : 1)
b <- c(42:82) # Indices para la simulacion 2 (Shift : 2)
c <- c(83:123) # Indices para la simulacion 3 (Shift : 3)
for (j in 1:3) {

  data <-  get(simulaciones[j])
  
  if(j==1){indxs <- a}else{if(j==2){indxs <-b}else{indxs <-c}}
  
  for (i in 1:41) {
    N_sim <- 100
    metricas <- data.frame(matrix(0,ncol=5,nrow = N_sim))
    names(metricas) <- c("ANI","TAP","FAP","MSE_X","MSE_Sc4")
    
    L <-unlist(data[,i])  # L=1
    
    idx_ANI <- seq(1,5*N_sim,5) # Indices para extraer los ANIs
    idx_TAP <- seq(2,5*N_sim,5) # Indices para extraer los TAPs
    idx_FAP <- seq(3,5*N_sim,5) # Indices para extraer los FAPs
    idx_MSE_X <- seq(4,5*N_sim,5) # Indices para extraer los MSE_Xs
    idx_MSE_Sc4 <- seq(5,5*N_sim,5) # Indices para extraer los MSE_Sc4s
    
    metricas$ANI <- L[idx_ANI]
    metricas$TAP <- L[idx_TAP]
    metricas$FAP <- L[idx_FAP]
    metricas$MSE_X <- L[idx_MSE_X]
    metricas$MSE_Sc4 <- L[idx_MSE_Sc4]
    
    estimaciones_metricas <- apply(metricas, 2, mean)
    datos_plot[indxs[i],c("ANI","TAP","FAP","MSE_X","MSE_Sc4")] <- estimaciones_metricas
  }
}

#save(datos_plot,file = "datos_plot.RData")





