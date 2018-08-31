

library(shiny)

# Define server logic for slider examples
shinyServer(function(input, output, session) {
  
  # initialize reactive values
  ts <- reactiveValues(ct=1,ds=NULL)


    Ds <- eventReactive(input$go,{
        time_max  <- input$time         # Maximum time steps
        N_initial <- input$N_initial    # Initial number of prey
        P_initial <- input$P_initial    # Initial number of predators
        lambda_N  <- input$lambda_N     # Prey growth rate
        lambda_P  <- input$lambda_P     # Predator decline rate
        attack    <- input$attack       # Attack rate of predators
        K         <- input$K            # Prey carrying capacity
        b         <- input$b            # Predator births per attack
        mY        <- input$dispY * 100  # Prey movement
        mP        <- input$dispP * 100  # Predator movement
        if(input$variable == "Numer"){   
            timestep  <- 1;
            N         <- rep(x=0, times=time_max);
            P         <- rep(x=0, times=time_max);
            N[1]      <- N_initial;
            P[1]      <- P_initial;
            while(timestep < time_max){
                N_growth        <- (lambda_N * N[timestep]) * (1 - (N[timestep] / K));
                N_attacked      <- attack * P[timestep] * N[timestep];
                N[timestep + 1] <- N[timestep] + N_growth - N_attacked;
                P_death         <- (lambda_P * P[timestep]);
                P_birth         <- b * attack * P[timestep] * N[timestep];
                P[timestep + 1] <- P[timestep] - P_death + P_birth;
                timestep        <- timestep + 1;
            }
            Abs <- cbind(N,P);
            FPO <- NULL;
            return(list(Abs,FPO));
        }else{
            timestep  <- 1;
            Kmax      <- round(K   * 1000);
            N_initial <- round(N_initial * 1000);
            P_initial <- round(P_initial * 1000);
            N         <- rep(x=0, times=time_max);
            P         <- rep(x=0, times=time_max);
            N[1]      <- N_initial;
            P[1]      <- P_initial;
            NPOS      <- NULL;
            PPOS      <- NULL;
            Np1  <- round(runif(n=N_initial,min=0,max=100),digits=0);
            Np2  <- round(runif(n=N_initial,min=0,max=100),digits=0);
            Np   <- cbind(Np1,Np2);
            Pext <- 1;
            if(P_initial > 0){
                Pp1  <- round(runif(n=P_initial,min=0,max=100),digits=0);
                Pp2  <- round(runif(n=P_initial,min=0,max=100),digits=0);
                Pp   <- cbind(Pp1,Pp2);
                Pext <- 0;
            }
            while(timestep < time_max){
                if(N[timestep] == 0){
                    break;
                }
                Nm1  <- round(runif(n=dim(Np)[1],min=-mY,max=mY),digits=0);
                Nm2  <- round(runif(n=dim(Np)[1],min=-mY,max=mY),digits=0);
                No   <- cbind(Nm1,Nm2);
                Np   <- Np + No;
                if(P_initial > 0 & Pext == 0){
                    if(is.matrix(Pp)==TRUE){
                        Po <- matrix(data=round(runif(n=length(Pp),min=-mP,max=mP)),ncol=2);
                        Pp <- Pp + Po;
                        newpred <- NULL;
                    }
                }
                for(i in 1:dim(Np)[1]){ # These create a torus landscape
                    if(Np[i,1] > 100){
                        Np[i,1] <- Np[i,1] - 100;
                    }
                    if(Np[i,2] > 100){
                        Np[i,2] <- Np[i,2] - 100;
                    }
                    if(Np[i,1] < 0){
                        Np[i,1] <- Np[i,1] + 100;
                    }
                    if(Np[i,2] < 0){
                        Np[i,2] <- Np[i,2] + 100;
                    }
                }
                if(P_initial > 0 & Pext == 0){
                    if(is.matrix(Pp)==TRUE){
                        for(i in 1:dim(Pp)[1]){ # These create a torus landscape
                            if(Pp[i,1] > 100){
                                Pp[i,1] <- Pp[i,1] - 100;
                            }
                            if(Pp[i,2] > 100){
                                Pp[i,2] <- Pp[i,2] - 100;
                            }
                            if(Pp[i,1] < 0){
                                Pp[i,1] <- Pp[i,1] + 100;
                            }
                            if(Pp[i,2] < 0){
                                Pp[i,2] <- Pp[i,2] + 100;
                            }
                        }
                        for(i in 1:dim(Pp)[1]){ # Explain why 1s and 2s work later.
                            for(j in 1:dim(Np)[1]){
                                if(abs(Pp[i,1]-Np[j,1])<=2 & abs(Pp[i,2]-Np[j,2])<=2){
                                    kill1  <- 0.5*attack > runif(n=1,min=0,max=1);
                                    kill2  <- 0.5*attack > runif(n=1,min=0,max=1);
                                    if(kill1 == TRUE | kill2 == TRUE){
                                        Np[j,] <- c(-200,-200);
                                        birth1   <- 0.5*b > runif(n=1,min=0,max=1);
                                        birth2   <- 0.5*b > runif(n=1,min=0,max=1);
                                        if(birth1==TRUE){ 
                                           Pp1     <- round(runif(n=1,min=0,max=100),digits=0);
                                           Pp2     <- round(runif(n=1,min=0,max=100),digits=0);
                                           newpred <- rbind(newpred,c(Pp1,Pp2));
                                        }
                                        if(birth2==TRUE){
                                           Pp1     <- round(runif(n=1,min=0,max=100),digits=0);
                                           Pp2     <- round(runif(n=1,min=0,max=100),digits=0);
                                           newpred <- rbind(newpred,c(Pp1,Pp2));
                                        }
                                    }
                                }
                            }
                        }
                        if(sum(Np[,1]==-200) > 0){ # Removes eaten prey
                            dead <- which(Np[,1]==-200);
                            Np   <- Np[-dead,];
                        }
                        Pp      <- rbind(Pp,newpred);
                        psurv1  <- rbinom(n=dim(Pp)[1], size=1, prob=(1 - 0.5*lambda_P));
                        psurv2  <- rbinom(n=dim(Pp)[1], size=1, prob=(1 - 0.5*lambda_P));
                        psurv   <- sum(psurv1 == 1 & psurv2 == 1);
                        if(psurv > 1){
                            pkeep   <- which(psurv1 == 1 & psurv2 == 1);
                            Pp      <- Pp[pkeep,];
                        }else{
                            Pext    <- 1;
                        }           
                    }
                }
                birthp <- runif(n=dim(Np)[1],min=0,max=1);
                births <- sum(birthp < (0.5*lambda_N));
                Np1    <- round(runif(n=births,min=0,max=100),digits=0);
                Np2    <- round(runif(n=births,min=0,max=100),digits=0);
                newNs  <- cbind(Np1,Np2);
                Np     <- rbind(Np,newNs);
                if(is.null(dim(Np)[1])){
                    break;
                }
                extr   <- dim(Np)[1] - Kmax;
                if(extr > 0){
                    mort <- extr / (Kmax + extr);
                    surv <- rbinom(n=dim(Np)[1], size=1, prob=(1-mort));
                    keep <- which(surv==1);
                    Np   <- Np[keep,];
                }
                NPOS            <- rbind(NPOS,cbind(rep(timestep,dim(Np)[1]),Np));
                N[timestep + 1] <- dim(Np)[1];
                if(P_initial > 0 & Pext == 0){
                    if(is.matrix(Pp) == TRUE){
                        PPOS            <- rbind(PPOS,cbind(rep(timestep,dim(Pp)[1]),Pp));
                        P[timestep + 1] <- dim(Pp)[1];
                    }else{
                        P[timestep + 1] <- 0;
                    }
                }
                timestep        <- timestep + 1;
            }
            Abs <- cbind(N,P);
            if(P_initial > 0){ 
                return(list(Abs,NPOS,PPOS));
            }else{
                return(list(Abs=Abs,NPOS=NPOS));
            }
        }


    })

    output$distPlot <- renderPlot({
        if(input$variable == "Numer"){   
            Ress <- Ds();
            Dens <- Ress[[1]];
            endd <- ts$ct;
            par(mar=c(5,5,1,1),lwd=2);
            plot(x=0,y=0,type="n",xlim=c(0,input$time),ylim=c(0,1.4),yaxt="n",
                xlab="Time step",ylab="Species density",cex.lab=1.8,cex.axis=1.8);
            axis(side=2,at=c(0,0.2,0.4,0.6,0.8),cex.axis=1.8);
            legend(x=0,y=1,legend=c("Prey","Predator"),col=c("black","red"),cex=1,
                   pch=15,horiz=TRUE);
            points(x=1:endd,y=Dens[1:endd,1],type="l",col="black");
            if(input$P_initial > 0){
                points(x=1:endd,y=Dens[1:endd,2],type="l",col="red");
            }
            polygon(x=c(0,input$time,input$time,0),y=c(1,1,1.4,1.4),
                    border="black",col="grey85");
            prd <- round(Dens[endd,1],digits=2);
            pdd <- round(Dens[endd,2],digits=2);
            if(input$P_initial > 0){
                text(x=input$time/2, y=1.3, cex=1.5,
                     labels=bquote(paste(.(prd) == N[.(ts$ct)] + (.(input$lambda_N)) 
                                     * N[.(ts$ct)] * (1-N[.(ts$ct)]/.(input$K)) 
                                     - (.(input$attack))*P[.(ts$ct)]*N[.(ts$ct)]  )));
                text(x=input$time/2, y=1.1, cex=1.5, col="red",
                     labels=bquote(paste(.(pdd) == P[.(ts$ct)] + .(input$attack) * 
                               (.(input$b)) * P[.(ts$ct)] * N[.(ts$ct)] -
                               (.(input$lambda_P)) * P[.(ts$ct)] )));
            }else{
                text(x=input$time/2, y=1.3, cex=1.5,
                     labels=bquote(paste(.(prd) == N[.(ts$ct)] + (.(input$lambda_N)) 
                               * N[.(ts$ct)] * (1-N[.(ts$ct)]/.(input$K)))));
            }
        }else{
            Ress <- Ds();
            Dens <- Ress[[1]] * 0.001;
            FPO  <- Ress[[2]];
            if(ts$ct == input$time){
                gen <- ts$ct - 1;
            }else{
                gen <- ts$ct;
            }
            FPO  <- FPO[FPO[,1]==gen,];
            if(length(Ress) > 2){
                PPO <- Ress[[3]];
                PPO <- PPO[PPO[,1]==gen,];
            }
            endd <- gen;
            par(mar=c(5,5,1,1),lwd=2);
            plot(x=0,y=0,type="n",xlim=c(0,input$time),ylim=c(0,1.4),yaxt="n",
                xlab="Time step",ylab="Species density",cex.lab=1.8,cex.axis=1.8);
            axis(side=2,at=c(0,0.2,0.4,0.6,0.8),cex.axis=1.8);
            legend(x=0,y=1,legend=c("Prey","Predator"),col=c("black","red"),cex=1,
                   pch=15,horiz=TRUE);
            points(x=1:endd,y=Dens[1:endd,1],type="l",col="black");
            if(input$P_initial > 0){
                points(x=1:endd,y=Dens[1:endd,2],type="l",col="red");
            }
            polygon(x=c(0,input$time,input$time,0),y=c(1,1,1.4,1.4),
                    border="black",col="tan",lwd=3);
            if(sum(FPO) > 0){
                xpos <- 1 + input$time * (FPO[,2] / (input$time + 5));
                ypos <- 0.39 * (FPO[,3] / input$time) + 1.01;
                points(x=xpos,ypos,pch=20);
            }
            if(length(Ress) > 2){
                xpos <- 1 + input$time * (PPO[,2] / (input$time + 5));
                ypos <- 0.39 * (PPO[,3] / input$time) + 1.01;
                points(x=xpos,ypos,pch=20,col="red");
            }
        }
    })


  observe({
    isolate({
    if(ts$ct < input$time){
        ts$ct <- ts$ct + 1;
    }
    }) 
   if (((isolate(ts$ct) < input$time)) & ((input$go > 0) | (input$finish > 0))){
      invalidateLater(400, session)
    }
  })

   observeEvent(input$reset, {
       ts$ct <- 1;
       seed  <- sample(x=1:100000,size=1);
       set.seed(seed);
   })

   observeEvent(input$finish, {
       ts$ct <- input$time;
   })

})




