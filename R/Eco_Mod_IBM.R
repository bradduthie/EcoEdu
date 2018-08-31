# ----------------------------------------------------------
# XXX HIGHLIGHT VARIABLE INPUT VALUES BELOW THEN HIT:
# XXX CONTROL + R (FOR WINDOWS)
# XXX CONTROL + ENTER (FOR MAC)
# XXX CONTROL + F8 (FOR LINUX)
# ----------------------------------------------------------
# VARIABLE INPUT
# ----------------------------------------------------------
time_max  <- 100;    # Maximum time steps
N_initial <- 0.4;    # Initial density of prey
P_initial <- 0.05;   # Initial density of predators
lambda_N  <- 1.2;    # Prey growth rate
lambda_P  <- 0.3;    # Predator decline rate
attack    <- 0.5;    # Attack rate of predators
K         <- 0.5;    # Prey density carrying capacity
b         <- 0.2;    # Predator births per attack
mY        <- 1;      # Prey movement
mP        <- 1;      # Predator movement
slowdown  <- 0.2;    # Slows the simulation time down
# ----------------------------------------------------------

# ----------------------------------------------------------
# BRING UP A NEW PLOTTING WINDOW
# ----------------------------------------------------------
par(mar=c(5,5,1,1),lwd=2);
plot(x=0,y=0,type="n",xlim=c(0,time_max),ylim=c(0,1.4),yaxt="n",
     xlab="Time step",ylab="Species density",cex.lab=1.8,cex.axis=1.8);
axis(side=2,at=c(0,0.2,0.4,0.6,0.8),cex.axis=1.8);
legend(x=0,y=1,legend=c("Prey","Predator"),col=c("black","red"),cex=1,pch=15,horiz=TRUE);
polygon(x=c(0,time_max-2,time_max-2,0),y=c(1,1,1.41,1.41), border="black",col="tan",lwd=3);
# ----------------------------------------------------------

# ----------------------------------------------------------
# MODEL INITIALISATION
# ----------------------------------------------------------
Land      <- 100 * 100;                 # Total size of the landscape
Kmax      <- round(K * Land);           # Total carrying capacity number
N_initial <- round(N_initial * Land);   # Total number of initial prey
P_initial <- round(P_initial * Land);   # Total number of initial predators
N         <- rep(x=0, times=time_max);  # Prey abundance over time
P         <- rep(x=0, times=time_max);  # Predator abundance over time
N[1]      <- N_initial;                 # Set initial prey abundance
P[1]      <- P_initial;                 # Set initial predator abundance
NPOS      <- NULL;                      # Will be prey positions on landscape
PPOS      <- NULL;                      # Will be predator positionson landscape
Np1  <- round(runif(n=N_initial,min=0,max=100),digits=0); # Random x position
Np2  <- round(runif(n=N_initial,min=0,max=100),digits=0); # Random y position
Np   <- cbind(Np1,Np2); # Prey x and y locations on the landscape
Pext <- 1; # Are predators extinct? Assume so at first
if(P_initial > 0){ # If predators are not initially extinct (start with some)
    Pp1  <- round(runif(n=P_initial,min=0,max=100),digits=0);
    Pp2  <- round(runif(n=P_initial,min=0,max=100),digits=0);
    Pp   <- cbind(Pp1,Pp2); # Predator x and y locations on the landscape
    Pext <- 0; # And clarify that they are not extinct
}
timestep  <- 1; # Start with the first time step
# ----------------------------------------------------------
# The big loop (simulating prey and predator interactions over time
# ----------------------------------------------------------
while(timestep < time_max){ # As long as the time step is below the maximum
    if(N[timestep] == 0){   # First, are there any prey left?
        break;  # If not, break out of this while loop and end the thing!
    } # Nm1 and Nm2 below are distance moved in x and y on landscape
    Nm1  <- round(runif(n=dim(Np)[1],min=-mY,max=mY),digits=0);
    Nm2  <- round(runif(n=dim(Np)[1],min=-mY,max=mY),digits=0);
    No   <- cbind(Nm1,Nm2); # Now have a 2 column table of each prey's movement
    Np   <- Np + No; # Add the prey's movement to their current positions
    if(P_initial > 0 & Pext == 0){ # If started with predators, & not extinct
        if(is.matrix(Pp)==TRUE){ # Small R hangup -- make sure is a matrix
            Po <- matrix(data=round(runif(n=length(Pp),min=-mP,max=mP)),ncol=2);
            Pp <- Pp + Po; # Makes the predators move too
            newpred <- NULL; # Note, if not matrix, above calc would fail in R
        } # Now we have moved the predators, if there were any
    }   # - This for loop below creates a `torus' landscape for the prey
    # So if prey go off one edge of the landscape, they end up on the other side
    for(i in 1:dim(Np)[1]){
        if(Np[i,1] > 100){ # Note: Check if over the landscape right edge
            Np[i,1] <- Np[i,1] - 100; # If so, bring back to left side
        }
        if(Np[i,2] > 100){ # Do the same for top-bottom, and other edges
            Np[i,2] <- Np[i,2] - 100;
        }
        if(Np[i,1] < 0){
            Np[i,1] <- Np[i,1] + 100;
        }
        if(Np[i,2] < 0){
            Np[i,2] <- Np[i,2] + 100;
        }
    } # The for loop keeping prey on the landscape is now complete
    if(P_initial > 0 & Pext == 0){ # If we have predators that aren't extinct
        if(is.matrix(Pp)==TRUE){ # Need this for R to work; are ways around it
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
            } # For loop keeping predators on the landscape is now complete
            # Things start to get tricky: Below check if predators encounter prey
            for(i in 1:dim(Pp)[1]){ # For each predator, sequentially
                # Note, defining an encounter as within TWO spaces of a prey, can change.
                enc  <- which(abs(Pp[i,1]-Np[,1])<=2 & abs(Pp[i,2]-Np[,2])<=2); 
                if(sum(enc) > 0){ # If at least 1 prey is 2 spaces away from i
                    odz  <- runif(n = 2*length(enc), min = 0, max = 1); # Get random #s
                    kill <- matrix(data = odz, ncol = 2); # Put above in 2 col matrix
                    kill <- (0.5*attack) > kill; # See if random < attack rate (1 if true)
                    dead <- which(apply(X=kill,MARGIN=1,FUN=sum)>0); # Which are 1 (and dead)?
                    if(sum(dead) > 0){ # If at least 1 prey died from the attack
                        Np[enc[dead],] <- c(-200,-200); # Change its position to 2 negatives
                        bodz           <- runif(n = 2*length(dead), min = 0, max = 1); #Birth?
                        bodz           <- matrix(data = bodz, ncol = 2); # Rand #s matrix
                        born           <- 0.5*b > bodz; # Returns 1 if rand # lower than 0.5*b
                        birth          <- sum(born); # Sum the number of 1s (total births)
                        if(birth > 0){ # If there was at least one predator birth
                            tbrn    <- birth; # use the count for the total here
                            Pp1     <- round(runif(n=tbrn,min=0,max=100),digits=0);
                            Pp2     <- round(runif(n=tbrn,min=0,max=100),digits=0);
                            newpred <- rbind(newpred,cbind(Pp1,Pp2)); # Newly born predators
                        } # Above 3 lines randomly place new predators on the landscape
                    }
               } # Now about to exit the loop of predator attacks for this generation
           } # Note that prey that have been eaten were given x & y locations of -200
           if(sum(Np[,1]==-200) > 0){ # If there were any prey eaten, remove below
               dead <- which(Np[,1]==-200); # Identifies who was eaten
               Np   <- Np[-dead,]; # Removes them from the table of prey
           } # Now have only the uneaten prey in the table Np
           Pp      <- rbind(Pp,newpred); # Add new predators to the table of predators
           psurv1  <- rbinom(n=dim(Pp)[1], size=1, prob=(1 - 0.5*lambda_P)); # Predator death
           psurv2  <- rbinom(n=dim(Pp)[1], size=1, prob=(1 - 0.5*lambda_P)); # 2 chances to die
           psurv   <- sum(psurv1 == 1 & psurv2 == 1); # Sum of preds surviving both chances
           if(psurv > 1){ # If at least 2 preds survived (assume need 2 to avoid extinction)
               pkeep   <- which(psurv1 == 1 & psurv2 == 1); # Retain the survivors
               Pp      <- Pp[pkeep,]; # Define new Pp table as survivors from the old Pp table
           }else{ # If one or fewer predators survived this generation
               Pext    <- 1; # Then the predators have gone extinct
           }
       } # Coming out of the big if statement (line 78) that requires predators to exist
    } # Below just consider birth and death of prey independent of predators
    birthp <- runif(n=dim(Np)[1],min=0,max=1); # Some random numbers to see if prey repr
    births <- sum(birthp < (0.5*lambda_N)); # If rand # lower than growth, reproduction
    if(births > 0){ # If some prey give birth, give new positions and attach to Np below
        Np1    <- round(runif(n=births,min=0,max=100),digits=0);
        Np2    <- round(runif(n=births,min=0,max=100),digits=0);
        newNs  <- cbind(Np1,Np2);
        Np     <- rbind(Np,newNs);
    }
    if(is.null(dim(Np)[1])){ # Check to see if prey have gone extinct
        break; # If so, break out of the big while loop and end the simulation
    }
    extr   <- dim(Np)[1] - Kmax; # How far above carrying capacity is the prey population?
    if(extr > 0){ # If it's over, increase all prey prob of death accordingly
        mort <- extr / (Kmax + extr);
        surv <- rbinom(n=dim(Np)[1], size=1, prob=(1-mort));
        keep <- which(surv==1);
        Np   <- Np[keep,];
    } # Now have retained survivors of density regulation
    N[timestep + 1] <- dim(Np)[1]; # Prey abundance in the whole population
    P[timestep + 1] <- dim(Pp)[1]; # Predator abundance in the whole population
    # ------- Below plots prey, predators, and their graphs.
    polygon(x=c(0,time_max-2,time_max-2,0),y=c(1,1,1.41,1.41), border="black",col="tan",lwd=3);
    xpos <- 1 + time_max * (Np[,1] / (time_max + 5));
    ypos <- 0.39 * (Np[,2] / time_max) + 1.01;
    points(x=xpos,ypos,pch=20,col="black");
    xpos <- 1 + time_max * (Pp[,1] / (time_max + 5));
    ypos <- 0.39 * (Pp[,2] / time_max) + 1.01;
    points(x=xpos,ypos,pch=20,col="red");
    points(x=1:timestep,y=N[1:timestep]/Land,type="l",col="black",lwd=2);
    points(x=1:timestep,y=P[1:timestep]/Land,type="l",col="red",lwd=2);
    #-------- End the plot, now just need to move onto the next time step
    timestep <- timestep + 1;
    Sys.sleep(slowdown); # Just slows the simulation down so you can see it.
} # Finally, close out of the big while loop when timstep > time_max
# ----------------------------------------------------------
# Note: After running through the above loop, you now have a record of prey abundances (N) and predator abundances (P) over time.
# ----------------------------------------------------------





