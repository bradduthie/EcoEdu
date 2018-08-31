Evolve <- function(A=0.5, a=0.5, AA=0, Aa=0, aa=0, generations=150,PopSize=100,wAA=1,wAa=1,waa=1,
				M.AA=0,M.Aa=0,M.aa=0,Graph.Alleles=TRUE,Graph.Genotypes=FALSE,Show.Data=FALSE,X11=TRUE){
	if(X11==1){
		X11(width=6,height=6) 	#First call up a new window in which to plot things
		}
	#Below is an error message that stops the simulation if allele frequencies don't sume to one
	if(A+a!=1){ #If the allele frequencies don't add up to one. The program barks at you
		plot(x=0,y=0,xlim=c(0,generations),ylim=c(0,1.2),type="n",xlab="Generation",ylab="Frequency")
			legend(x=0,y=0.75,"ALLELE FREQUENCIES NEED TO ADD UP TO ONE",cex=1,bty="n")
		stop("ALLELE FREQUENCIES NEED TO ADD UP TO ONE")
		} #This puts a warning in the command line and plot that allele frequencies need to sum to 1.
	#Below is an error message that stops the simulation if relative survival is greater than one
	if(wAA > 1 || wAa > 1 || waa > 1){ #If a relative survival value is greater than one
		plot(x=0,y=0,xlim=c(0,generations),ylim=c(0,1.2),type="n",xlab="Generation",ylab="Frequency")
			legend(x=0,y=0.75,"A RELATIVE SURVIVAL VALUE CAN'T BE GREATER THAN ONE",cex=1,bty="n")
		stop("A RELATIVE SURVIVAL VALUE CAN'T BE GREATER THAN ONE")
		} #This puts a warning on the command line and in the plot that relative survival isn't > 1
	#Below stops the simulation if the population size is too large.
	#Over 1,000,000 is unnecessary and will slow down the simulation.
	if(PopSize>1000000) {
		plot(x=0,y=0,xlim=c(0,generations),ylim=c(0,1.2),type="n",xlab="Generation",ylab="Frequency")
			legend(x=0,y=0.75,"A POPULATION SIZE OVER ONE MILLION WILL TAKE WAY TOO LONG!",cex=1,bty="n")
		stop("A POPULATION SIZE OVER ONE MILLION WILL TAKE WAY TOO LONG!")
		}	#This puts a warning on the command line and in the plot that over 1000000 is too high
	if(AA==0 && Aa==0 && aa==0){ #If individual counts aren't set
		AA <- A*A*PopSize #Then figure out the number of indivuduals
		Aa <- 2*A*a*PopSize #by using HW times the population size
		aa <- a*a*PopSize
		}else
			{ #Else, if individual conts are set, then find the frequencies.
			A <- ((2*AA)+Aa)/(2*sum(AA,Aa,aa))
			a <- ((2*aa)+Aa)/(2*sum(AA,Aa,aa))
			}
	if((AA+Aa+aa)!=PopSize){ #If the total number of individuals doesn't equal the PopSize input
		Popsize <- AA+Aa+aa #Override the PopSize input. The real PopSize is the sum of individuals
		}	 
	FreqA <- Freqa <- NULL #Set a null for later allele frequencies
	FreqAA <- FreqAa <- Freqaa <- NULL #Set a null for genotype frequencies
	#Stick down a plot that will be used later.
	plot(x=0,y=0,xlim=c(0,generations),ylim=c(0,1.2),type="n",xlab="Generation",ylab="Frequency")
	if(Graph.Alleles==1 && Graph.Genotypes==0){ #If we only are graphing allele frequencies
		legend(x=0,y=1.2,fill=c("red","blue"),legend=expression(F(A),F(a))) #Add this legend
		}
	if(Graph.Alleles==0 && Graph.Genotypes==1){ #If we only are graphing genotype frequencies
		legend(x=40,y=1.2,fill=c("orange","green","purple"),legend=expression(F(AA),F(Aa),F(aa)),horiz=TRUE)
		}	#Add the one above
	if(Graph.Alleles==1 && Graph.Genotypes==1){ #If we are graphing both
		legend(x=0,y=1.2,fill=c("red","blue"),legend=expression(F(A),F(a))) #Add this legend
		legend(x=generations*0.25,y=1.2,fill=c("orange","green","purple"), #Note the split here.
		legend=expression(F(AA),F(Aa),F(aa)),horiz=TRUE)
		} #And the legend above as well.
	for(j in 1:generations){ #For each generation
		#We start out by looking at how drift affects frequencies.We get the new allele frequency 
		#of A using a binomial distribution. The below simulates 2*PopSize (assume diploid) 
		#alleles from a binomial distribution. Alles are either zero, or one, and the ones are 
		#summed up to give us the new frequency of A after drift.
		A <- sum(rbinom(n=2*PopSize,size=1,prob=A))/(2*PopSize)
		a <- 1-A #And the frequency of a is calculated also
		#Below models natural selection for changing the allele frequencies
		#The equations for calculating new allele frequencies can be found in:
		#Halliburton 2004, p. 136. Note wAA, wAa, and waa are relative fitnesses
		#We will need to divide by wbar to get the allele frequencies to equal one
		wbar <- A*A*wAA + 2*A*a*wAa + a*a*waa
		#Below is the new frequency of A after selection
		A <- (A*A*wAA + A*a*wAa)/wbar
		a <- 1-A #Again, a is calculated by subtracting from one.
		#Genotype frequencies are calculated below
		AA <- A*A*wAA/wbar #Note that they should all sum to one.
		Aa <- 2*A*a*wAa/wbar
		aa <- a*a*waa/wbar
		#After drift and selection, we move on to migration:
		M.AA.Pr <- M.AA/PopSize #Figure out the proportion of AA migrants
		M.Aa.Pr <- M.Aa/PopSize #Figure out the proportion of Aa migrants
		M.aa.Pr <- M.aa/PopSize #Figure out the proportion of aa migrants
		AA <- (M.AA.Pr+AA)/(M.AA.Pr+AA+M.Aa.Pr+Aa+M.aa.Pr+aa) #New AA proportion
		Aa <- (M.Aa.Pr+Aa)/(M.AA.Pr+AA+M.Aa.Pr+Aa+M.aa.Pr+aa) #New Aa proportion
		aa <- (M.aa.Pr+aa)/(M.AA.Pr+AA+M.Aa.Pr+Aa+M.aa.Pr+aa) #New aa proportion
		#Below figures out the new allele frequencies after migration
		A <- ((2*AA)+Aa)/(2*sum(AA,Aa,aa)) #For the A allele
		a <- ((2*aa)+Aa)/(2*sum(AA,Aa,aa)) #For the a allele
		#Below sticks the frequencies and genotypes for this generation in a vector
		#This isn't necessary for the loop, but it allows us to print them out, if desired.
		FreqA[j] <- A; Freqa[j] <- a
		FreqAA[j] <- AA; FreqAa[j] <- Aa; Freqaa[j] <- aa
		#Below plots the points for both alleles, if Graph.Alleles=TRUE
		if(Graph.Alleles==1){
			points(x=j,y=A,pch=20,col="red")
			points(x=j,y=a,pch=20,col="blue")
			}
		#Below plots the points for genotypes
		if(Graph.Genotypes==1){	
			points(x=j,y=AA,pch=20,col="orange")
			points(x=j,y=Aa,pch=20,col="green")
			points(x=j,y=aa,pch=20,col="purple")
			}
		}#Below replots all of the points again, but connects them nicely for easy viewing.
		if(Graph.Alleles==1){ #The connection is done with the type command: "l" for "line"
			points(x=1:generations,y=FreqA,pch=20,col="red",type="l")
			points(x=1:generations,y=Freqa,pch=20,col="blue",type="l")
			}
		if(Graph.Genotypes==1){	
			points(x=1:generations,y=FreqAA,pch=20,col="orange",type="l")
			points(x=1:generations,y=FreqAa,pch=20,col="green",type="l")
			points(x=1:generations,y=Freqaa,pch=20,col="purple",type="l")
		}
	Frequencies <- cbind(FreqA,Freqa,FreqAA,FreqAa,Freqaa)#Bind the frequencies together.
	colnames(x=Frequencies) <- c("F(A)","F(a)","F(AA)","F(Aa)","F(aa)") #Column names
	if(Show.Data==1){ #If you've asked to show the data in the program.
		return(Frequencies) #Spits out the frequencies shown in the graph over time.
		}
	}	