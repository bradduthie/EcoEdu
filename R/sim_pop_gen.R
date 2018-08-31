#' Simulate allele and genotype frequencies at a single locus
#' 
#' This function simulates population genetics at a single locus by modelling the change in allele frequencies over generations as a consequence of genetic drift, migration, and natural selection. Arguments allow users to specify starting allele or genotype frequencies, and different genotype fitnesses or population sizes.
#' 
#'@param A The initial frequency of the 'A' allele
#'@param a The initial frequency of the 'a' allele
#'@param AA The initial number of individuals of the 'AA' genotype
#'@param Aa The initial number of individuals of the 'Aa' genotype
#'@param aa The initial number of individuals of the 'aa' genotype
#'@param generations The number of simulations for which the simulation will run
#'@param PopSize The size of the population being simulated (this will remain constant over generations)
#'@param w_AA The relative fitness of the 'AA' genotype
#'@param w_Aa The relative fitness of the 'Aa' genotype
#'@param w_aa The relative fitness of the 'aa' genotype
#'@param M_AA The number of individuals of genotype 'AA' migrating into the population each generation
#'@param M_Aa The number of individuals of genotype 'Aa' migrating into the population each generation
#'@param M_aa The number of individuals of genotype 'aa' migrating into the population each generation
#'@param plot_alleles Whether or not allele frequencies should be plotted in each generation
#'@param plot_genotypes Whether or not genotype frequencies should be plotted in each generation
#'@param print_data Whether or not the allele and genotype frequencies should be returned in the form of a table as output
#'@param new_plot Whether a new window should be opened to plot
#'@return A plot showing allele frequencies, genotype frequencies, or both, plotted over generations. Output can also include data on allele frequencies and genotype frequencies in each generation.
#'@examples
#'sim_pop_gen(A = 0.7, a = 0.3);
#'@export
sim_pop_gen <- function(A              = 0.5,   # Freq. of the the 'A' allele
                        a              = 0.5,   # Freq. of the the 'a' allele
                        AA             = 0,     # Number of the 'AA' genotype
                        Aa             = 0,     # Number of the 'AA' genotype
                        aa             = 0,     # Number of the 'aa' genotype
                        generations    = 150,   # Number of generations
                        PopSize        = 100,   # Population size
                        w_AA           = 1,     # 'AA' relative fitness
                        w_Aa           = 1,     # 'Aa' relative fitness
                        w_aa           = 1,     # 'aa' relative fitness
				        M_AA           = 0,     # Migrating 'AA' genotypes
				        M_Aa           = 0,     # Migrating 'Aa' genotypes
				        M_aa           = 0,     # Migrating 'aa' genotypes
				        plot_alleles   = TRUE,  # Show allele freq. on plot
				        plot_genotypes = FALSE, # Show genotypes on plot
				        print_data     = FALSE, # Print allele & genotype freqs
				        new_plot       = TRUE   # Add a new plot
				        ){
    # First call up a new window in which to plot allele/genotype frequencies
	if(new_plot == TRUE){
		X11(width = 6, height = 6);
	}
	# Error message if initial allele frequencies don't sum to one
	if(A + a != 1){ 
		plot(x = 0, y = 0, xlim = c(0, generations), ylim=c(0, 1.2), type="n",
		     xlab = "Generation", ylab = "Frequency");
		legend(x = 0, y = 0.75, "Allele frequencies must sum to 1", bty = "n");
		stop("Allele frequencies must sum to 1");
	}
	# Error message if relative survival is greater than one
	if(w_AA > 1 || w_Aa > 1 || w_aa > 1){ 
		plot(x = 0, y = 0, xlim = c(0, generations), ylim = c(0, 1.2), 
		     type = "n", xlab = "Generation", ylab = "Frequency");
		legend(x = 0, y = 0.75,
		       "A relative survival value cannot be greater than 1", bty = "n");
		stop("A relative survival value cannot be greater than 1");
	}
	#Below stops the simulation if the population size is too large.
	if(PopSize > 1000000){
		plot(x = 0, y = 0, xlim = c(0, generations), ylim = c(0, 1.2),
		     type = "n", xlab = "Generation", ylab = "Frequency");
		legend(x = 0, y = 0.75, "Population size cannot be higher than 1000000",
		       bty="n");
		stop("Population size cannot be higher than 1000000");
	}
	if(AA == 0 && Aa == 0 && aa == 0){ # If genotypes aren't set
		AA <- A * A * PopSize;         # Then set to HWE
		Aa <- 2 * A * a * PopSize; 
		aa <- a * a * PopSize;
	} else{ # If genotypes are set, over-ride allele frequencies
		A <- ( (2 * AA) + Aa ) / ( 2 * sum(AA, Aa, aa) );
		a <- ( (2 * aa) + Aa ) / ( 2 * sum(AA, Aa, aa) );
	}
	if( (AA + Aa + aa) != PopSize ){ # Over-ride PopSize if sum of genos differs
		Popsize <- AA + Aa + aa;
	}
    FreqA  <- NULL;
    Freqa  <- NULL;
    FreqAA <- NULL;
	FreqAA <- NULL;
    FreqAa <- NULL;
    Freqaa <- NULL;
	# Initialise a plot that will be used later.
	plot(x = 0, y = 0, xlim = c(0, generations), ylim = c(0, 1.2),
	     type = "n", xlab = "Generation", ylab = "Frequency");
	if(plot_alleles == 1 && plot_genotypes == 0){ 
		legend(x = 0, y = 1.2, fill = c("red", "blue"), 
		       legend = expression(F(A),F(a)) );
	}
	if(plot_alleles == 0 && plot_genotypes == 1){
		legend(x = 40, y = 1.2, fill = c("orange", "green", "purple"),
		       legend = expression(F(AA),F(Aa),F(aa)), horiz = TRUE );
	}
	if(plot_alleles == 1 && plot_genotypes == 1){
		legend(x = 0, y = 1.2, fill = c("red", "blue"),
		       legend = expression(F(A),F(a)) );
		legend(x = generations * 0.25, y = 1.2, 
		       fill = c("orange", "green", "purple"), 
		       legend = expression(F(AA),F(Aa),F(aa)), horiz = TRUE);
	}
	for(j in 1:generations){ 
		A <- sum( rbinom(n = 2*PopSize, size = 1, prob = A) ) / (2 * PopSize);
		a <- 1 - A;
		wbar <- (A * A * w_AA) + (2 * A * a * w_Aa) + a*a*w_aa;
		# Below is the new frequency of A after selection
		A <- ((A * A * w_AA) + (A * a * w_Aa)) / wbar;
		a <- 1 - A;
		#Genotype frequencies are calculated below
		AA <- (A * A * w_AA) / wbar;
		Aa <- (2 * A * a * w_Aa) / wbar;
		aa <- (a * a * w_aa) / wbar;
		#After drift and selection, we move on to migration:
		M_AA_Pr <- M_AA / PopSize; # Proportion of AA migrants
		M_Aa_Pr <- M_Aa / PopSize; # Proportion of Aa migrants
		M_aa_Pr <- M_aa / PopSize; # Proportion of aa migrants
		AA <- (M_AA_Pr + AA) / (M_AA_Pr + AA + M_Aa_Pr + Aa + M_aa_Pr + aa); 
		Aa <- (M_Aa_Pr + Aa) / (M_AA_Pr + AA + M_Aa_Pr + Aa + M_aa_Pr + aa);
		aa <- (M_aa_Pr + aa) / (M_AA_Pr + AA + M_Aa_Pr + Aa + M_aa_Pr + aa);
		#Below figures out the new allele frequencies after migration
		A <- ( (2 * AA) + Aa) /(2 * sum(AA, Aa, aa)); 
		a <- ( (2 * aa) + Aa) /(2 * sum(AA, Aa, aa));
		FreqA[j]  <- A; 
		Freqa[j]  <- a;
		FreqAA[j] <- AA; 
		FreqAa[j] <- Aa; 
		Freqaa[j] <- aa;
		#Below plots the points for both alleles, if plot_alleles=TRUE
		if(plot_alleles == TRUE){
			points(x = j, y = A, pch = 20, col = "red");
			points(x = j, y = a, pch = 20, col = "blue");
		}
		#Below plots the points for genotypes
		if(plot_genotypes == TRUE){	
			points(x = j, y = AA, pch = 20, col = "orange");
			points(x = j, y = Aa, pch = 20, col = "green");
			points(x = j, y = aa, pch = 20, col = "purple");
			}
		}
		if(plot_alleles == TRUE){
			points(x = 1:generations, y = FreqA, pch = 20, col = "red", 
			       type = "l");
			points(x = 1:generations, y = Freqa, pch = 20, col = "blue",
			       type = "l");
			}
		if(plot_genotypes == TRUE){	
			points(x = 1:generations, y = FreqAA, pch = 20, col = "orange",
			       type = "l");
			points(x = 1:generations, y = FreqAa, pch = 20, col = "green",
			       type = "l");
			points(x = 1:generations, y = Freqaa, pch = 20, col = "purple",
			       type = "l");
		}
	Frequencies <- cbind(FreqA, Freqa, FreqAA, FreqAa, Freqaa);
	colnames(x = Frequencies) <- c("F(A)", "F(a)", "F(AA)", "F(Aa)", "F(aa)");
	if(print_data == TRUE){ 
		return(Frequencies);
	}
}