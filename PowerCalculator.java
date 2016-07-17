package gb.esac.periodogram;

import Jama.Matrix;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister64;
import cern.jet.stat.Probability;
import gb.esac.tools.BasicStats;
import gb.esac.tools.DataUtils;
import org.apache.log4j.Logger;

public final class PowerCalculator {

    private static Logger logger  = Logger.getLogger(PowerCalculator.class);

    public static double[] getCorrectedPowerComponentsForThisHarmonic(double[] times, double period, int harmonic) {
	//  This is a generalised version of the Modified Rayleigh power for any number of harmonics
	//  The formalism is identical but omega is replaced by harmonic*omega (kOmega)

	int nevents = times.length;
	//   Implementation of equations (5) to (9)  
	double tfirst = times[0];
	double tlast = times[nevents-1];
	double duration = tlast - tfirst;
	double omega = 2*Math.PI/period;
	double kOmega = harmonic*omega;
	double expectedMeanCos = (Math.sin(kOmega*tlast) - Math.sin(kOmega*tfirst)) / (kOmega*duration);
	double expectedMeanSin = (Math.cos(kOmega*tfirst) - Math.cos(kOmega*tlast)) / (kOmega*duration);
	double crossTerm = ( Math.sin(kOmega*tlast)*Math.cos(kOmega*tlast) - Math.sin(kOmega*tfirst)*Math.cos(kOmega*tfirst) ) / (2*kOmega*duration);
	double varianceOfCos = (0.5 + crossTerm - Math.pow(expectedMeanCos, 2))/nevents;
	double varianceOfSin = (0.5 - crossTerm - Math.pow(expectedMeanSin, 2))/nevents;
	double SinSqrdTerm = ( Math.pow(Math.sin(kOmega*tlast), 2) - Math.pow(Math.sin(kOmega*tfirst), 2) ) / (2*kOmega*duration);
	double covarianceOfCosSin = (SinSqrdTerm - expectedMeanCos*expectedMeanSin)/nevents;
	//logger.info(expectedMeanCos+"\t"+expectedMeanSin+"\t"+crossTerm);
	//   Calculate C and S  
	double meanCos = 0;
	double meanSin = 0;
	for ( int i=0; i < nevents; i++ ) {
	    double tOverP = times[i]/period;
	    double phase = tOverP - Math.floor(tOverP);
	    phase *= 2*Math.PI;
	    meanCos += Math.cos(harmonic*phase)/nevents;
	    meanSin += Math.sin(harmonic*phase)/nevents;
	}
	//   Calculate (C - <C>) and (S - <S>) 
	double cosDiff = meanCos - expectedMeanCos;
	double sinDiff = meanSin - expectedMeanSin;
	//logger.info(cosDiff+"\t"+sinDiff);		

	//   Get power from cross terms  
	// 	double power = //Math.pow(meanCos,2) + Math.pow(meanSin,2) 
	// 	    + 2*(meanCos*expectedMeanCos + meanSin*expectedMeanSin) 
	// 	    - (Math.pow(expectedMeanCos,2) + Math.pow(expectedMeanSin,2));
	// 	power *= 2*nevents;
		
	//   Construct Matrices to calculate power  
	double[][] cosAndSinMatrixValues = new double[][] { {cosDiff} , {sinDiff} };
	Matrix cosAndSinMatrix = new Matrix(cosAndSinMatrixValues);
	double[][] covMatrixValues = new double[][] 
	    { {varianceOfCos, covarianceOfCosSin} , 
	      {covarianceOfCosSin, varianceOfSin} 
	    };
	//logger.info("Cov Matrix:");
	// logger.info("{"+varianceOfCos+", "+covarianceOfCosSin+"}");
	//logger.info("{"+covarianceOfCosSin+", "+varianceOfCos+"}");
	//logger.info();
	Matrix covarianceMatrix = new Matrix(covMatrixValues);
	//System.out.print("Cov");
	//covarianceMatrix.print(2, 10);
	Matrix inverseCovMatrix = covarianceMatrix.inverse();
	//inverseCovMatrix.print(2, 10);
	Matrix transposeCosAndSinMatrix = cosAndSinMatrix.transpose();
	Matrix powerMatrix = transposeCosAndSinMatrix.times(inverseCovMatrix).times(cosAndSinMatrix);
	double power = powerMatrix.det();
	return new double[] {power, meanCos, meanSin, expectedMeanCos, expectedMeanSin, varianceOfCos, varianceOfSin, covarianceOfCosSin};
    }

    public static double getCorrectedPowerForThisHarmonic(double[] times, double period, int harmonic) {
	return getCorrectedPowerComponentsForThisHarmonic(times, period, harmonic)[0];
    }

    public static double getModRayleighPower(double[] times, double period) {
	//  The Modified Rayleigh Power is computed according to 
	//     K. J. Orford et al. 1996, in Astroparticle Physics 4 (1996) 235-239
	//  for the fundamental harmonic, 
	//  and was generalised to any harmonic by me.
 	//  See method getCorrectedPowerForThisHarmonic above for implementation details
	int harmonic = 1;
 	return getCorrectedPowerForThisHarmonic(times, period, harmonic);
    }

    public static double getModRayleighPower(double[] times, double[] rates, double period, double mean) {
	int harmonic = 1;
	double[] errors = new double[rates.length];
	for ( int i=0; i < rates.length; i++ ) errors[i] = 1.0;
	return getModRayleighPower(times, rates, errors, period, mean, harmonic);
    }
	
    public static double getModRayleighPower(double[] times, double[] rates, double[] errors, double period, double mean) {
	int harmonic = 1;
	return getModRayleighPower(times, rates, errors, period, mean, harmonic);
    }

    public static double getModRayleighPower(double[] times, double[] rates, double[] errors, double period, double mean, int harmonic) {
	// logger.info("Calculating modified Rayleigh power on binned time series");
	// logger.info("  Period = "+period);
	// logger.info("  Harmonic = "+harmonic);
	// logger.info("  Number of bins = "+times.length);
	// logger.warn("  Times assumed to start at 0");

	// TIMES ASSUMED TO START AT ZERO

	int k = harmonic;
	double c = 0;
	double s = 0;
	double sumOfCSqrd = 0;
	double sumOfSSqrd = 0;
	double cov = 0;
	for ( int i=1; i < times.length; i++ ) {
	    double tOverP = times[i]/period;
	    double phase = tOverP - Math.floor(tOverP);
	    phase *= 2*Math.PI;
	    double cosPhi = Math.cos(k*phase);
	    double sinPhi = Math.sin(k*phase);
	    double meanSubRate = rates[i] - mean;
	    double weight = meanSubRate/Math.pow(errors[i],2);
	    //logger.info(rates[i]+" "+errors[i]);
	    c += cosPhi*weight;
	    s += sinPhi*weight;
	    sumOfCSqrd += Math.pow(cosPhi*weight, 2);
	    sumOfSSqrd += Math.pow(sinPhi*weight, 2);
	}
	double power = Math.pow(c, 2)/sumOfCSqrd + Math.pow(s, 2)/sumOfSSqrd;
	return power;
    }
	
    public static double getLombPower(double[] times, double[] rates, double period, double mean, double variance) {
	// logger.info("Calculating Lomb power on binned time series");
	// logger.info("  Period = "+period);
	// logger.info("  Number of bins = "+times.length);
	// logger.warn("  Times assumed to start at 0");
	// TIMES ASSUMED TO BE RESET TO ZERO	
	// Note: mean and variance are given as arguments to minimise loops in each method call
	//       especially because the method is generally called many times for the same time series
	//       in order to contruct the periodogram, but the mean and variance only need to be computed once.
	double sumOfSin2wt = 0;
	double sumOfCos2wt = 0;
	double w = 2*Math.PI/period;
	//   Calculate tau for the given test frequency  
	for ( int i=0; i < times.length; i ++ ) {
	    sumOfSin2wt += Math.sin(2*w*times[i]);
	    sumOfCos2wt += Math.cos(2*w*times[i]);
	}
	double tau = Math.atan(sumOfSin2wt/sumOfCos2wt)/(2*w);
	//   Calculate the cos and sin components of power
	double sumOfWeightedCos = 0;
	double sumOfWeightedSin = 0;
	double sumOfSqrdCos = 0;
	double sumOfSqrdSin = 0;
	for ( int i=0; i < times.length; i ++ ) {
	    sumOfWeightedCos += (rates[i] - mean)*Math.cos(w*(times[i] - tau));
	    sumOfWeightedSin += (rates[i] - mean)*Math.sin(w*(times[i] - tau));
	    sumOfSqrdCos += Math.pow(Math.cos(w*(times[i] - tau)), 2);
	    sumOfSqrdSin += Math.pow(Math.sin(w*(times[i] - tau)), 2);
	}
	//   Calculate Lomb-Scargle Power 
	double norm = 1./(2*variance);  //  Standard norm. from Num. Recipes
	//norm *= 2; // Scale by a factor of 2 to get Leahy normalisation
	double power = norm * (Math.pow(sumOfWeightedCos, 2)/sumOfSqrdCos + Math.pow(sumOfWeightedSin, 2)/sumOfSqrdSin);
	return power;
    }

    public static double getModifiedZ2Power(double[] times, double period, int nHarm) {
	// logger.info("Calculating modified Z2 power on unbinned time series (arrival times)");
	// logger.info("  Period = "+period);
	// logger.info("  Number of harmonics = "+nHarm);
	// logger.info("  Number of events = "+times.length);
	// logger.warn("  Times assumed to start at 0");
	// TIMES ASSUMED TO BE RESET TO ZERO
	double z2 = 0;
	for ( int k=1; k <= nHarm; k++ ) {
	    z2 += getCorrectedPowerForThisHarmonic(times, period, k);
	}
	// logger.info("Power = "+z2);
	return z2;
    }

    public static double getZ2Power(double[] times, double period, int nHarm) {
	// logger.info("Calculating classical Z2 power on unbinned time series (arrival times)");
	// logger.info("  Period = "+period);
	// logger.info("  Number of harmonics = "+nHarm);
	// logger.info("  Number of events = "+times.length);
	// logger.warn("  Times assumed to start at 0");
	// TIMES ASSUMED TO BE RESET TO ZERO
	double c = 0;
	double s = 0;
	double z2 = 0;
	double nTot = (new Double(times.length)).doubleValue();
	for ( int k=0; k < nHarm; k++ ) {
	    for ( int i=0; i < nTot; i++ ) {
		double tOverP = times[i]/period;
		double phase = tOverP - Math.floor(tOverP);
		phase *= 2*Math.PI*(k+1);
		c += Math.cos(phase);
		s += Math.sin(phase);
	    }
	    z2 += (2/nTot)*(Math.pow(c, 2) + Math.pow(s, 2));
	}
	// logger.info("Power = "+z2);
	return z2;
    }

    public static double[] getZ2statsFromPhases(double[] phases, int nHarm) {
	double nTot = (new Double(phases.length)).doubleValue();
	double[] c = new double[nHarm]; 
	double[] s = new double[nHarm];
	double z2 = 0;
	for ( int k=0; k < nHarm; k++ ) {
	    for ( int i=0; i < nTot; i++ ) {
		double phase = 2*Math.PI*phases[i]*(k+1);
		c[k] += Math.cos(phase);
		s[k] += Math.sin(phase);
	    }
	    z2 += Math.pow(c[k], 2) + Math.pow(s[k], 2);
	}
	z2 *= 2/nTot;
	double[] z2stats = new double[] {c[0], s[0], z2};
	return z2stats;
    }
	
    public static double[] getZ2stats(double[] times, double period, int nHarm) {
	// TIMES ASSUMED TO BE RESET TO ZERO
	double nTot = (double) times.length;		
	//  Determine number of cycles and number of events in last cycle
	double cycles = (times[times.length-1] - times[0])/period;
	int ncycles = (int) Math.floor(cycles);
	double remainder = cycles - ncycles;
	double timeAtLastCycle = period*ncycles + times[0];
	//int eventsInLastCycle = 0;
	//for ( int i=0; i < times.length; i++ ) {
	//if ( times[i] > timeAtLastCycle ) eventsInLastCycle++;

	    //logger.info(times[i]+" "+timeAtLastCycle+" "+eventsInLastCycle);
	//}
	//double nEventsUsed = nTot - eventsInLastCycle;
	//logger.info(cycles+" "+ncycles+" "+timeAtLastCycle+" "+nEventsUsed);
	double[] c = new double[nHarm]; 
	double[] s = new double[nHarm];
	double z2 = 0;
		
	//  CASE 1
	//  Construct phases using all events. Each array element is the sum on 1 harm
	for ( int k=0; k < nHarm; k++ ) {
	    for ( int i=0; i < nTot; i++ ) {
		double tOverP = times[i]/period;
		double phase = tOverP - Math.floor(tOverP);
		phase *= 2*Math.PI*(k+1);
		c[k] += Math.cos(phase);
		s[k] += Math.sin(phase);
	    }
	    z2 += (2/nTot)*(Math.pow(c[k], 2) + Math.pow(s[k], 2));
	}
		
	// //  CASE 2
	// //  Construct phases excluding the events in the last partial cycle
	// for ( int k=0; k < nHarm; k++ ) {
 	//     for ( int i=0; i < nEventsUsed; i++ ) {
	// 	double tOverP = times[i]/period;
	// 	double phase = tOverP - Math.floor(tOverP);
	// 	phase *= 2*Math.PI*(k+1);
	// 	c[k] += Math.cos( phase );
	// 	s[k] += Math.sin( phase );
	//     }
	//     z2 += (2/nEventsUsed)*(Math.pow(c[k], 2) + Math.pow(s[k], 2));
	// }
		
	// //  CASE 3
	// //  Construct phases by using the first events to complete the last partial cycle
	// for ( int k=0; k < nHarm; k++ ) {
	//     for ( int i=0; i < nTot; i++ ) {
	// 	if ( i < nEventsUsed ) {
	// 	    phase = times[i]/period; 
	// 	}
	// 	else {
	// 	    phase = times[i - nEventsUsed]/period; 
	// 	}
	// 	phase -= Math.floor(phase);
	// 	phase *= 2*Math.PI;
	// 	c[k] += Math.cos((k+1)*phase);
	// 	s[k] += Math.sin((k+1)*phase);
	//     }
	//     z2 += (2/nTot)*(Math.pow(c[k], 2) + Math.pow(s[k], 2));
	// }
		
	// 	//  Normalize alpha[0] and beta[0] ( first harmonic )
	// 	double alpha_0 = alpha[0]*Math.sqrt(2/nTot);
	// 	double beta_0 = beta[0]*Math.sqrt(2/nTot);
		
	//  Return alpha, beta and z2
	double[] z2stats = new double[] {c[0], s[0], z2};
	return z2stats;
    }
	
    public static double[] getZ2stats(double[] times, double[] rate, double[] error, double period, int nHarm) {
		
	// TIMES ASSUMED TO BE RESET TO ZERO	
		
	//   Determine number of cycles and number of events in last cycle  
	double cycles = (times[times.length-1] - times[0])/period;
	int ncycles = (new Double(Math.floor(cycles))).intValue();
	double remainder = cycles - ncycles;
	double timeAtLastCycle = period*ncycles + times[0];
	double nTot = (new Double(times.length)).doubleValue();
	//int eventsInLastCycle = 0;
	// 	for ( int i=0; i < times.length; i++ ) {
	// 	    if ( times[i] > timeAtLastCycle ) eventsInLastCycle++;
	// 	    //logger.info(times[i]+" "+timeAtLastCycle+" "+eventsInLastCycle);
	// 	}
	//int nEventsUsed = (new Double(nTot - eventsInLastCycle)).intValue();
	//logger.info(cycles+" "+ncycles+" "+timeAtLastCycle+" "+nEventsUsed);
		
	//  Define and initialize alpha and beta
	double[] c = new double[nHarm];  // k cos coefficients for k harm
	double[] s = new double[nHarm];    // k sin coefficients for k harm
	double[] cErr = new double[nHarm]; 
	double[] sErr = new double[nHarm]; 
	//  Determine mean rate
	double mean = BasicStats.getWMean(rate, error);
	//  Construct phases, alpha and beta using all events
	double z2 = 0;
	int nUsed = 0;
	for ( int k=0; k < nHarm; k++ ) {
	    for ( int i=0; i < nTot; i++ ) {
		double tOverP = times[i]/period;
		double phase = tOverP - Math.floor(tOverP);
		phase *= 2*Math.PI*(k+1);
		if ( ! Double.isNaN(rate[i]) ) {
		    c[k] += rate[i]*Math.cos(phase);
		    s[k] += rate[i]*Math.sin(phase);
		    cErr[k] += Math.pow(error[i]*Math.cos(phase), 2);
		    sErr[k] += Math.pow(error[i]*Math.sin(phase), 2);
		    nUsed++;
		}
	    }
	    z2 += Math.pow(c[k], 2)/cErr[k] + Math.pow(s[k], 2)/sErr[k];
	}
	//  Return c[0], s[0] and z2
	double[] z2stats = new double[] {Math.sqrt(c[0]), Math.sqrt(s[0]), z2};
	return z2stats;
    }
	
    public static double getZ2prob(double[] times, double period, int nHarm) {
	double z2stats[] = getZ2stats(times, period, nHarm);
	double z2 = z2stats[2];
	double dof = nHarm*2;
	double z2Prob = Probability.chiSquareComplemented(dof, z2);
	return z2Prob;
    }
	
	
    public static double getH_value(double[] times, double period) {
	// 	//  Determine optimum number of harmonics
	// 	double[] z2_values = new double[20]; 
	// 	double[] h_values = new double[20];
	// 	for ( int nHarm=0; nHarm < 20; nHarm++ ) {
	// 	    z2_values[nHarm] = getZ2_value(times, period, nHarm+1);
	// 	    h_values[nHarm] = z2_values[nHarm] - 4*(nHarm+1) + 4;
	// 	    //logger.info("nHarm = "+(nHarm+1)+
	// 	    //	       ": Z2 = "+z2_values[nHarm]+
	// 	    //	       "\t H = "+h_values[nHarm]);
	// 	}
	// 	double maxH = MinMax.getMax(h_values);
	// 	for ( int i=0; i < h_values.length; i++ ) logger.info(i+" "+h_values[i]+" "+maxH);
	// 	//int m = Arrays.binarySearch(h_values, maxH);
	// 	//logger.info(m+1);
	// 	return h_values[0];
		
	int m = 3;
	double[] z2stats = getZ2stats(times, period, m);
	double z2_value = z2stats[2];
	double h_value = z2_value - 4*m + 4;
	return h_value;
    }
	
    public static double getHprob(double[] times, double period, int nHarm) {
	double[] z2stats = getZ2stats(times, period, nHarm);
	double z2 = z2stats[2];
	double h_value = z2 - 4*nHarm + 4;
	//  Calculate probability associated with that H-value
	double a = 0, b = 0, c = 0;
	if ( h_value <= 23 ) {
	    a = 0.9999755;  b = -0.39802;  c = 0;
	}
	else if ( h_value > 23 && h_value <= 50 ) {
	    a = 1.210597;  b = -0.45901;  c = 0.00229;
	}
	else { a = 0;  b = 0;  c = 0; }
	double prob = a*Math.exp(b*h_value + c*Math.pow(h_value, 2));
	return prob;
    }
}