package gb.esac.periodogram;

import Jama.Matrix;
import gb.esac.eventlist.EventList;


public class SingleFrequencyPowerEstimator {

    public static void main(String[] args) throws Exception {

	if ( args.length != 3 ) {
	    System.out.println("Usage: java SingleFrequencyPowerEstimator evlist frequency harmonic");
	    System.exit(-1);
	}
	String filename = args[0];
	double frequency = (new Double(args[1])).doubleValue();
	int harmonic = (new Integer(args[2])).intValue();

	EventList evlist = new EventList(filename);
	double period = 1d/frequency;
	double power = getCorrectedPowerForThisHarmonic(evlist.getArrivalTimes(), period, harmonic);
	System.out.println(frequency+"\t"+power);
    }

    
    public static double getCorrectedPowerForThisHarmonic(double[] times, double period, int harmonic) {

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

		
	//   Get the phase corresponding to each time  
	double[] phases = getPhases(times, period);

		
	//   Calculate C and S  
	double meanCos = 0;
	double meanSin = 0;
	for ( int i=0; i < nevents; i++ ) {
	    meanCos += Math.cos(harmonic*2*Math.PI*phases[i])/nevents;
	    meanSin += Math.sin(harmonic*2*Math.PI*phases[i])/nevents;
	}
		
		
	//   Calculate [C - E(C)] and [S - E(S)] 
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
	double[][] covMatrixValues = new double[][] { {varianceOfCos, covarianceOfCosSin} , {covarianceOfCosSin, varianceOfSin} };
	//logger.info("Cov Matrix:");
	// logger.info("{"+varianceOfCos+", "+covarianceOfCosSin+"}");
	//logger.info("{"+covarianceOfCosSin+", "+varianceOfCos+"}");
	//logger.info();
	Matrix covarianceMatrix = new Matrix(covMatrixValues);
	Matrix inverseCovMatrix = covarianceMatrix.inverse();
	Matrix transposeCosAndSinMatrix = cosAndSinMatrix.transpose();
				
	//   Calulate the Modified Rayleigh Power by matrix multiplication  	
	Matrix powerMatrix = transposeCosAndSinMatrix.times(inverseCovMatrix).times(cosAndSinMatrix);
	double power = powerMatrix.det();
		
	return power;

    }

    public static double[] getPhases(double[] times, double period) {

	double[] phases = new double[times.length];
	double tOverP = 0;
	for ( int i=0; i < times.length; i++ ) {
	    tOverP = times[i]/period;
	    phases[i] = tOverP - Math.floor(tOverP);
	}
	return phases;
    }


}
