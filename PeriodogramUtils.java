package gb.esac.periodogram;

import java.text.DecimalFormat;
import java.util.Arrays;

import cern.colt.list.DoubleArrayList;
import gb.esac.tools.Converter;
import gb.esac.tools.LeastSquaresFitter;
import org.apache.log4j.Logger;


/**
 * Class <code>Periodogram</code> 
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC
 * @version 1.0  (last modified: Aug 2015)
 *
 */
public final class PeriodogramUtils {

    private static Logger logger  = Logger.getLogger(PeriodogramUtils.class);
    private static DecimalFormat exp = new DecimalFormat("0.0000#E00");

    /**
     * The method <code>getFFTFrequencies</code> returns the independent frequencies (IF) for the specifies duration and number of bins. The Fourier spacing is defined by df=1/duration. The number of independent Fourier spacings (IFS) in the interval between nuMin and nuMax is given by nIFS = (nuMax-nuMin)/df = duration*(nuMax-nuMin). The minium and maximum frequencies that can be tested are respectively nuMin = 1/duration and the Nyquist frequency given by nuMax = 1/2dt, where dt is the sampling frequency (bin time) of the time series. Given that dt = duration/nbins, we can also write nuMax = nbins/2*duration, and thus only need to specify one of the two quantities nbins and dt. The number of the IFS is then nIFS = duration*(nbins/2*duration - 1/duration) = nbins/2 -1. Therefore the number of IF is nbins/2 given that nuMin and nuMax are respectively the lower and the upper bounds of the testable frequency range.
     *
     * @param nbins an <code>int</code> value that specifies the number of bins in the time series.
     * @param duration a <code>double</code> value that specifies the duration of the data set.
     * @return a <code>double[]</code> value given the independent Fourier frequencies for the duration and number of bins.
     * @exception PeriodogramException if the number of bins is not a power of 2.
     */
    public static double[] getFFTFrequencies(int nbins, double duration) throws PeriodogramException {
	//  Check that nLCBins is a power of 2 >= nOldBins
	double n1 = Math.floor(Math.log10(nbins)/Math.log10(2));
	double n2 = Math.log10(nbins)/Math.log10(2);
	if ( n1 != n2 ) {
	    new PeriodogramException("nbins = "+nbins+" is not a power of 2");
	}
	//  Define the frequencies 
	double nuMin = 1/duration;
 	int nFreqs = nbins/2;
 	double[] frequencies = new double[nFreqs];
	for ( int i=0; i < nFreqs; i++ ) {
	    frequencies[i] = nuMin*(i+1);
	}
	return frequencies;
    }

    /**
     * The method <code>getFourierFrequencies</code> returns the Fourier frequencies in a frequency interval using the specified oversampling. There is oversampling between the nuMin and nuMax that define the lower and upper bounds of the testable frequency range. This means that for the full Fourier range with nuMin=1/duration and nuMax = 1/2dt, the number of test frequencies will be given by nIF * oversampling - oversampling. 
     *
     * @param nuMin a <code>double</code> value that defines the minimum frequency
     * @param nuMax a <code>double</code> value that defines the maximum frequency
     * @param duration a <code>double</code> value that defines the total duration of the data. It is used to set frequency step given by 1/dutation.
     * @param oversampling a <code>int</code> value that defines oversampling factor: the number of test frequencies per IFS.
     * @return a <code>double[]</code> that are the Fourier frequencies between nuMin and nuMax with oversampling.
     */
    public static double[] getFourierFrequencies(double nuMin, double nuMax, double duration, int sampling) {
	logger.info("Defining test frequencies");
	logger.info("  Min = "+nuMin);
	logger.info("  Max = "+nuMax);
	logger.info("  Sampling factor = "+sampling);
	DoubleArrayList freqs = new DoubleArrayList();
	double df = 1d/(duration*sampling);
	double min = 1d/duration;
	// This is required in order to have a total number of test frequencies of nIFS*sampling
	double max = nuMax + (sampling-1)*df;
	double testFreq = min;
	// Skip acceptable test frequencies that are below the requested nuMin
	while ( testFreq < nuMin-0.5*df ) {
	    testFreq += df;
	}
	while ( testFreq < max+0.5*df ) {
	    freqs.add(testFreq);
	    testFreq += df;
	}
	freqs.trimToSize();
	logger.info("  Number of frequencies = "+freqs.size());
	return freqs.elements();
    }

    public static double[] getFourierFrequencies(double nuMin, double nuMax, double df) {
	int nFreqs = (int) Math.round((nuMax-nuMin)/df);
	double[] frequencies = new double[nFreqs];
	for ( int i=0; i < nFreqs; i++ ) {
	    frequencies[i] = nuMin + i*df;
	}
	return frequencies;
    }

    /**
     * The method <code>getFourierPeriods</code> returns the periods that correspond to the frequencies returned by getFourierTestFrequencies, but in reverse order so that the periods are ascending.
     *
     * @param pmin a <code>double</code> value
     * @param pmax a <code>double</code> value
     * @param duration a <code>double</code> value
     * @param sampling a <code>int</code> value
     * @return a <code>double[]</code> that are the test periods ordered from smallest to largest
     */
    public static double[] getFourierPeriods(double pmin, double pmax, double duration, int sampling) {
	double numin = 1/pmax;
	double numax = 1/pmin;
	double[] testFreqs = getFourierFrequencies(numin, numax, duration, sampling); 
	int nFreqs = testFreqs.length;
	double[] testPeriods = new double[nFreqs];
	for ( int i=0; i < nFreqs; i++ ) {
	    testPeriods[i] = 1/testFreqs[nFreqs-1-i];
	}
	return testPeriods;
    }
    
    /**
     * The method <code>getPhases</code> uses the specified period to convert each time to its corresponding phase between 0 and 1.
     *
     * @param times a <code>double[]</code> that specifyies the arrival times.
     * @param period a <code>double</code> that specifies the period used to calculate the phase of each arrival time.
     * @return a <code>double[]</code> that gives the phases between 0 and 1.
     */
    public static double[] getPhases(double[] times, double period) {
	double[] phases = new double[times.length];
	double tOverP = 0;
	for ( int i=0; i < times.length; i++ ) {
	    tOverP = times[i]/period;
	    phases[i] = tOverP - Math.floor(tOverP);
	}
	Arrays.sort(phases);
	return phases;
    }

    // public static double[] getPhases(double[] times, double period, double phaseShift) {
    // 	double[] phases = getPhases(times, period);
    // 	return shiftPhases(phases, phaseShift);
    // }

    public static double[] getPhases(double[] times, double[] periods) throws PeriodogramException {
	if ( times.length != periods.length ) {
	    throw new PeriodogramException("Input arrays must have the same number of elements");
	}
	double[] phases = new double[times.length];
	double tOverP = 0;
	for ( int i=0; i < times.length; i++ ) {
	    tOverP = times[i]/periods[i];
	    phases[i] = tOverP - Math.floor(tOverP);
	}
	Arrays.sort(phases);
	return phases;
    }

    public static double[] getPhases(double[] times, double slope, double intercept) throws PeriodogramException {
	double[] periods = new double[times.length];
	for ( int i=0; i < times.length; i++ ) {
	    periods[i] = slope*times[i] + intercept;
	}
	return getPhases(times, periods);
    }

    public static double[] getPhases(double[] times, double tZero, double nu, double nuDot, double nuDotDot) {
	double[] phases = new double[times.length];
	for ( int i=0; i < times.length; i++ ) {
	    double timeDiff = times[i] - tZero;
	    phases[i] = nu*timeDiff + 1/2*nuDot*Math.pow(timeDiff, 2) + 1/6*nuDotDot*Math.pow(timeDiff, 3);
	    if ( phases[i] < 0 ) {
		double shift = Math.ceil(-phases[i]);
		phases[i] += shift;
	    }
	    else {
		double shift = Math.floor(phases[i]);
		phases[i] -= shift;
	    }
	}
	Arrays.sort(phases);
	return phases;	
    }

    public static double[] getPhases(double[] times, double tZero, double nu, double nuDot, double nuDotDot, double phaseShift) {
	double[] phases = getPhases(times, tZero, nu, nuDot, nuDotDot);
	return shiftPhases(phases, phaseShift);
    }

    public static double[] shiftPhases(double[] phases, double phaseShift) {
	double[] shiftedPhases = new double[phases.length];
	for ( int i=0; i < phases.length; i++ ) {
	    shiftedPhases[i] = phases[i] + phaseShift;
	    if ( shiftedPhases[i] < 0 ) shiftedPhases[i] += 1;
	    if ( shiftedPhases[i] >= 1 ) shiftedPhases[i] -= 1;	    
	}
	Arrays.sort(shiftedPhases);
	return shiftedPhases;
    }

    public static int nIFS(double duration, double nuMin, double nuMax) {
	return (new Double(Math.floor(duration*(nuMax - nuMin)))).intValue();
    }


    public static double[] getIFSEdgesInFreqSpace(double nuMin, double nuMax, double duration, double sampling, double ifsOffset) {
	//  Determine the number of IFS in interval
	int nIFS = (new Double(Math.round(duration*(nuMax - nuMin)))).intValue();
	int nFreqs = (new Double((nIFS+1)*sampling)).intValue();
	//  Define and initialise array of test frequencies
	double[] ifsEdgesInFreqSpace = new double[nFreqs+1];
	double freq = 1/(sampling*duration);
	int n = 0;
	int i = 2;
	while ( n < nFreqs ) {
	    freq = i/(sampling*duration) + ifsOffset/duration;
	    if ( freq >= nuMin && freq <= nuMax ) {
		ifsEdgesInFreqSpace[n] = freq + 0.5/(sampling*duration);
		//System.out.println(n+" "+freq+" "+ifsEdgesInFreqSpace[n]);
		n++;
		i++;
	    }
	    else i++;
	}
	ifsEdgesInFreqSpace[n] = ifsEdgesInFreqSpace[n-1] + 1/(sampling*duration);
	//System.out.println(n+" "+ifsEdgesInFreqSpace[n]); 
	return ifsEdgesInFreqSpace;
    }

    public static double[] getIFSEdgesInPeriodSpace(double pmin, double pmax, double duration, double sampling, double ifsOffset) {
	double numin = 1/pmax;
	double numax = 1/pmin;
	double[] ifsEdgesInFreqSpace = getIFSEdgesInFreqSpace(numin, numax, duration, sampling, ifsOffset);
	int nEdges = ifsEdgesInFreqSpace.length;
	double[] ifsEdgesInPeriodSpace = new double[nEdges];
	for ( int i=0; i < nEdges; i++ ) {
	    ifsEdgesInPeriodSpace[i] = 1/ifsEdgesInFreqSpace[nEdges-1 - i];
	    //System.out.println(i+" "+ifsEdgesInPeriodSpace[i]);
	}
	return ifsEdgesInPeriodSpace;
    }

    /**
     * Method <code>getIntegratedPower</code> calculated the integral of the Periodogram: sum of binWidths[i]*power[i].
     *
     * @param nuStart a <code>double</code> that defines the lower integration limit
     * @param nuStop a <code>double</code> that defines the upper integration limit
     * @return a <code>double</code> the integral of the Periodogram between nuStart and nuStop
     */
    public static double getIntegratedPower(Periodogram periodogram, double nuStart, double nuStop) {
	double integral = 0;
	// Skip all frequencies before nuStart
	int i=0;
	double[] freqs = periodogram.getFreqs();
	while ( freqs[i] < nuStart ) 
	    i++;
	// Integrate up to nuStop
	int nBins = periodogram.nBins();
	double[] binWidths = periodogram.getBinWidths();
	double[] powers = periodogram.getPowers();
	while ( i < nBins && freqs[i] <= nuStop ) {
	    integral += binWidths[i]*powers[i];
	    i++;
	}
	return integral;
    }

    public static double getIntegratedPower(Periodogram periodogram) {
	return getIntegratedPower(periodogram, periodogram.nuMin(), periodogram.nuMax());
    }

    /**
     * The method <code>fitPowerLawInLogSpace</code> performs a least squares fit to the periodogram in log-log space, i.e., a line y=mx+b where y if the power, and x is the frequency.
     *
     * @return a <code>double[]</code> that contains 4 numbers: the best fit slope, error on slope, ordinate and error on ordinate in this order.
     */
    public static double[] fitPowerLawInLogSpace(Periodogram psd) {
	logger.info("Fitting power-law in log space");
	double[] freqsInLogSpace = Converter.lin2logSpace(psd.getFreqs());
	double[] powsInLogSpace =  Converter.lin2logSpace(psd.getPowers());
	return LeastSquaresFitter.leastSquaresFitLine(freqsInLogSpace, powsInLogSpace);
    }

    /**
     * The method <code>fitPowerLawInLinearSpace</code> performs a least squares fit to the peridogram in linear space, i.e., a power law y=n*x^alpha.
     *
     * @return a <code>double[]</code> that contains 2 numbers: the best fit index and the normalization constant in this order.
     */
    public static double[] fitPowerLawInLinearSpace(Periodogram psd) {
	logger.info("Fitting power-law in linear space");
	double[] fitResult = LeastSquaresFitter.leastSquaresFitPowerLaw(psd.getFreqs(), psd.getPowers());
	double index = fitResult[0];
	double normalization = fitResult[1];
	return new double[] {index, normalization};
    }
}
