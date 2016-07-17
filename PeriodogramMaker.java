package gb.esac.periodogram;

//import gb.esac.binner.Binner;
//import gb.esac.io.AsciiDataFileWriter;
import cern.colt.list.DoubleArrayList;
import gb.esac.binner.BinningException;
import gb.esac.eventlist.EventList;
import gb.esac.likelihood.ExponentialLikelihood;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesException;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesResampler;
import gb.esac.timeseries.TimeSeriesUtils;
import gb.esac.tools.BasicStats;
import gb.esac.tools.Complex;
import gb.esac.tools.ComplexNumbers;
import gb.esac.tools.FFT;
import gb.esac.tools.MinMax;
import gb.esac.tools.MyFFT;
import org.apache.log4j.Logger;


/**
 * <code>PeriodogramMaker</code> is a factory class that makes four basic different kinds of periodograms:
 <ul>
 <li><code>FFTPeriodogram</code> (<a href="http://adsabs.harvard.edu/abs/2002nrc..book.....P">Press 2002<a>), 
 <li><code>ModifiedRayleighPeriodogram</code> (<a href="http://adsabs.harvard.edu/abs/1996APh.....4..235O">Orford 1996, Astroparticle Physics, 4, 235<a>), 
 <li><code>ModifiedRayleighPeriodogram (kth order)</code> 
 <li><code>ModifiedZPeriodogram</code> 
 <li><code>LikelihoodPeriodogram</code> 
 <li><code>RayleighPeriodogram</code> (<a href="http://adsabs.harvard.edu/abs/1983ApJ...272..256L">Leahy et al. 1983, ApJ, 272, 256<a>) and 
 <li><code>ZPeriodogram</code> (<a href="http://adsabs.harvard.edu/abs/1983A%26A...128..245B">Buccheri et al. 1983, A&A, 128, 245<a>).
 </ul>
 For the <code>FFTPeriodogram</code>, there is a large number of factory methods that allow to make periodograms 
with or without applying a smoothing window, and with or without oversampling.
 *
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, European Space Astronomy Centre, Villanueva de la Canada (Madrid), Spain
 * @version (last modified) Dec 2015
 */
public final class PeriodogramMaker {
    private static Logger logger  = Logger.getLogger(PeriodogramMaker.class);

    //  Core factory method for FFTPeriodogram
    public static FFTPeriodogram makeUnnormalizedWindowedFFTPeriodogram(double[] countsOrRates, double duration, String windowName, int samplingFactor) throws PeriodogramException {
	logger.info("Making UnnormalizedWindowedFFTPeriodogram");
	logger.info("  Window function: "+windowName);
	logger.info("  Sampling factor: "+samplingFactor);
	logger.warn("  Treatment assumes uniform sampling with equal bins");
	if ( samplingFactor > 1 ) {
	    double powerOfTwo = Math.log(samplingFactor)/Math.log(2);
	    double integerPart = Math.floor(powerOfTwo);
	    double diff = powerOfTwo - integerPart;
	    if ( diff != 0 ) {
		throw new PeriodogramException("Sampling factor for FFTPeriodogram must be a power of 2: "+samplingFactor+" is not.");
	    }
	}
	//  Subtract mean	
 	int nDataBins = countsOrRates.length;
	double[] binHeights = new double[nDataBins];
    	double avgOfInputData = BasicStats.getMean(countsOrRates);
	for ( int i=0; i < nDataBins; i++ ) {
	    binHeights[i] = countsOrRates[i] - avgOfInputData;
	}
	logger.info("  Sum of intensities BEFORE mean-subtraction = "+BasicStats.getSum(countsOrRates)+" (after = "+BasicStats.getSum(binHeights)+")");
	logger.info("  Sum of squared intensities AFTER mean-subtraction = "+BasicStats.getSumOfSquares(binHeights));
	//  Apply the smoothing window
	WindowFunction windowFunction = null;
	try {
	    windowFunction = new WindowFunction(windowName);
	}
	catch ( WindowFunctionException e ) {
	    throw new PeriodogramException("Cannot construct window function ("+windowName+")", e);
	}
	double[] windowedData = windowFunction.apply(binHeights);
	// NO WINDOW
	windowedData = binHeights;
  	//  Define number of bins as a power-of-two
  	double n = Math.log(nDataBins)/Math.log(2);
  	double exp = Math.ceil(n);
  	int nPowerOfTwoBins = (int) Math.pow(2, exp);
	//  Keep the original number of bins
	//int nPowerOfTwoBins = nDataBins;
	//
 	//  Padding screws up the frequencies because the step is smaller,
	//  but it can be fixed: remains to be fixed
 	//
  	int nBins = samplingFactor*nPowerOfTwoBins;
  	double[] paddedData = padWithZeros(windowedData, nBins);
	// NO PADDING
	paddedData = windowedData;

 	//  Define test frequencies
 	double timeBinWidth = duration/nPowerOfTwoBins;
	double nuMin = 1d/duration;
	double nuMax = 1d/(2*timeBinWidth);
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, duration, samplingFactor);
	//  Do the FFT
	//  Using FFT.java
	logger.info("Calculating the FFT of input data");
   	Complex[] binHeightsForFFT = Complex.realDoubleArrayToComplexArray(paddedData);
   	Complex[] resultOfFFT = FFT.fft(binHeightsForFFT);
  	double[] power = Complex.normSquared(resultOfFFT);
	//  Using  MyFFT.java
 	// int nn = paddedData.length;
  	// double[] resultOfFFT = MyFFT.fft(ComplexNumbers.myComplex(paddedData), nn, +1);
  	// //double[] resultOfFFT = MyFFT.fft(ComplexNumbers.myComplex(smoothedData), nn, +1);
  	// double[] power = ComplexNumbers.getPower(resultOfFFT);
 	// //  Correct for different normalization compared to FFT.java
 	// for ( int i=0; i < power.length; i++ ) {
 	//     power[i] *= nn*nn;
 	// }
	//  Using  jTransform
//   	DoubleFFT_1D jtransformFFT = new DoubleFFT_1D(paddedData.length);
//   	jtransformFFT.realForward(paddedData);
//   	double[] power = new double[paddedData.length];
//   	for ( int i=0; i < paddedData.length/2; i++ ) {
//   	    power[i] = paddedData[2*i]*paddedData[2*i] + paddedData[2*i+1]*paddedData[2*i+1];
//   	}

	//  Drop first terms and second half of power spectrum corresponding to negative frequencies
	int size = power.length/2;
	double[] pow = new double[size];
	for ( int i=0; i < size; i++ ) {
	    pow[i] = power[i+samplingFactor];
	}
// 	int i=0;
// 	int j=0;
// 	while ( i < size ) {
// 	    while ( i <= samplingFactor ) {
// 		pow[j] = power[i+samplingFactor+1];
// 		i++;
// 		j++;
// 	    }
// 	    j++;
// 	}
	// Keep only the good physically meaningful frequencies
	DoubleArrayList goodFreqs = new DoubleArrayList();
	DoubleArrayList goodPowers = new DoubleArrayList();
	int i=0;
	double min = nuMin;
	double max = nuMax + testFreqs.length*Math.ulp(nuMax);
	while ( testFreqs[i] < min ) i++;
	while ( i < testFreqs.length && testFreqs[i] <= max ) {
	    if ( !Double.isNaN(pow[i]) ) {
		goodPowers.add(pow[i]);
		goodFreqs.add(testFreqs[i]);
	    }
	    i++;
	}
	goodPowers.trimToSize();
	goodFreqs.trimToSize();
	if ( goodPowers.size() == 0 ) {
	    throw new PeriodogramException("All power values are NaN: Cannot construct Periodogram");
	}
	//  Construct the un-normalised FFTPeriodogram
	return new FFTPeriodogram(goodFreqs.elements(), goodPowers.elements(), samplingFactor);
    }

    /**
     * <code>makeUnnormalizedWindowedRateFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static FFTPeriodogram makeUnnormalizedWindowedRateFFTPeriodogram(TimeSeries timeSeries, String windowName, int samplingFactor) throws PeriodogramException, BinningException {
	//  We must ensure uniform sampling and equal sized bins
	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(timeSeries);
	TimeSeries ts_resamp = TimeSeriesMaker.makeTimeSeries(ts);
	if ( ! ts.binWidthIsConstant() ) {
	    double minBinWidth = timeSeries.minBinWidth();
	    ts_resamp = TimeSeriesResampler.resample(ts, minBinWidth);
	    //  This is not a good strategy because it degrades frequency resolution
	    // int k=0;
	    // while ( ts_resamp.thereAreGaps() ) {
	    // 	k++;
	    // 	logger.warn("There are still gaps: resampling with larger bin width");
	    // 	double binWidth = minBinWidth + (k*0.1)*minBinWidth;
	    // 	ts_resamp = TimeSeriesResampler.resample(ts, binWidth);
	    // }
	}
	if ( ts_resamp.thereAreGaps() ) {
	    ts_resamp = TimeSeriesUtils.fillGaps(ts_resamp);
	}
	return makeUnnormalizedWindowedFFTPeriodogram(ts_resamp.getRates(), ts_resamp.duration(), windowName, samplingFactor);
    }

    /**
     * <code>makeUnnormalizedWindowedFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static FFTPeriodogram makeUnnormalizedWindowedFFTPeriodogram(TimeSeries timeSeries, String windowName, int samplingFactor) throws PeriodogramException, BinningException {
	//  We must ensure uniform sampling and equal sized bins
	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(timeSeries);
	if ( ! ts.binWidthIsConstant() ) {
	    double minBinWidth = timeSeries.minBinWidth();
	    ts = TimeSeriesResampler.resample(ts, minBinWidth);
	}
	return makeUnnormalizedWindowedFFTPeriodogram(ts.getBinHeights(), ts.duration(), windowName, samplingFactor);
    }

    //  Oversampled Windowed FFTPeriodogram

    /**
     * <code>makeOversampledWindowedFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledWindowedFFTPeriodogram(TimeSeries timeSeries, String windowName, String normName, int samplingFactor) throws PeriodogramException, BinningException  {
	FFTPeriodogram basicPeriodogram = makeUnnormalizedWindowedFFTPeriodogram(timeSeries, windowName, samplingFactor);
	return applyNormalization(basicPeriodogram, timeSeries, windowName, normName, "counts");
    }

    /**
     * <code>makeOversampledWindowedFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesFileException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static FFTPeriodogram makeOversampledWindowedFFTPeriodogram(EventList evlist, String windowName, String normName, int samplingFactor) throws TimeSeriesException, PeriodogramException, BinningException {
	return makeOversampledWindowedFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), windowName, normName, samplingFactor);
    }

    //  Windowed FFTPeriodogram

    /**
     * <code>makeWindowedFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeWindowedFFTPeriodogram(TimeSeries timeSeries, String windowName, String normName) throws PeriodogramException, BinningException {
	int sampling = 1;
	return makeOversampledWindowedFFTPeriodogram(timeSeries, windowName, normName, sampling);
    }

    /**
     * <code>makeWindowedFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesFileException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static FFTPeriodogram makeWindowedFFTPeriodogram(EventList evlist, String windowName, String normName) throws TimeSeriesException, PeriodogramException, BinningException {
	return makeWindowedFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), windowName, normName);
    }

    /**
     * <code>makeWindowedFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeWindowedFFTPeriodogram(TimeSeries timeSeries, String windowName) throws PeriodogramException, BinningException {
	String normName = "leahy";
	return makeWindowedFFTPeriodogram(timeSeries, windowName, normName);
    }

    /**
     * <code>makeWindowedFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param windowName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesFileException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static FFTPeriodogram makeWindowedFFTPeriodogram(EventList evlist, String windowName) throws TimeSeriesException, PeriodogramException {
	String normName = "leahy";
	return makeWindowedFFTPeriodogram(evlist, windowName);
    }

    //   Oversampled Plain FFTPeriodogram (using a Rectangular window)

    /**
     * <code>makeOversampledPlainFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param normName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledPlainFFTPeriodogram(TimeSeries timeSeries, String normName, int samplingFactor) throws PeriodogramException, BinningException {
	String windowName = "rectangular";
	return makeOversampledWindowedFFTPeriodogram(timeSeries, windowName, normName, samplingFactor);
    }

    /**
     * <code>makeOversampledPlainFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param normName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesFileException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledPlainFFTPeriodogram(EventList evlist, String normName, int samplingFactor) throws TimeSeriesException, PeriodogramException, BinningException {
	String windowName = "rectangular";
	return makeOversampledPlainFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), normName, samplingFactor);
    }

    /**
     * <code>makeOversampledPlainFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledPlainFFTPeriodogram(TimeSeries timeSeries, int samplingFactor) throws PeriodogramException, BinningException {
	String normName = "leahy";
	return makeOversampledPlainFFTPeriodogram(timeSeries, normName, samplingFactor);
    }

    /**
     * <code>makeOversampledPlainFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesFileException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledPlainFFTPeriodogram(EventList evlist, int samplingFactor) throws TimeSeriesException, PeriodogramException, BinningException {
	return makeOversampledPlainFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), samplingFactor);
    }

    //  Plain FFT Periodogram

    /**
     * <code>makePlainFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainFFTPeriodogram(TimeSeries timeSeries, String normName) throws PeriodogramException, BinningException {
	int sampling = 1;
	return makeOversampledPlainFFTPeriodogram(timeSeries, normName, sampling);
    }

    /**
     * <code>makePlainFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesFileException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainFFTPeriodogram(EventList evlist, String normName) throws TimeSeriesException, PeriodogramException, BinningException {
	return makePlainFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), normName);
    }

    /**
     * <code>makePlainFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainFFTPeriodogram(TimeSeries timeSeries) throws PeriodogramException, BinningException {
	String normName = "leahy";
	return makePlainFFTPeriodogram(timeSeries, normName);
    }

    /**
     * <code>makePlainFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainFFTPeriodogram(EventList evlist) throws TimeSeriesException, PeriodogramException, BinningException {
	return makePlainFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist));
    }

    //  Make windowed Rate FFTPeriodogram

    /**
     * <code>makeOversampledWindowedRateFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledWindowedRateFFTPeriodogram(TimeSeries timeSeries, String windowName, String normName, int samplingFactor) throws PeriodogramException, BinningException {
	FFTPeriodogram basic = makeUnnormalizedWindowedRateFFTPeriodogram(timeSeries, windowName, samplingFactor);
	return applyNormalization(basic, timeSeries, windowName, normName, "rates");
    }

    /**
     * <code>makeOversampledWindowedRateFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeOversampledWindowedRateFFTPeriodogram(EventList evlist, String windowName, String normName, int samplingFactor) throws TimeSeriesException, PeriodogramException, BinningException {
	return makeOversampledWindowedRateFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), windowName, normName, samplingFactor);
    }

    /**
     * <code>makeWindowedRateFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeWindowedRateFFTPeriodogram(TimeSeries timeSeries, String windowName, String normName) throws PeriodogramException, BinningException {
	int sampling = 1;
	return makeOversampledWindowedRateFFTPeriodogram(timeSeries, windowName, normName, sampling);
    }

    /**
     * <code>makeWindowedRateFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param windowName a <code>String</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makeWindowedRateFFTPeriodogram(EventList evlist, String windowName, String normName) throws TimeSeriesException, PeriodogramException, BinningException {
	int sampling = 1;
	return makeWindowedRateFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), windowName, normName);
    }

    /**
     * <code>makePlainRateFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainRateFFTPeriodogram(TimeSeries timeSeries, String normName) throws PeriodogramException, BinningException {
	String windowName = "rectangular";
	return makeWindowedRateFFTPeriodogram(timeSeries, windowName, normName);
    }

    /**
     * <code>makePlainRateFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param normName a <code>String</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainRateFFTPeriodogram(EventList evlist, String normName) throws TimeSeriesException, PeriodogramException, BinningException {
	String windowName = "rectangular";
	return makeWindowedRateFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), windowName, normName);
    }

    /**
     * <code>makePlainRateFFTPeriodogram</code>
     *
     * @param timeSeries a <code>TimeSeries</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainRateFFTPeriodogram(TimeSeries timeSeries) throws PeriodogramException, BinningException {
	String normName = "leahy";
	return makePlainRateFFTPeriodogram(timeSeries, normName);
    }

    /**
     * <code>makePlainRateFFTPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @return a <code>FFTPeriodogram</code> value
     * @exception TimeSeriesException if an error occurs
     * @exception PeriodogramException if an error occurs
     */
    public static  FFTPeriodogram makePlainRateFFTPeriodogram(EventList evlist) throws TimeSeriesException, PeriodogramException, BinningException {
	String normName = "leahy";
	return makePlainRateFFTPeriodogram(TimeSeriesMaker.makeTimeSeries(evlist), normName);
    }

    //  Factory methods for ModifiedRayleighPeriodogram (including the General Modified Rayleigh Periodogram)

    private static double[] checkNuMinAndNuMax(double duration, double dtMin, double nuMin, double nuMax) {
	//  Check nuMin
	double min = 1/duration;
	if ( nuMin < min ) {
	    logger.warn("  Specified nuMin < 1/duration: "+nuMin+" < "+min+". Resetting to 1/duration.");
	    nuMin = min;
	}
	//  Check nuMax
	double nyquistFrequency = 1/(2*dtMin);
	//nyquistFrequency -= Math.ulp(nyquistFrequency);
	double max = nyquistFrequency;
	if ( nuMax > max ) {
	    logger.warn("  Specified nuMax > Nyquist frequency = 1/(2*dtMin): "+nuMax+" > "+max+". Resetting to Nyquist frequency.");
	    nuMax = nyquistFrequency;
	}
	return new double[] {nuMin, nuMax};
    }

    private static int checkSamplingFactor(int samplingFactor) {
	int sampling = samplingFactor;
	if ( samplingFactor < 1 ) {
	    logger.warn("  Specified sampling factor < 1. Resetting to 1.");   
	    sampling = 1;
	}
	return sampling;
    }

    public static LikelihoodPeriodogram makeLikelihoodPeriodogram(Periodogram dataPeriodogram, double[] modelPowers) throws PeriodogramException {
	double[] freqs = dataPeriodogram.getFreqs();
	double[] dataPowers = dataPeriodogram.getPowers();
	if ( dataPowers.length != modelPowers.length ) {
	    throw new PeriodogramException("Unequal number of data and model power values");
	}
	ExponentialLikelihood expL = new ExponentialLikelihood();
	double[] inverseLikelihoods = new double[dataPowers.length];
	for ( int i=0; i < dataPowers.length; i++ ) {
	    double likelihood = expL.getLogLikelihood(modelPowers[i], dataPowers[i]);
	    inverseLikelihoods[i] = -likelihood;
	    //double likelihood = expL.getLikelihood(modelPowers[i], dataPowers[i]);
	    //inverseLikelihoods[i] = 1/likelihood;
	}
	return new LikelihoodPeriodogram(freqs, inverseLikelihoods, dataPeriodogram.samplingFactor());
    }

    public static LikelihoodPeriodogram makeLikelihoodPeriodogram(Periodogram dataPeriodogram, Periodogram modelPeriodogram) throws PeriodogramException {
	double[] dataFreqs = dataPeriodogram.getFreqs();
	double[] dataPowers = dataPeriodogram.getPowers();
	double[] modelFreqs = modelPeriodogram.getFreqs();
	double[] modelPowers = modelPeriodogram.getPowers();
	if ( dataPowers.length != modelPowers.length ) {
	    throw new PeriodogramException("Unequal number of data and model power values");
	}
	ExponentialLikelihood expL = new ExponentialLikelihood();
	double[] inverseLikelihoods = new double[dataPowers.length];
	for ( int i=0; i < dataPowers.length; i++ ) {
	    double likelihood = expL.getLogLikelihood(modelPowers[i], dataPowers[i]);
	    inverseLikelihoods[i] = -likelihood;
	    //double likelihood = expL.getLikelihood(modelPowers[i], dataPowers[i]);
	    //inverseLikelihoods[i] = 1/likelihood;
	}
	return new LikelihoodPeriodogram(dataFreqs, inverseLikelihoods, dataPeriodogram.samplingFactor());
    }

    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(EventList evlist) throws PeriodogramException {
	int harmonic = 1;
	return makeModifiedRayleighPeriodogram(evlist, harmonic);
    }

    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param harmonic an <code>int</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(EventList evlist, int harmonic) throws PeriodogramException {
	double duration = evlist.duration();
	double nuMin = 1/duration;
	double nyquistFrequency = 2*evlist.minEventSpacing();
	double effectiveNyquistFrequency = 2*evlist.meanRate();
	effectiveNyquistFrequency += Math.ulp(effectiveNyquistFrequency);
	double nuMax = effectiveNyquistFrequency;
	int samplingFactor = 1;
	return makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, samplingFactor, harmonic);
    }

    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(EventList evlist, double nuMin, double nuMax, int samplingFactor) throws PeriodogramException {
	int harmonic = 1;
	return makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, samplingFactor, harmonic);
    }

    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @param harmonic an <code>int</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(EventList evlist, double nuMin, double nuMax, int samplingFactor, int harmonic) throws PeriodogramException {
	logger.info("Making ModifiedRayleighPeriodogram from EventList");
	logger.info("  harmonic = "+harmonic);
	logger.info("  sampling = "+samplingFactor);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	logger.info("  nuMin = "+nuMin);
	logger.info("  nuMax = "+nuMax);
	double duration = evlist.duration();
	double dtMin = evlist.minEventSpacing();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	//  Define the test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	return makeModifiedRayleighPeriodogram(evlist, testFreqs, samplingFactor, harmonic);
    }

    public static ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(EventList evlist, double[] testFreqs, int samplingFactor) throws PeriodogramException {
	int harmonic = 1;
	return makeModifiedRayleighPeriodogram(evlist, testFreqs, samplingFactor, harmonic);
    }

    public static ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(EventList evlist, double[] testFreqs, int samplingFactor, int harmonic) throws PeriodogramException {
	double[] powers = calculateModifiedRayleighPowers(evlist, testFreqs, harmonic);
	return new ModifiedRayleighPeriodogram(testFreqs, powers, samplingFactor, harmonic);
    }

    public static double[] calculateModifiedRayleighPowers(EventList evlist, double[] testFreqs) {
	int harmonic = 1;
	return calculateModifiedRayleighPowers(evlist, testFreqs, harmonic);
    }

    public static double[] calculateModifiedRayleighPowers(EventList evlist, double[] testFreqs, int harmonic) {
	logger.info("Calculating modified Rayleigh powers");
	double[] arrivalTimes = evlist.getArrivalTimes();
	int nTrials = testFreqs.length;
	double[] powers = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {		
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getCorrectedPowerForThisHarmonic(arrivalTimes, period, harmonic);
	}
	return powers;
    }

    public static double[][] calculateModifiedRayleighPowerComponents(EventList evlist, double[] testFreqs, int harmonic) {
	double[] arrivalTimes = evlist.getArrivalTimes();
	int nTrials = testFreqs.length;
	double[] powers = new double[nTrials];
	double[] meansCos = new double[nTrials];
	double[] meansSin = new double[nTrials];
	double[] expectedMeansCos = new double[nTrials];
	double[] expectedMeansSin = new double[nTrials];
	double[] variancesOfCos = new double[nTrials];
	double[] variancesOfSin = new double[nTrials];
	double[] covariancesOfCosSin = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {		
	    double period = 1.0/testFreqs[i];
	    double[] components = PowerCalculator.getCorrectedPowerComponentsForThisHarmonic(arrivalTimes, period, harmonic);
	    powers[i] = components[0];
	    meansCos[i] = components[1];
	    meansSin[i] = components[2];
	    expectedMeansCos[i] = components[3];
	    expectedMeansSin[i] = components[4];
	    variancesOfCos[i] = components[5];
	    variancesOfSin[i] = components[6];
	    covariancesOfCosSin[i] = components[7];
	}
	return new double[][] {powers, meansCos, meansSin, expectedMeansCos, expectedMeansSin, variancesOfCos, variancesOfSin, covariancesOfCosSin};
    } 
    
    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     */
    public static  ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(TimeSeries lc, int samplingFactor) {
	double nuMin = 1/lc.duration();
	double dt = lc.minBinWidth();
	dt = lc.maxBinWidth();
	double nyquistFrequency = 1/(2*dt);
	return makeModifiedRayleighPeriodogram(lc, nuMin, nyquistFrequency, samplingFactor);
    }

    public static  ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(TimeSeries lc) {
	return makeModifiedRayleighPeriodogram(lc, 1);
    }


    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     */
    public static  ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(TimeSeries lc, double nuMin, double nuMax, int samplingFactor) {
	int harmonic = 1;
	return makeModifiedRayleighPeriodogram(lc, nuMin, nuMax, samplingFactor, harmonic);
    }

    /**
     * <code>makeModifiedRayleighPeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @param harmonic an <code>int</code> value
     * @return a <code>ModifiedRayleighPeriodogram</code> value
     */
    public static  ModifiedRayleighPeriodogram makeModifiedRayleighPeriodogram(TimeSeries lc, double nuMin, double nuMax, int samplingFactor, int harmonic) {
	logger.info("Making ModifiedRayleighPeriodogram from TimeSeries");
	logger.info("  harmonic (k) = "+harmonic);
	logger.info("  sampling = "+samplingFactor);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	logger.info("  nuMin = "+nuMin);
	logger.info("  nuMax = "+nuMax);
	double duration = lc.duration();
	double dtMin = lc.minBinWidth();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	//  Define test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Get rates
	double[] binCentres = lc.getBinCentres();
	double[] rates = lc.getRates();
	double[] errors = new double[rates.length];
	if ( lc.errorsAreSet() ) {
	    logger.info("Time series errors are set: using defined uncertainties on rates");
	    errors = lc.getErrorsOnRates();
	}
	else {
	    logger.info("Time series errors are not set: using unweighted power calculation (errors[i]=1.0)");
	    for ( int i=0; i < errors.length; i++ ) {
		errors[i] = 1.0;
	    }
	}
	double meanRate = lc.meanRate();
	//   Fill gaps
 	// double[] filledRates = DataUtils.fillDataGaps(rates);
 	// doublep[ filledErrors = DataUtils.fillDataGaps(errors);
 	// meanRate = BasicStats.getMean(filledRates);
	//   Calculate power
	logger.info("Calculating modified Rayleigh powers");
	double[] powers = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getModRayleighPower(binCentres, rates, errors, period, meanRate, harmonic);
	    //powers[i] = PowerCalculator.getModRayleighPower(binCentres, filledRates, filledErrors, period, meanRate, harmonic);
	}
	return new ModifiedRayleighPeriodogram(testFreqs, powers, samplingFactor, harmonic);
    }

    /**
     * <code>makeLombScarglePeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>LombScarglePeriodogram</code> value
     */
    public static LombScarglePeriodogram makeLombScarglePeriodogram(TimeSeries lc, double nuMin, double nuMax, int samplingFactor) {
	logger.info("Making LombScarglePeriodogram with sampling factor "+samplingFactor);
	//  Check min and max frequencies
	double duration = lc.duration();
	double dtMin = lc.minBinWidth();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	//  Define test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Get rates info
	double[] binCentres = lc.getBinCentres();
	double[] rates = lc.getRates();
	double meanRate = lc.meanRate();
	double varianceInMeanSubtractedRates = BasicStats.getVariance(lc.getMeanSubtractedRates());
	//  Calculate powers
	logger.info("Calculating Lomb powers");
	double[] powers = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getLombPower(binCentres, rates, period, meanRate, varianceInMeanSubtractedRates);
	}
	return new LombScarglePeriodogram(testFreqs, powers, samplingFactor);
    }

    public static LombScarglePeriodogram makeLombScarglePeriodogram(TimeSeries ts) {
	double nuMin = 1/ts.duration();
	double dt = MinMax.getMin(ts.getBinWidths());
	//System.out.println(dt+" s");
	double nyquistFrequency = 1/(2*dt);
	//System.out.println(nyquistFrequency+" Hz");
	nyquistFrequency -= Math.ulp(nyquistFrequency);
	int samplingFactor = 1;
	return makeLombScarglePeriodogram(ts, nuMin, nyquistFrequency, samplingFactor);
    }

    /**
     * <code>makeRayleighPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @return a <code>RayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static RayleighPeriodogram makeRayleighPeriodogram(EventList evlist) throws PeriodogramException {
	double duration = evlist.duration();
	double nuMin = 1/duration;
	double effectiveNyquistFrequency = 2*evlist.meanRate();
	effectiveNyquistFrequency -= Math.ulp(effectiveNyquistFrequency);
	int samplingFactor = 1;
	return makeRayleighPeriodogram(evlist, nuMin, effectiveNyquistFrequency, samplingFactor);
    }

    /**
     * <code>makeRayleighPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>RayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  RayleighPeriodogram makeRayleighPeriodogram(EventList evlist, double nuMin, double nuMax, int samplingFactor) throws PeriodogramException {
	logger.info("Making RayleighPeriodogram from EventList with sampling factor of "+samplingFactor);
	double duration = evlist.duration();
	double dtMin = evlist.minEventSpacing();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	//  Define the test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Calculate powers
	double[] arrivalTimes = evlist.getArrivalTimes();
	double[] powers = new double[nTrials];
	int nHarmonics = 1;
	for ( int i=0; i < nTrials; i++ ) {		
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getZ2stats(arrivalTimes, period, nHarmonics)[2];
	}
	return new RayleighPeriodogram(testFreqs, powers, samplingFactor);
    }

    /**
     * <code>makeRayleighPeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @return a <code>RayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  RayleighPeriodogram makeRayleighPeriodogram(TimeSeries lc) throws PeriodogramException {
	double nuMin = 1/lc.duration();
	double dt = MinMax.getMin(lc.getBinWidths());
	double nyquistFrequency = 1/(2*dt);
	nyquistFrequency -= Math.ulp(nyquistFrequency);
	int samplingFactor = 1;
	return makeRayleighPeriodogram(lc, nuMin, nyquistFrequency, samplingFactor);
    }

    /**
     * <code>makeRayleighPeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @return a <code>RayleighPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  RayleighPeriodogram makeRayleighPeriodogram(TimeSeries lc, double nuMin, double nuMax, int samplingFactor) throws PeriodogramException {
	logger.info("Making RayleighPeriodogram from TimeSeries with sampling factor of "+samplingFactor);
	double duration = lc.duration();
	double dtMin = lc.minBinWidth();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	//  Define test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Calculate powers
	double[] binCentres = lc.getBinCentres();
	double[] rates = lc.getRates();
	double[] errors = lc.getErrorsOnRates();
	double[] powers = new double[nTrials];
	int nHarmonics = 1;
	for ( int i=0; i < nTrials; i++ ) {
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getZ2stats(binCentres, rates, errors, period, nHarmonics)[2];
	}
	return new RayleighPeriodogram(testFreqs, powers, samplingFactor);
    }

    /**
     * <code>makeModifiedZPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @param nHarmonics an <code>int</code> value
     * @return a <code>ModifiedZPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  ModifiedZPeriodogram makeModifiedZPeriodogram(EventList evlist, double nuMin, double nuMax, int samplingFactor, int nHarmonics) throws PeriodogramException {
	double duration = evlist.duration();
	double dtMin = evlist.minEventSpacing();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	//  Define test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Calculate powers
	double[] arrivalTimes = evlist.getArrivalTimes();
	double[] powers = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getModifiedZ2Power(arrivalTimes, period, nHarmonics);
	}
	return new ModifiedZPeriodogram(testFreqs, powers, samplingFactor, nHarmonics);
    }

    /**
     * <code>makeZPeriodogram</code>
     *
     * @param evlist an <code>EventList</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @param nHarmonics an <code>int</code> value
     * @return a <code>ZPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  ZPeriodogram makeZPeriodogram(EventList evlist, double nuMin, double nuMax, int samplingFactor, int nHarmonics) throws PeriodogramException {
	double duration = evlist.duration();
	double dtMin = evlist.minEventSpacing();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	//  Define test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Calculate powers
	double[] arrivalTimes = evlist.getArrivalTimes();
	double[] powers = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getZ2stats(arrivalTimes, period, nHarmonics)[2];
	}
	return new ZPeriodogram(testFreqs, powers, samplingFactor, nHarmonics);
    }

    /**
     * <code>makeZPeriodogram</code>
     *
     * @param lc a <code>TimeSeries</code> value
     * @param nuMin a <code>double</code> value
     * @param nuMax a <code>double</code> value
     * @param samplingFactor an <code>int</code> value
     * @param nHarmonics an <code>int</code> value
     * @return a <code>ZPeriodogram</code> value
     * @exception PeriodogramException if an error occurs
     */
    public static  ZPeriodogram makeZPeriodogram(TimeSeries lc, double nuMin, double nuMax, int samplingFactor, int nHarmonics) throws PeriodogramException {
	double duration = lc.duration();
	double dtMin = lc.minBinWidth();
	double[] nuMinAndNuMax = checkNuMinAndNuMax(duration, dtMin, nuMin, nuMax);
	int sampling = checkSamplingFactor(samplingFactor);
	samplingFactor = sampling;
	//  Define test frequencies
	double[] testFreqs = PeriodogramUtils.getFourierFrequencies(nuMinAndNuMax[0], nuMinAndNuMax[1], duration, samplingFactor);
	int nTrials = testFreqs.length;
	//  Calculate powers
	double[] binCentres = lc.getBinCentres();
	double[] rates = lc.getRates();
	double[] errors = lc.getErrorsOnRates();
	double[] powers = new double[nTrials];
	for ( int i=0; i < nTrials; i++ ) {
	    double period = 1.0/testFreqs[i];
	    powers[i] = PowerCalculator.getZ2stats(binCentres, rates, errors, period, nHarmonics)[2];
	}
	return new ZPeriodogram(testFreqs, powers, samplingFactor, nHarmonics);
    }

    private static double[] padWithZeros(double[] data, int nBins) {
	double[] paddedData = new double[nBins];
	int nZeros = nBins - data.length;
	int nZerosBefore = (int) Math.floor(nZeros/2d);
	int nZerosAfter = nZeros - nZerosBefore;
	for ( int i=0; i < nZerosBefore; i++ ) {
	    paddedData[i] = 0;
	}
	for ( int i=0; i < data.length; i++ ) {
	    paddedData[i+nZerosBefore] = data[i];
	}
	for ( int i=0; i < nZerosAfter; i++ ) {
	    paddedData[i+nZerosBefore+data.length] = 0;
	}
	return paddedData;
    }

    private static String[] normNames = new String[] {"Leahy", "Miyamoto (rms^2)", "Variance", "Leahy-like"};
    public static void printAvailableNormalizations() {
	logger.info("Available normalizations are:");
	for ( int i=0; i < normNames.length; i++ ) {
	    logger.info("  "+normNames[i]);
	}
    }

    private static FFTPeriodogram applyNormalization(FFTPeriodogram basicPer, TimeSeries ts, String windowName, String normName, String inputDataType) throws PeriodogramException {
	logger.info("Applying '"+normName+"' normalization (for '"+inputDataType+"' data type)");
	//  Calculate normalization factors
	double duration = ts.duration();
	double freqBinWidth = basicPer.binWidth();
	double[] rawPowers = basicPer.getPowers();
	double sumOfRawPowers = BasicStats.getSum(rawPowers);
	double avgRawPower = sumOfRawPowers/rawPowers.length;
	int n = ts.nBins();
	double varNorm = 2d/(n*n*freqBinWidth);
	double leahyLikeNorm = 2d/avgRawPower;
	double sumOfSquaredIntensities = 0;
	if ( inputDataType.equals("counts") ) {
	    sumOfSquaredIntensities = BasicStats.getSumOfSquares(ts.getMeanSubtractedBinHeights());
	}
	else if ( inputDataType.equals("rates") ) {
	    sumOfSquaredIntensities = BasicStats.getSumOfSquares(ts.getMeanSubtractedRates());
	}
	else {
	    throw new PeriodogramException("Input data type must be: 'counts' or 'rates'");
	}
	double leahyNorm = 2d/sumOfSquaredIntensities;
	double rmsNorm = leahyNorm/ts.meanRate();


	//  Define the normalization
	double norm = 0;
	if ( normName.equalsIgnoreCase("leahy") ) {
	    norm = leahyNorm;
	}
	else if ( normName.equalsIgnoreCase("miyamoto") ) {
	    norm = rmsNorm;
	}
	else if ( normName.equalsIgnoreCase("variance") ) {
	    norm = varNorm;
	}
	else if ( normName.equalsIgnoreCase("leahy-like") ) {
	    norm = leahyLikeNorm;
	}
	else {
	    printAvailableNormalizations();
	    throw new PeriodogramException("Unknown normalization ("+normName+")");
	}
	//  Apply the normalization correction for the window function
	WindowFunction windowFunction = null;
	try { 
	    windowFunction = new WindowFunction(windowName);
	}
	catch ( WindowFunctionException e ) {
	    throw new PeriodogramException("Cannot construct window function", e);
	}
	double[] function = windowFunction.getFunction(ts.nBins());
 	double sumOfSquaredWeights = BasicStats.getSumOfSquares(function);
 	double windowNorm = function.length/sumOfSquaredWeights;
	norm *= windowNorm;
	//  Return the normalized FFTPeriodogram
	return (FFTPeriodogram) basicPer.scale(norm);
    }
}
