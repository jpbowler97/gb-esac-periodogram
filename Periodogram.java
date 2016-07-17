package gb.esac.periodogram;

import java.io.IOException;
import java.util.Arrays;

import gb.esac.binner.BinningException;
import gb.esac.tools.BasicStats;
import gb.esac.tools.DataUtils;
import gb.esac.tools.MinMax;
import org.apache.log4j.Logger;


/**
 * Class <code>Periodogram</code> 
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC
 * @version 1.0  (last modified: Aug 2015)
 *
 */
public abstract class Periodogram {

    // Class variable
    private static Logger logger  = Logger.getLogger(Periodogram.class);

    //  Instance variables
    double nuMin;
    double nuMax;
    int nBins;
    double[] freqs;
    int samplingFactor;
    double[] powers;
    double avgPower;
    double maxPower;
    double freqAtMaxPower;
    double[] binWidths;
    double[] halfBinWidths;
    double binWidth;
    double[] binEdges;
    boolean binWidthIsConstant;

    /**
     * Method <code>setFreqs</code> defines the frequencies of the Periogram. 
     * This should be set third. 
     * Here are also defined the instance variables: nBins, nuMin, nuMax
     * @param freqs a <code>double[]</code> value
     */
    protected void setFreqs(double[] freqs) {
	this.freqs = Arrays.copyOf(freqs, freqs.length);
	this.nBins = freqs.length;
	this.nuMin = this.freqs[0];
	this.nuMax = this.freqs[this.nBins-1];
	logger.info("Periodogram has "+this.nBins+" frequency bins");
	logger.info("  nuMin = "+this.nuMin);
	logger.info("  nuMax = "+this.nuMax);
    }

    protected void setSamplingFactor(int samplingFactor) {
	this.samplingFactor = samplingFactor;
    }
    
    /**
     * Method <code>setBinWidth</code> defines a constant frequency binWidth. 
     * Usually Periodograms have constant binWidth. 
     * When the Periodogram is rebinned, then a logarithmic binning is often used, 
     * which leads to variable binWidhts. 
     * Here are also defined the halfBinWidths and binEdges, 
     * and the boolean instance variance this.binWidthIsConstant is set to true.
     * @param binWidth a <code>double</code> value
     */
    protected void setBinWidth(double binWidth) {
    	// defines arrays for binned data with constant bin width
	this.binWidth = binWidth;
	this.binWidthIsConstant = true;
	//  Define the binWidths, halfBinWidths and binEdges
	this.binWidths = new double[this.nBins];
	this.halfBinWidths = new double[this.nBins];
	this.binEdges = new double[2*this.nBins];
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binWidths[i] = this.binWidth;
	    this.halfBinWidths[i] = 0.5*this.binWidth;
	    this.binEdges[2*i] = this.freqs[i] - this.halfBinWidths[i];
	    this.binEdges[2*i+1] = this.freqs[i] + this.halfBinWidths[i];
	}
    }

    
    /**
     * Method <code>setBinWidths</code> defines a binWidth for each bin of the periodogram. 
     `Here are also defined the halfBinWidths and binEdges, and the boolean instance variable 
     this.binWidthIsConstant is set to true if the variance of the binwidths < 1e-6, 
     and set to false otherwise.
     * @param binWidths a <code>double[]</code> value
     */
    protected void setBinWidths(double[] binWidths) {
    	// defines arrays for binned data with non-constant bin width
	this.binWidths = Arrays.copyOf(binWidths, binWidths.length);
	//  Check if binning is constant
	double var = BasicStats.getVariance(binWidths);
	if ( var < 1e-6) {
	    binWidthIsConstant = true;
	    this.binWidth = binWidths[0];
	}
	else {
	    binWidthIsConstant = false;
	}
	//  set halfBinWidths and binEdges
	this.halfBinWidths = new double[nBins];
	this.binEdges = new double[2*nBins];	
	for ( int i=0; i < nBins; i++ ) {
	    this.halfBinWidths[i] = 0.5*this.binWidths[i];
	    this.binEdges[2*i] = freqs[i] - this.halfBinWidths[i];
	    this.binEdges[2*i+1] = freqs[i] + this.halfBinWidths[i];
	}
    }

    protected void setPowers(double[] powers) {
	this.powers = Arrays.copyOf(powers, powers.length);
	this.maxPower = MinMax.getMax(this.powers);
	this.freqAtMaxPower = this.freqs[DataUtils.getIndex(this.maxPower, this.powers)];
	this.avgPower = BasicStats.getMean(this.powers);
	logger.info("Avg power = "+this.avgPower);
	logger.info("Max power = "+this.maxPower+" (at f = "+this.freqAtMaxPower+")"); 
    }

    // ABSTRACT METHODS 
    abstract Periodogram modifyFreqs(double[] newFreqs);
    abstract Periodogram modifyPowers(double[] newPowers);
    abstract Periodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers);

    // PUBLIC
    public double nuMin() { return this.nuMin; }
    public double nuMax() { return this.nuMax; }
    public int nBins() { return this.nBins; }
    public double[] getFreqs() { return Arrays.copyOf(this.freqs, this.freqs.length); };
    public int samplingFactor() { return this.samplingFactor; }
    public double[] getPowers() { return Arrays.copyOf(this.powers, this.powers.length); }
    public double getMaxPower() { return this.maxPower; }
    public double getFreqAtMaxPower() { return this.freqAtMaxPower; }
    public double getAvgPower() { return this.avgPower; }
    public boolean binWidthIsConstant() { return this.binWidthIsConstant; }
    public double binWidth() throws PeriodogramException {
	if ( this.binWidthIsConstant ) {
	    return this.binWidth;
	}
	else {
	    throw new PeriodogramException("Bin width is not constant. Use getBinWidths() instead.");
	}
    }
    public double[] getBinWidths() { return Arrays.copyOf(this.binWidths, this.binWidths.length); }
    public double[] getHalfBinWidths() { return Arrays.copyOf(this.halfBinWidths, this.halfBinWidths.length); }
    public double[] getBinEdges() { return Arrays.copyOf(this.binEdges, this.binEdges.length); }

    public double getPowerAt(double frequency) {
	int k = DataUtils.getClosestIndexInSortedData(frequency, this.freqs);
	return this.powers[k];
    }

    public double getIntegratedPower(double nuStart, double nuStop) {
	return PeriodogramUtils.getIntegratedPower(this, nuStart, nuStop);
    }

    public double getIntegratedPower() {
	return getIntegratedPower(this.nuMin, this.nuMax);
    }

    public Periodogram add(double constant) {
    	// shifts periodogram vertically by constant
	double[] newPowers = new double[this.nBins()];
	for ( int i=0; i < this.nBins(); i++ ) {
	    newPowers[i] = this.powers[i] + constant;
	}
	return modifyPowers(newPowers);
    }

    public Periodogram add(Periodogram periodogram) throws PeriodogramException {
    	// adds the data from 2 periodograms (requires combatability)
	if ( this.nBins != periodogram.nBins() ) {
	    throw new PeriodogramException("Cannot combine periodograms: different number of bins ("+this.nBins+" != "+periodogram.nBins()+").");
	}
	logger.warn("This operation combining two periodograms assumes all frequencies are identical. If this is not true, use AggregatePeriodogram.");
	double[] p = periodogram.getPowers();
	double[] f = periodogram.getFreqs();
	double[] newPowers = new double[this.nBins];
	for ( int i=0; i < this.nBins; i++ ) {
	    newPowers[i] = p[i] + this.powers[i];
	}
	return modifyPowers(newPowers);
    }

    public Periodogram subtract(double constant) {
	double minusConstant = -constant;
	return add(minusConstant);
    }

    public Periodogram subtract(Periodogram periodogram) throws PeriodogramException {
    	// splits data into 2 periodograms (requires compatability)
	if ( this.nBins != periodogram.nBins() ) {
	    throw new PeriodogramException("Cannot combine periodograms: different number of bins ("+this.nBins+" != "+periodogram.nBins()+").");
	}
	logger.warn("This operation combining two periodograms assumes all frequencies are identical. If this is not true, use AggregatePeriodogram.");
	double[] p = periodogram.getPowers();
	double[] f = periodogram.getFreqs();
	double[] newPowers = new double[this.nBins];
	for ( int i=0; i < this.nBins; i++ ) {
	    newPowers[i] = this.powers[i] - p[i];
	}
	return modifyPowers(newPowers);
    }

    public Periodogram scale(double constant) {
    	// scales data by constant
	logger.info("Scaling periodogram by factor: "+constant);
	double[] newPowers = new double[this.nBins()];
	for ( int i=0; i < this.nBins(); i++ ) {
	    newPowers[i] = this.powers[i] * constant;
	}
	return modifyPowers(newPowers);	
    }

    public Periodogram scale(double[] scalingFactors) {
    	// scales data by scalingFactors
	logger.info("Scaling each power by different scaling factor");
	double[] newPowers = new double[this.nBins()];
	for ( int i=0; i < this.nBins(); i++ ) {
	    newPowers[i] = this.powers[i] * scalingFactors[i];
	}
	return modifyPowers(newPowers);
    }

    public Periodogram rebin(double rebinFactor, String binningType) throws BinningException {
 	return PeriodogramBinner.rebin(this, rebinFactor, binningType);
    }

    public Periodogram dropFirstFrequency() {
	double[] newFreqs = Arrays.copyOfRange(this.freqs, 1, freqs.length);
	double[] newPowers = Arrays.copyOfRange(this.powers, 1, powers.length);
	return modifyFreqsAndPowers(newFreqs, newPowers);
    }
    
    // writeAsQDP methods
    public void writeAsQDP(Periodogram[] psdArray, String filename) throws IOException {
	PeriodogramWriter.writeAsQDP(this, psdArray, filename);
    }

    public void writeAsQDP(String filename) {
	PeriodogramWriter.writeAsQDP(this, filename);
    }

    public void writeAsQDP(String[] header, String filename) {
	PeriodogramWriter.writeAsQDP(this, header, filename);
    }

    public void writeAsQDP(double[] function, String filename) {
	PeriodogramWriter.writeAsQDP(this, function, filename);
    }

    public void writeAsQDP(double[] function1, double[] function2, String filename) {
	PeriodogramWriter.writeAsQDP(this, function1, function2, filename);
    }    

    public void writeAsQDP(double[] func1, double[] func2, String dataLab, String lab1, String lab2, String panel, String filename) {
	PeriodogramWriter.writeAsQDP(this, func1, func2, lab1, lab2, dataLab, panel, filename);
    }

    public void writeAsQDP(double[] func1, double[] func2, double[] func3, String lab1, String lab2, String lab3, String filename) {
	PeriodogramWriter.writeAsQDP(this, func1, func2, func3, lab1, lab2, lab3, filename);
    }


}
