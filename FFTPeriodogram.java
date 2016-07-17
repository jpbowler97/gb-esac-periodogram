package gb.esac.periodogram;

import org.apache.log4j.Logger;


/**
 * Describe class <code>FFTPeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, European Space Astronomy Centre, SRE-O, Villanueva de la Canada (Madrid), Spain
 * @version 1.0 (Aug 2015)
 */
public class FFTPeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(FFTPeriodogram.class);
    private String normalization;
    private int samplingFactor;

    /**
     * Creates a new <code>FFTPeriodogram</code> instance. 
     * This default empty constructor can only be called from within the class.
     *
     */
    private FFTPeriodogram() {

    }

    FFTPeriodogram(double[] freqs, double[] powers) {
	this(freqs, powers, 1);
    }

    FFTPeriodogram(double[] freqs, double[] powers, int samplingFactor) {
	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
	this.normalization = "Not normalized";
    }

    FFTPeriodogram(double[] freqs, double[] powers, String normalization) {
	this(freqs, powers, 1, normalization);
    }
    
    FFTPeriodogram(double[] freqs, double[] powers, int samplingFactor, String normalization) {
	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
	this.normalization = normalization;
    }

    FFTPeriodogram(double[] freqs, double[] binWidths, double[] powers, String normalization) {
	this(freqs, binWidths, powers, 1, normalization);
    }

    FFTPeriodogram(double[] freqs, double[] binWidths, double[] powers, int samplingFactor, String normalization) {
	this(freqs, powers, normalization);
	setBinWidths(binWidths);
	setSamplingFactor(samplingFactor);
    }


    FFTPeriodogram modifyFreqs(double[] newFreqs) {
	return new FFTPeriodogram(newFreqs, this.powers, this.normalization);
    }

    FFTPeriodogram modifyPowers(double[] newPowers) {
	return new FFTPeriodogram(this.freqs, newPowers, this.normalization);
    }

    FFTPeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new FFTPeriodogram(newFreqs, newPowers, this.normalization);
    }


    public String normalization() {
	return this.normalization;
    }


}
