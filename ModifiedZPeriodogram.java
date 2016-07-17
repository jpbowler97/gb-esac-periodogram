package gb.esac.periodogram;


import org.apache.log4j.Logger;


/**
 * Describe class <code>ModifiedZPeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (Oct 2013, ESAC)
 */
public class ModifiedZPeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(ModifiedZPeriodogram.class);

    private int nHarmonics;

    private ModifiedZPeriodogram() {

    }

    ModifiedZPeriodogram(double[] freqs, double[] powers, int samplingFactor, int nHarmonics) {

	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
	this.nHarmonics = nHarmonics;
    }

    ModifiedZPeriodogram modifyFreqs(double[] newFreqs) {
	return new ModifiedZPeriodogram(newFreqs, this.powers, this.samplingFactor, this.nHarmonics);
    }

    ModifiedZPeriodogram modifyPowers(double[] newPowers) {
	return new ModifiedZPeriodogram(this.freqs, newPowers, this.samplingFactor, this.nHarmonics);
    }

    ModifiedZPeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new ModifiedZPeriodogram(newFreqs, newPowers,  this.samplingFactor, this.nHarmonics);
    }

    public int nHarmonics() {
	return this.nHarmonics;
    }
}
