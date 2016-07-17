package gb.esac.periodogram;


import org.apache.log4j.Logger;


/**
 * Describe class <code>ZPeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (Nov 2010, ESAC)
 */
public class ZPeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(ZPeriodogram.class);

    private int nHarmonics;

    private ZPeriodogram() {

    }

    ZPeriodogram(double[] freqs, double[] powers, int samplingFactor, int nHarmonics) {

	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
	this.nHarmonics = nHarmonics;
    }

    ZPeriodogram modifyFreqs(double[] newFreqs) {

	return new ZPeriodogram(newFreqs, this.powers, this.samplingFactor, this.nHarmonics);
    }

    ZPeriodogram modifyPowers(double[] newPowers) {

	return new ZPeriodogram(this.freqs, newPowers, this.samplingFactor, this.nHarmonics);
    }

    ZPeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {

	return new ZPeriodogram(newFreqs, newPowers,  this.samplingFactor, this.nHarmonics);
    }

    public int nHarmonics() {

	return this.nHarmonics;
    }
}
