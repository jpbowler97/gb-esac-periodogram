package gb.esac.periodogram;


import org.apache.log4j.Logger;


/**
 * Describe class <code>LombScarglePeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (May 2014, ESAC)
 */
public class LombScarglePeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(LombScarglePeriodogram.class);

    private LombScarglePeriodogram() {
    
    }

    LombScarglePeriodogram(double[] freqs, double[] powers,  int samplingFactor) {

	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
    }

    LombScarglePeriodogram modifyFreqs(double[] newFreqs) {
	return new LombScarglePeriodogram(newFreqs, this.powers, this.samplingFactor);
    }

    LombScarglePeriodogram modifyPowers(double[] newPowers) {
	return new LombScarglePeriodogram(this.freqs, newPowers, this.samplingFactor);
    }

    LombScarglePeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new LombScarglePeriodogram(newFreqs, newPowers, this.samplingFactor);
    }

}
