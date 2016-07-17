package gb.esac.periodogram;


import org.apache.log4j.Logger;


/**
 * Describe class <code>LikelihoodPeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (Apr 2014, ESAC)
 */
public class LikelihoodPeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(LikelihoodPeriodogram.class);

    private int samplingFactor;

    private LikelihoodPeriodogram() {
    
    }

    LikelihoodPeriodogram(double[] freqs, double[] powers,  int samplingFactor) {

	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
    }

    LikelihoodPeriodogram modifyFreqs(double[] newFreqs) {
	return new LikelihoodPeriodogram(newFreqs, this.powers, this.samplingFactor);
    }

    LikelihoodPeriodogram modifyPowers(double[] newPowers) {
	return new LikelihoodPeriodogram(this.freqs, newPowers, this.samplingFactor);
    }

    LikelihoodPeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new LikelihoodPeriodogram(newFreqs, newPowers, this.samplingFactor);
    }


}
