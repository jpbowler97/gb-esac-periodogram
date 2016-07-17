package gb.esac.periodogram;


import org.apache.log4j.Logger;


/**
 * Describe class <code>RayleighPeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (Nov 2010, ESAC)
 * @version 2.0 (Oct 2013, ESAC)
 */
public class RayleighPeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(RayleighPeriodogram.class);

    private RayleighPeriodogram() {

    }

    RayleighPeriodogram(double[] freqs, double[] powers, int samplingFactor) {

	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);

    }

    RayleighPeriodogram modifyFreqs(double[] newFreqs) {
	return new RayleighPeriodogram(newFreqs, this.powers, this.samplingFactor);
    }

    RayleighPeriodogram modifyPowers(double[] newPowers) {
	return new RayleighPeriodogram(this.freqs, newPowers, this.samplingFactor);
    }

    RayleighPeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new RayleighPeriodogram(newFreqs, newPowers, this.samplingFactor);
    }

}
