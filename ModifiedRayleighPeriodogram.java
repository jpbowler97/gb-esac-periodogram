package gb.esac.periodogram;


import org.apache.log4j.Logger;


/**
 * Describe class <code>ModifiedRayleighPeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (Nov 2010, ESAC)
 * @version 2.0 (Oct 2013, ESAC) 
 * @version 3.0 (Dec 2013, ESAC) # Added harmonic definition and constructor
 */
public class ModifiedRayleighPeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(ModifiedRayleighPeriodogram.class);

    private int harmonic;

    private ModifiedRayleighPeriodogram() {
    
    }

    ModifiedRayleighPeriodogram(double[] freqs, double[] powers,  int samplingFactor) {
	this(freqs, powers, samplingFactor, 1);
    }

    ModifiedRayleighPeriodogram(double[] freqs, double[] powers,  int samplingFactor, int harmonic) {
	setFreqs(freqs);
	setPowers(powers);
	setBinWidth(freqs[1] - freqs[0]);
	setSamplingFactor(samplingFactor);
	this.harmonic = harmonic;
    }

    ModifiedRayleighPeriodogram modifyFreqs(double[] newFreqs) {
	return new ModifiedRayleighPeriodogram(newFreqs, this.powers, this.samplingFactor);
    }

    ModifiedRayleighPeriodogram modifyPowers(double[] newPowers) {
	return new ModifiedRayleighPeriodogram(this.freqs, newPowers, this.samplingFactor);
    }

    ModifiedRayleighPeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new ModifiedRayleighPeriodogram(newFreqs, newPowers, this.samplingFactor);
    }

    public int harmonic() {
	return this.harmonic;
    }


}
