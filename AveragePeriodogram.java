package gb.esac.periodogram;

import gb.esac.io.AsciiDataFileWriter;
import java.io.IOException;
import java.util.Arrays;
import org.apache.log4j.Logger;


/**
 * Describe class <code>AveragePeriodogram</code> here.
 *
 * @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (Nov 2010, ESAC)
 */
public class AveragePeriodogram extends Periodogram {

    private static Logger logger  = Logger.getLogger(AveragePeriodogram.class);

    private double[] errors;

    private AveragePeriodogram() {
    
    }

    AveragePeriodogram(double[] freqs, double[] powers,  double[] errors) {

	setFreqs(freqs);
	setPowers(powers);
	setErrors(errors);
	setBinWidth(freqs[1] - freqs[0]);
    }

    //  These two methods relating to erros only exist in this sub-class of Periodogram
    private void setErrors(double[] errors) {
	this.errors = Arrays.copyOf(errors, errors.length);
    }

    public double[] getErrors() {
	return Arrays.copyOf(this.errors, this.errors.length);
    }


    //  These three must be implemented because they are abstract in Periodogram
    AveragePeriodogram modifyFreqs(double[] newFreqs) {
	return new AveragePeriodogram(newFreqs, this.powers, this.errors);
    }

    AveragePeriodogram modifyPowers(double[] newPowers) {
	return new AveragePeriodogram(this.freqs, newPowers, this.errors);
    }

    AveragePeriodogram modifyFreqsAndPowers(double[] newFreqs, double[] newPowers) {
	return new AveragePeriodogram(newFreqs, newPowers, this.errors);
    }


    //  These one are added in order to allow to modify the errors
    AveragePeriodogram modifyErrors(double[] newErrors) {
	return new AveragePeriodogram(this.freqs, this.powers, newErrors);
    }

    AveragePeriodogram modifyPowersAndErrors(double[] newPowers, double[] newErrors) {
	return new AveragePeriodogram(this.freqs, newPowers, newErrors);
    }


    
    /**
     * Method <code>writeAsQDP</code> writes the periodogram as ASCI in QDP format.
     *
     * @param filename a <code>String</code> value
     */
    public void writeAsQDP(String filename) {


	//  Define the header
	String[] header1 = new String[] {
		"DEV /XS",
		"READ SERR 2",
		"LAB T", "LAB F",
		"TIME OFF",
		"LINE STEP",
		"LOG ON",
		"LW 4", "CS 1.5",
		"LAB X Frequency (Hz)",
		"LAB Y Power",
		"VIEW 0.2 0.1 0.8 0.9",
		"SKIP SINGLE",
		"!"
	};

	String[] header2 = new String[] {
		"DEV /XS",
		"READ SERR 1 2",
		"LAB T", "LAB F",
		"TIME OFF",
		"LINE STEP",
		"LOG ON",
		"LW 4", "CS 1.5",
		"LAB X Frequency (Hz)",
		"LAB Y Power",
		"VIEW 0.2 0.1 0.8 0.9",
		"SKIP SINGLE",
		"!"
	};

	//  Compute the half error bar
	double[] halfErrors = new double[errors.length];
	for ( int i=0; i<errors.length; i++ ) {
	    halfErrors[i] = errors[i]/2d;
	}


	AsciiDataFileWriter out = null;	
	try {

	    out = new AsciiDataFileWriter(filename);
	    if ( this.binWidthIsConstant ) {
		out.writeData(header1, freqs, powers, halfErrors);
	    }
	    else {
		out.writeData(header2, freqs, halfBinWidths, powers, halfErrors);
	    }
	}
	catch ( IOException ex ) {

	    System.out.println("Cannot write to file "+filename+"\n"+ex);
	    System.exit(-1);
	}

	logger.info("Periodogram written as QDP to file "+filename);
	
    }


}
