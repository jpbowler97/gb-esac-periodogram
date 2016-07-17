package gb.esac.periodogram;

import java.text.DecimalFormat;
import org.apache.log4j.Logger;


public class TestPeriodogramUtils {

    private static Logger logger  = Logger.getLogger(TestPeriodogramUtils.class);

    private static DecimalFormat exp = new DecimalFormat("0.0000#E00");
    private static DecimalFormat index = new DecimalFormat("0000");

    public static void main(String[] args) throws PeriodogramException {

	// System.out.println(Math.ulp(0));
	// System.out.println(1e13+0.75*Math.ulp(1e13));
	// System.out.println(1e13+Math.ulp(1e13));
	// System.exit(0);

	int nbins = 1024;
	double duration = 9998;

	//  Test method getFFTFrequencies
	logger.info("Testing static method getFFTFrequencies(nbins="+nbins+", duration="+duration+")");
	double[] fftFreqs = PeriodogramUtils.getFFTFrequencies(nbins, duration);
	assert (fftFreqs.length == nbins/2);
	assert (fftFreqs[0] == 1/duration);
	assert (fftFreqs[fftFreqs.length-1] == nbins/(2*duration));
	int i=0; 
	while ( i < fftFreqs.length ) {
	    logger.info(index.format(i)+"\t"+exp.format(fftFreqs[i]));
	    i++;
	}
	logger.info("Test complete\n");
	
	//  Test method getFourierFrequencies with oversampling
	double nuMin = 1/duration;
	double dt = duration/nbins;
	double nuMax = 1/(2*dt);
	double[] freqs = null;
	logger.info("Testing static method getFourierFrequencies(nuMin="+nuMin+", nuMax="+nuMax+", duration="+duration+", oversampling)");
	for ( int j=1; j <= 1; j++ ) {
	    int oversampling = j;
	    logger.info("Oversampling factor = "+oversampling);
	    freqs = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, duration, oversampling);
	    assert (freqs.length == (fftFreqs.length*oversampling));
	    assert (freqs[0] == 1/duration);
	    double df = 1/(duration*oversampling);
	    assert (freqs[freqs.length-1] == nbins/(2*duration) + (oversampling-1)*df);
	    i=0; 
	    while ( i < freqs.length ) {
		System.out.println(index.format(i)+"\t"+exp.format(freqs[i]));
		i++;
	    }
	    logger.info("Done oversampling = "+oversampling);
	}
	logger.info("Test complete\n");

	System.exit(-1);
	
	//  Test method getFourierTestPeriods with oversampling
	double pMin = 1/nuMax;
	double pMax = 1/nuMin;
	double[] periods = null;
	logger.info("Testing static method getFourierPeriods(pMin="+pMin+", pMax="+pMax+", duration="+duration+", oversampling)");
	for ( int j=1; j <= 3; j++ ) {
	    int oversampling = j;
	    logger.info("Oversampling factor = "+oversampling);
	    periods = PeriodogramUtils.getFourierPeriods(pMin, pMax, duration, oversampling);
	    assert (periods.length == freqs.length);
	    assert (periods[0] == 1/freqs[freqs.length-1]);
	    assert (periods[periods.length-1] == 1/freqs[0]);
	    i=0; 
	    while ( i < periods.length ) {
		System.out.println(index.format(i)+"\t"+exp.format(periods[i]));
		i++;
	    }
	    logger.info("Done oversampling = "+oversampling);
	}
	logger.info("Test complete");


    }

}
