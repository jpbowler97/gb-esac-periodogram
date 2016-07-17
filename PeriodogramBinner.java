package gb.esac.periodogram;

import gb.esac.binner.BinningException;
import gb.esac.binner.BinningUtils;
import gb.esac.binner.Rebinner;
import gb.esac.binner.Resampler;
import org.apache.log4j.Logger;

final class PeriodogramBinner {

    private static Logger logger  = Logger.getLogger(PeriodogramBinner.class);

    static Periodogram rebin(Periodogram psd, double rebinFactor, String binningType) throws BinningException {
	if ( rebinFactor <= 1 ) {
	    throw new BinningException("The rebinning factor must be > 1");
	}
	logger.info("Rebinning periodogram");
	double xmin = psd.binEdges[0];
	double xmax = psd.binEdges[psd.binEdges.length-1];	
	logger.info("xmin = "+xmin);
	logger.info("xmax = "+xmax);
	//  Define oldBinWidth 
	int nOldBins = psd.nBins();
	double oldBinWidth = 0;
	logger.info("Number of bins = "+nOldBins);
	try {
	    oldBinWidth = psd.binWidth();
	}
	catch ( PeriodogramException e ) {
	    throw new BinningException("Cannot rebin Periodogram: bin width is not constant");
	}
	//  Define newBinWidth 
	double newBinWidth = oldBinWidth*rebinFactor;
	int nNewBins = (int) Math.floor(nOldBins*(oldBinWidth/newBinWidth));
	//  Do the rebinning
	double[] newFreqs = null;
	double[] newPowers = null;
	if ( binningType.equals("lin") ) {
	    logger.info("Number of new bins = "+nNewBins);
	    double[] newBinEdges = BinningUtils.getBinEdges(xmin, xmax, nNewBins);
	    newFreqs = BinningUtils.getBinCentresFromBinEdges(newBinEdges);
   	    newPowers = Resampler.resample(psd.getPowers(), psd.getBinEdges(), newBinEdges);
	}
	else if ( binningType.equals("log") ) {
	    double[] newBinEdges = BinningUtils.getBinEdgesInLogSpace(xmin, xmax, nNewBins);
	    newFreqs = BinningUtils.getBinCentresFromBinEdges(newBinEdges);
   	    newPowers = Resampler.resample(psd.getPowers(), psd.getBinEdges(), newBinEdges);
	}
	else if ( binningType.equals("papadakis") ) {
	    int nPoints = (int) Math.ceil(rebinFactor);
	    double[][] papadakis = Rebinner.rebinAsPapadakis(psd.getFreqs(), psd.getPowers(), nPoints);
	    newFreqs = papadakis[0];
	    newPowers = papadakis[1];
	}
	else {
	    throw new BinningException("Binning type must be: 'lin', 'log', or 'papadakis'");
	}
	return psd.modifyFreqsAndPowers(newFreqs, newPowers);
    }


}
