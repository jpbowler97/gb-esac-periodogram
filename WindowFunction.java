package gb.esac.periodogram;

import cern.jet.math.Bessel;
import gb.esac.tools.BasicStats;
import java.util.Arrays;
import org.apache.log4j.Logger;

public class WindowFunction {

    //  See comparison table for different window functions at http://www.physik.uni-wuerzburg.de/~praktiku/Anleitung/Fremde/ANO14.pdf
    //  More details on Wikipedia at http://en.wikipedia.org/wiki/Window_function
    //  and comparison plot showing the leakage functions http://en.wikipedia.org/wiki/File:Window_function_(comparsion).png
    

    private static Logger logger  = Logger.getLogger(WindowFunction.class);

    private static String windowType;
    private static double[] windowFunction;
    private static String[] windowNames = new String[] {
	"Hann", "Hamming", "Cosine", "Lanczos", "Bartlett", "Gaussian", "Kaiser3", "Lanczos", "Bartlett-Hann", "Blackman", "Welch", "Parzen", "Rectangular"
    };


    //  Constructors
    private WindowFunction() {}
  
    public WindowFunction(String windowType) throws WindowFunctionException {
	setType(windowType);
    }


    //  Methods
    private static void setType(String type) throws WindowFunctionException {
	
	int i=0;
	try {
	    while ( !type.equalsIgnoreCase(windowNames[i]) ) {
		i++;
	    }
	}
	catch ( ArrayIndexOutOfBoundsException e ) {
	    printAvailableFunctions();
	    throw new WindowFunctionException("Unknown window type: "+type);
	}
	windowType = type;
    }

    public static void printAvailableFunctions() {

	logger.info("Available window functions are:");
	for ( int i=0; i < windowNames.length; i++ ) {
	    logger.info("  "+windowNames[i]);
	}
	logger.info("More info see Wikipedia at http://en.wikipedia.org/wiki/Window_function");
	logger.info("And comparison between windows at http://en.wikipedia.org/wiki/File:Window_function_(comparsion).png");
    }

    public static String[] getAvailableFunctions() {
	return Arrays.copyOf(windowNames, windowNames.length);
    }

    public static String getType() {
	return windowType;
    }

    public static double[] getFunction(int nBins) {
	computeFunction(nBins);
	return Arrays.copyOf(windowFunction, windowFunction.length);
    }

    public static double[] apply(double[] data) {
	computeFunction(data.length);
	double[] windowedData = new double[data.length];
	
	for ( int i=0; i < data.length; i++ ) {
	    windowedData[i] = data[i]*windowFunction[i];
	}
	return windowedData;
    }

    public static double[] apply(double[] data, double[] t, double duration) throws WindowFunctionException {
	if ( data.length != t.length ) {
	    throw new WindowFunctionException("Data array and time array must be of equal length.");
	}
	computeFunction(t, duration);
	double[] windowedData = new double[data.length];
	for ( int i=0; i < data.length; i++ ) {
	    windowedData[i] = data[i]*windowFunction[i];
	}
	return windowedData;
    }

    private static void computeFunction(int nBins) {

	windowFunction = new double[nBins];

	//  Compute
	double nMinusOne = nBins - 1;
	if ( windowType.equalsIgnoreCase("Hann") ) {
	    for ( int i=0; i < nBins; i++ ) {
		//windowFunction[i] = 0.5* (1 - Math.cos(2*Math.PI*i/nMinusOne));  //  standard 
		windowFunction[i] = 0.5* (1 - Math.cos(2*Math.PI*(i+1)/(nBins+1)));  //  better (see discussion on window functions in Wikipedia)
	    }	    
	}
	else if ( windowType.equalsIgnoreCase("Hamming") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 0.54 - 0.46*Math.cos(2*Math.PI*i/nMinusOne);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Cosine") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = Math.sin(Math.PI*i/nMinusOne);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Lanczos") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = Math.sin(2d*i/nMinusOne -1)/(2d*i/nMinusOne -1);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Bartlett") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 2d/nMinusOne * (nMinusOne/2d - Math.abs(i - nMinusOne/2d) );
	    }
	}
	else if ( windowType.equalsIgnoreCase("Gaussian") ) {
	    double sigma = 0.5;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = Math.exp(-0.5*( (i - nMinusOne/2d )/(sigma*nMinusOne/2d)));
	    }
	}
	else if ( windowType.equalsIgnoreCase("Kaiser3") ) {
	    double alpha = 3;
	    double alphaPI = alpha*Math.PI;
	    for ( int i=0; i < nBins; i++ ) {
		double argTop = alphaPI*Math.sqrt(1 - Math.pow(2*i/nMinusOne - 1, 2));
		windowFunction[i] = Bessel.i0(argTop)/Bessel.i0(alphaPI);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Lanczos" ) ) {
	    for ( int i=0; i < nBins; i++ ) {
		double arg = 2*i/nMinusOne - 1;
		windowFunction[i] = Math.sin(arg)/arg;
	    }
	}
	else if ( windowType.equalsIgnoreCase("Bartlett-Hann") ) {
	    double a0 = 0.62;
	    double a1 = 0.48;
	    double a2 = 0.38;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = a0 - a1*Math.abs(i/nMinusOne - 0.5) - a2*Math.cos(2*Math.PI*i/nMinusOne);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Blackman") ) {
	    double alpha = 0.16;
	    double a0 = 0.5*(1- alpha);
	    double a1 = 0.5;
	    double a2 = 0.5*alpha;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = a0 - a1*Math.cos(2*Math.PI*i/nMinusOne) + a2*Math.cos(4*Math.PI*i/nMinusOne);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Welch") ) {
	    for ( int i=0; i < nBins; i++ ) {
		double kernel = (i - 0.5*nMinusOne) / (0.5*(nBins+1));
		windowFunction[i] = 1 - kernel*kernel;
	    }
	}
	else if ( windowType.equalsIgnoreCase("Parzen") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 1 - Math.abs( (i - 0.5*nMinusOne)/(0.5*(nBins+1)) );
	    }
	}
	else if ( windowType.equalsIgnoreCase("Rectangular") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 1d;
	    }
	}
    }


    private static void computeFunction(double[] t, double duration) {

	int nBins = t.length;
	double nMinusOne = nBins - 1;
	windowFunction = new double[nBins];

	//  Compute
	if ( windowType.equalsIgnoreCase("Hann") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 0.5* (1 - Math.cos(2*Math.PI*t[i]/duration));  //  better (see discussion on window functions in Wikipedia)
	    }	    
	}
	else if ( windowType.equalsIgnoreCase("Hamming") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 0.54 - 0.46*Math.cos(2*Math.PI*t[i]/duration);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Cosine") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = Math.sin(Math.PI*t[i]/duration);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Lanczos") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = Math.sin(2d*t[i]/duration -1)/(2d*t[i]/duration -1);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Bartlett") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 2d/duration * (duration/2d - Math.abs(t[i] - duration/2d) );
	    }
	}
	else if ( windowType.equalsIgnoreCase("Gaussian") ) {
	    double sigma = 0.5;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = Math.exp(-0.5*( (t[i] - duration/2d )/(sigma*duration/2d)));
	    }
	}
	else if ( windowType.equalsIgnoreCase("Kaiser3") ) {
	    double alpha = 3;
	    double alphaPI = alpha*Math.PI;
	    for ( int i=0; i < nBins; i++ ) {
		double argTop = alphaPI*Math.sqrt(1 - Math.pow(2*t[i]/duration - 1, 2));
		windowFunction[i] = Bessel.i0(argTop)/Bessel.i0(alphaPI);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Lanczos" ) ) {
	    for ( int i=0; i < nBins; i++ ) {
		double arg = 2*t[i]/duration - 1;
		windowFunction[i] = Math.sin(arg)/arg;
	    }
	}
	else if ( windowType.equalsIgnoreCase("Bartlett-Hann") ) {
	    double a0 = 0.62;
	    double a1 = 0.48;
	    double a2 = 0.38;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = a0 - a1*Math.abs(i/duration - 0.5) - a2*Math.cos(2*Math.PI*t[i]/duration);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Blackman") ) {
	    double alpha = 0.16;
	    double a0 = 0.5*(1- alpha);
	    double a1 = 0.5;
	    double a2 = 0.5*alpha;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = a0 - a1*Math.cos(2*Math.PI*i/duration) + a2*Math.cos(4*Math.PI*t[i]/duration);
	    }
	}
	else if ( windowType.equalsIgnoreCase("Welch") ) {
	    double effectiveDeltaT = (t[t.length-1] - t[0])/t.length;
	    for ( int i=0; i < nBins; i++ ) {
		double kernel = (t[i] - 0.5*duration) / (0.5*(duration+effectiveDeltaT));
		windowFunction[i] = 1 - kernel*kernel;
	    }
	}
	else if ( windowType.equalsIgnoreCase("Parzen") ) {
	    double effectiveDeltaT = (t[t.length-1] - t[0])/t.length;
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 1 - Math.abs( (t[i] - 0.5*duration)/(0.5*(duration+effectiveDeltaT)) );
	    }
	}
	else if ( windowType.equalsIgnoreCase("Rectangular") ) {
	    for ( int i=0; i < nBins; i++ ) {
		windowFunction[i] = 1d;
	    }
	}
    }


}
