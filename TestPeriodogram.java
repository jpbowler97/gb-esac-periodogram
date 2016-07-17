package gb.esac.periodogram;

import gb.esac.eventlist.EventList;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.montecarlo.WhiteNoiseGenerator;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesUtils;
import gb.esac.tools.DataUtils;
import gb.esac.tools.BasicStats;
import gb.esac.io.AsciiDataFileReader;

public class TestPeriodogram {


    public static void main(String[] args) throws Exception {

	String filename = "/Users/gbelanger/Documents/astroData/xmm/mkn421//rgslccorr/0560983301/RGS/Mkn421.0560983301.RGS.20.lc.fits";
	// filename = "/Users/gbelanger/Documents/astroData/xmm/mkn421//rgslccorr/0411083201/RGS/Mkn421.0411083201.RGS.20.lc.fits";
	// filename = "/Users/gbelanger/Documents/astroData/xmm/mkn421//rgslccorr/0411082501/RGS/Mkn421.0411082501.RGS.20.lc.fits";
 	//filename = "flare_clean.txt";
// 	filename = "/Users/gbelanger/Documents/astroData/xmm/1340/filt_m1_EVLI_SgrA.fits";
// 	filename = "/Users/gbelanger/Documents/astroData/xmm/1340/filt_m2_EVLI_SgrA.fits";
 	//filename = "/Users/gbelanger/Documents/astroData/xmm/1340/filt_pn_EVLI_SgrA.fits";
	//filename = "/Users/gbelanger/Documents/astroData/xmm/866/pn-flare-times.txt";
	//filename = "simEvlist.fits";

	//filename = "crab_pn.ds";
	//EventList evlist = new EventList(filename);

	// int nBins = 1024;
	//TimeSeries ts = TimeSeriesMaker.makeTimeSeries(evlist, nBins);
	//ts = TimeSeriesMaker.makeTimeSeries(ts.getBinCentres(), ts.getRates(), ts.getErrorsOnRates());
	
	//ts = TimeSeriesUtils.fillGaps(ts);
  	//FFTPeriodogram fft = PeriodogramMaker.makePlainFFTPeriodogram(ts);
	//ModifiedRayleighPeriodogram r2 = PeriodogramMaker.makeModifiedRayleighPeriodogram(ts);
	//fft.writeAsQDP(r2.getPowers(), "fft_r2.qdp");
	// double nuMin = 1/evlist.duration();
	// double nuMax = 1e-2;
	// nuMin = 29.7734;
	// nuMax = 29.77365;
	// int nHarm = 10;
	// int sampling = 41;
	//ModifiedZPeriodogram z2 = PeriodogramMaker.makeModifiedZPeriodogram(evlist, nuMin, nuMax, sampling, nHarm);
	//z2 = (ModifiedZPeriodogram) z2.scale(1d/nHarm);
	//double[] locAndWidth = DataUtils.getLocationAndFWHMOfPeak(z2.getFreqs(), z2.getPowers());
	//System.out.println("Peak location = "+locAndWidth[0]+"\t FWHM = "+locAndWidth[1]);

	double[] peakLocations = new double[5];
	double[] peakHeights = new double[5];
	double[] peakHWHMs = new double[5];

	//ModifiedRayleighPeriodogram r2 = PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, sampling, 1);	
	//double[] locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(r2.getFreqs(), r2.getPowers());
	AsciiDataFileReader in = new AsciiDataFileReader("psd_k1.qdp");
	double[] locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(in.getDblCol(0), in.getDblCol(1));
	System.out.println("Peak location = "+locAndWidth[0]+"\t Peak height = "+locAndWidth[1]+"\t HWHM = "+locAndWidth[2]/2);
	int i=0;
	peakLocations[i] = locAndWidth[0];
	peakHeights[i] = locAndWidth[1];
	peakHWHMs[i] = locAndWidth[2]/2;

	// ModifiedRayleighPeriodogram r2_2 = PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, sampling, 2);
	// locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(r2_2.getFreqs(), r2_2.getPowers());
	in = new AsciiDataFileReader("psd_k2.qdp");
	locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(in.getDblCol(0), in.getDblCol(1));
	System.out.println("Peak location = "+locAndWidth[0]+"\t Peak height = "+locAndWidth[1]+"\t HWHM = "+locAndWidth[2]/2);
	i++;
	peakLocations[i] = locAndWidth[0];
	peakHeights[i] = locAndWidth[1];
	peakHWHMs[i] = locAndWidth[2]/2;

	// ModifiedRayleighPeriodogram r2_3 = PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, sampling, 3);
	// locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(r2_3.getFreqs(), r2_3.getPowers());
	in = new AsciiDataFileReader("psd_k3.qdp");
	locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(in.getDblCol(0), in.getDblCol(1));
	System.out.println("Peak location = "+locAndWidth[0]+"\t Peak height = "+locAndWidth[1]+"\t HWHM = "+locAndWidth[2]/2);
	i++;
	peakLocations[i] = locAndWidth[0];
	peakHeights[i] = locAndWidth[1];
	peakHWHMs[i] = locAndWidth[2]/2;

	// ModifiedRayleighPeriodogram r2_4 = PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, sampling, 4);
	// locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(r2_4.getFreqs(), r2_4.getPowers());
	in = new AsciiDataFileReader("psd_k4.qdp");
	locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(in.getDblCol(0), in.getDblCol(1));
	System.out.println("Peak location = "+locAndWidth[0]+"\t Peak height = "+locAndWidth[1]+"\t HWHM = "+locAndWidth[2]/2);
	i++;
	peakLocations[i] = locAndWidth[0];
	peakHeights[i] = locAndWidth[1];
	peakHWHMs[i] = locAndWidth[2]/2;

	// ModifiedRayleighPeriodogram r2_5 = PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, sampling, 5);
	// locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(r2_5.getFreqs(), r2_5.getPowers());
	in = new AsciiDataFileReader("psd_k5.qdp");
	locAndWidth = DataUtils.getLocationHeightAndFWHMOfPeak(in.getDblCol(0), in.getDblCol(1));
	System.out.println("Peak location = "+locAndWidth[0]+"\t Peak height = "+locAndWidth[1]+"\t HWHM = "+locAndWidth[2]/2);
	i++;
	peakLocations[i] = locAndWidth[0];
	peakHeights[i] = locAndWidth[1];
	peakHWHMs[i] = locAndWidth[2]/2;
	double[] weights = new double[5];
	for ( i=0; i<5; i++ ) {
	    weights[i] = peakHeights[i]/Math.pow(peakHWHMs[i],2);
	}
	double wMeanFreq = BasicStats.getWMean_weights(peakLocations, weights);
	double uncertainty = BasicStats.getErrOnWMean(peakHWHMs);
	System.out.println("Resulting Weighted Mean is: "+wMeanFreq+" +/- "+ uncertainty);


	//ModifiedRayleighPeriodogram r2 = PeriodogramMaker.makeModifiedRayleighPeriodogram(ts, nuMin, nuMax, sampling, nHarm);
	//z2.writeAsQDP(r2_2.getPowers(), r2_3.getPowers(), "mod_z2_r2.qdp");
	//r2.writeAsQDP(r2_2.getPowers(), r2_3.getPowers(), "mod_r2.qdp");

	//double[] wz2 = PowerCalculator.getWeightedModifiedZ2Power(DataUtils.resetToZero(evlist.getArrivalTimes()), 0.0335868, nHarm);
	//System.out.println(wz2[0]+"\t"+wz2[1]);

	// filename = "evlist_pulsed.qdp";
	// evlist = new EventList(filename);
	// ts = TimeSeriesMaker.makeTimeSeries(evlist, nBins);
  	// fft = PeriodogramMaker.makePlainFFTPeriodogram(ts);
	// r2 = PeriodogramMaker.makeModifiedRayleighPeriodogram(ts);
	// fft.writeAsQDP(r2.getPowers(), "fft2.qdp");


	//fft = PeriodogramMaker.makeFFTPeriodogram(TimeSeriesUtils.kalmanFilter(ts));
	//fft.writeAsQDP("fft-kalman.qdp");

// 	double duration = ts.duration();
// 	double meanRate = ts.meanRate();
// 	int nBins = ts.nBins();
// 	double alpha = 2.6;
// 	EventList evlist = new EventList(RedNoiseGenerator.generateArrivalTimes(meanRate, duration, alpha));
// 	//EventList evlist = new EventList(RedNoiseGenerator.generateArrivalTimes(meanRate, duration, alpha, nBins));
// 	nBins = 128;
// 	ts = TimeSeriesMaker.makeTimeSeries(evlist, nBins);

// 	fft = PeriodogramMaker.makeFFTPeriodogram(ts);
// 	fft.writeAsQDP("fft-sim.qdp");

//  	fft = PeriodogramMaker.makeFFTPeriodogram(TimeSeriesUtils.kalmanFilter(ts));
//  	fft.writeAsQDP("fft-sim-kalman.qdp");

  	// double nuMin = 1e-4;
 	// double nuMax = 1e-2;
	//RayleighPeriodogram psd = PeriodogramMaker.makeRayleighPeriodogram(evlist, nuMin, nuMax, 1);
	//ZPeriodogram psd = PeriodogramMaker.makeZPeriodogram(evlist, nuMin, nuMax, 20, 5);
	//ModifiedZPeriodogram psd = PeriodogramMaker.makeModifiedZPeriodogram(evlist, nuMin, nuMax, 20, 5);
  	//ModifiedRayleighPeriodogram psd = PeriodogramMaker.makeModifiedRayleighPeriodogram(evlist, nuMin, nuMax, 20);
  	//ModifiedRayleighPeriodogram psd = PeriodogramMaker.makeModifiedRayleighPeriodogram(ts, nuMin, nuMax, 10);
  	//psd.writeAsQDP("z2.qdp");

// 	int nBins = ts.nBins();
// 	WindowFunction smoother = new WindowFunction("Blackman");
// 	double[] windowFunction = smoother.getFunction(nBins);
// 	double[] windowed_counts = smoother.apply(ts.getBinHeights());
// 	TimeSeries new_ts = TimeSeriesMaker.makeTimeSeries(ts.getBinEdges(), windowed_counts);
// 	ray = PeriodogramMaker.makeModifiedRayleighPeriodogram(new_ts, nuMin, nuMax, 10);
//   	ray.writeAsQDP("ray-psd-windowed.qdp");


    }

}
