package gb.esac.periodogram;


import gb.esac.eventlist.EventList;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import java.util.Arrays;



public class TestAggregatePeriodogram {


	public static void main(String[] args) throws Exception {

	    double df = 1e-4;
	    AggregatePeriodogram avg = new AggregatePeriodogram();

// 	    double[] freqs = new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
// 	    double[] powers = new double[]{2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
// 	    ModifiedRayleighPeriodogram psd1 = new ModifiedRayleighPeriodogram(freqs, powers, 1);
// 	    avg.add(psd1);

// 	    freqs = new double[]{5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
// 	    powers = new double[]{3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
// 	    ModifiedRayleighPeriodogram psd2 = new ModifiedRayleighPeriodogram(freqs, powers, 1);
// 	    avg.add(psd2);
	    
	    int n=20;
	    double mean = 5;
	    double[] durations = new double[] {1e3, 1e4};
	    double duration = 1e3;
	    double alpha = 2;
	   
	    for ( int i=0; i < durations.length; i++ ) {
		duration = durations[i];		
		for ( int j=0; j < n; j++ ) {
		    EventList evlist = new EventList(RedNoiseGenerator.generateArrivalTimes(mean, duration, alpha));
		    double binTime = 1d;
		    TimeSeries ts = TimeSeriesMaker.makeTimeSeries(evlist, binTime);
		    //Periodogram p = PeriodogramMaker.makeModifiedRayleighPeriodogram(ts);
		    Periodogram p = PeriodogramMaker.makePlainFFTPeriodogram(ts, "leahy");
		    avg.add(p);
		}
	    }
	    //AveragePeriodogram avgPsd = avg.aggregate(Arrays.asList(new FFTPeriodogram[] {fft1, fft2}));
	    avg.getPeriodogram().writeAsQDP("avg.qdp");


	}


}
