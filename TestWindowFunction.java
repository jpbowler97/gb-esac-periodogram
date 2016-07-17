package gb.esac.periodogram;

import gb.esac.eventlist.EventList;
import gb.esac.montecarlo.WhiteNoiseGenerator;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;

public class TestWindowFunction {


    public static void main(String[] args) throws Exception {

	double duration = 5e2;
	double mean = 10;
	double[] t = WhiteNoiseGenerator.generateArrivalTimes(mean, duration);
	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(new EventList(t));
	
	FFTPeriodogram psd = PeriodogramMaker.makePlainFFTPeriodogram(ts);
	psd.writeAsQDP("psd-smooth.qdp");
	psd = PeriodogramMaker.makePlainFFTPeriodogram(ts, "leahy");
	psd.writeAsQDP("psd.qdp");
	
    }

}
