package gb.esac.periodogram;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class AggregatePeriodogram {

    List<Frequency> frequencies = new ArrayList<Frequency>();

    Map<Frequency, AveragedPower> data = new HashMap<Frequency, AveragedPower>();
        
    public AggregatePeriodogram() {}
    
    public AggregatePeriodogram(final double df) {
        freqFactory = new FrequencyBinFactory() {
		@Override
		public Frequency createFrequency(double f) {
		    return new FrequencyBin(f, df);
		}
	    };
    }

    FrequencyBinFactory freqFactory = new FrequencyBinFactory() {
	    @Override
	    public Frequency createFrequency(double f) {
		return new Frequency(f);
	    }
	};

    double[] getFreqs() {
        return getPeriodogram().getFreqs();
    }
    
    double[] getPowers() {
        return getPeriodogram().getPowers();
    }
        
    public AveragePeriodogram getPeriodogram() {

        int n = frequencies.size();
        double freq[] = new double[n];
        double pow[] = new double[n];
        double error[] = new double[n];
        
        Collections.sort(frequencies, new Comparator<Frequency>() {
			     
			     @Override
			     public int compare(Frequency f1, Frequency f2) {
				 return Double.compare(f1.f0, f2.f0);
			     }
			 });
        
	int i = 0;
	for (Frequency f: frequencies) {
	    freq[i] = f.f0;
	    double[] avgNerr = data.get(f).avgPowerAndError(); 
	    pow[i] = avgNerr[0];
	    error[i] = avgNerr[1];
	    i++;
	}
	return new AveragePeriodogram(freq, pow, error);
    }
    
    public void add(Periodogram periodogram) {
	double[] freqs = periodogram.getFreqs();
	double[] powers = periodogram.getPowers();
	for (int i =0; i < freqs.length; i++) {
	    add(freqs[i], powers[i]);
	}
    }

    private void add(double freq, double power) {
	Frequency f = null;
	for (Frequency fi: frequencies) {
	    if (fi.covers(freq)) {
		f = fi;
		break;
	    }
	}
        
	AveragedPower avg = null;
	if (f == null) {
	    f = freqFactory.createFrequency(freq);
	    frequencies.add(f);
	    avg = new AveragedPower();
	    data.put(f,avg);
	}
	else {
	    f.add(freq);
	    avg = data.get(f);
	}
	avg.add(power);
    }
    
    
    private interface FrequencyBinFactory {
	Frequency createFrequency(double f); 
    }

    
    static class Frequency {

	protected double f0;
	private List<Double> f = new ArrayList<Double>(1);
        
	public Frequency(double freq) {
	    f0 = freq;
	    f.add(freq);
	}
        
	Frequency add(double freq){
	    f.add(freq);
	    return this;
	}
        
	boolean covers(double f) {
	    return Math.max(f0, f) - Math.min(f0, f) < 1E-4;
	}
    }
    
    static class FrequencyBin extends Frequency {

	double df;
        
	public FrequencyBin(double freq, double df) {
	    super(freq);
	    this.df = df;
	    f0 = computeIndex(freq);
	}
        
	@Override
	boolean covers(double f) {
	    double bi = computeIndex(f);
	    return f0<= bi && (bi < f0+df);
	}
        
	private double computeIndex(double f) {
	    return f - (f % df);
	}
    }
    
    static class AveragedPower {

	List<Double> powers = new LinkedList<Double>();

	AveragedPower add(double power){
	    powers.add(power);
	    return this;
	}

	private double[] avgPowerAndError(){
	    
	    double[] runningAvgAndVar = getRunningAvgAndVar(powers);
	    int n = powers.size();
	    // based on the exponential pdf: 
	    //   mu=avg, sigma=avg, so error on mean is avg/sqrt(n)
	    double avg = runningAvgAndVar[0];
	    double err = runningAvgAndVar[0]/Math.sqrt(n);
	    return new double[]{avg, err};
	}
        
	double getPower(){
	    return avgPowerAndError()[0];
	}
        
	double getError(){
	    return avgPowerAndError()[1];
	}
    }

    private static double[] getRunningAvgAndVar(List<Double> data) {

        double sum = data.get(0);
	
        int n = 1;
        double ave = sum/n;
        double var = Math.pow((data.get(0) - ave), 2);
	
        double runAve = ave;
        double runVar = 0;

	if ( data.size() > 1 ) {
	    for ( int i=1;  i < data.size(); i++ ) {
		double datum = data.get(i);
		if ( ! Double.isNaN(datum) ) {
		    runAve = ave + (datum - ave)/(n+1);
		    runVar = var + (datum - runAve)*(datum - ave);
		    sum += datum;
		    n++;
		    ave = sum/n;
		    var += Math.pow((datum - ave),2);
		    
		}
	    }
	    runVar /= (n-1);
	}
	
        return new double[] {runAve, runVar};
    }

    public static AveragePeriodogram aggregate(Collection<? extends Periodogram> c) {
        double df = Double.MAX_VALUE;
        for (Periodogram pg: c) {
	    df = Math.min(df,pg.getFreqs()[0]);
        }
        
        AggregatePeriodogram agr = new AggregatePeriodogram(df);
        for (Periodogram pg: c) {
	    agr.add(pg);
        }
        
        return agr.getPeriodogram();
    }
}

