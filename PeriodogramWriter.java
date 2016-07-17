package gb.esac.periodogram;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import gb.esac.io.AsciiDataFileWriter;
import org.apache.log4j.Logger;


final class PeriodogramWriter {

    private static Logger logger  = Logger.getLogger(PeriodogramWriter.class);

    public static void writeAsQDP(Periodogram main, String[] header, String filename) {
	double[] freqs = main.getFreqs();
	double[] powers = main.getPowers();	
	AsciiDataFileWriter out = null;
	try {
	    out = new AsciiDataFileWriter(filename);
	    out.writeData(header, freqs, powers);
	}
	catch ( IOException ex ) {
	    logger.error("Cannot write to file "+filename+"\n"+ex);
	    System.exit(-1);
	}
	logger.info("Periodogram written (as QDP) to "+filename);
    }

    public static void writeAsQDP(Periodogram main, double[] function, String filename) {
	double[] freqs = main.getFreqs();
	double[] powers = main.getPowers();	
	AsciiDataFileWriter out = null;
	try {
	    out = new AsciiDataFileWriter(filename);
	    out.writeData(header1, freqs, powers, function);
	}
	catch ( IOException ex ) {
	    logger.error("Cannot write to file "+filename+"\n"+ex);
	    System.exit(-1);
	}
	logger.info("Periodogram written (as QDP) to "+filename);
    }

    public static void writeAsQDP(Periodogram main, double[] func1, double[] func2, String filename) {
	double[] freqs = main.getFreqs();
	double[] powers = main.getPowers();	
	AsciiDataFileWriter out = null;
	try {
	    out = new AsciiDataFileWriter(filename);
	    out.writeData(header1, freqs, powers, func1, func2);
	}
	catch ( IOException ex ) {
	    logger.error("Cannot write to file "+filename+"\n"+ex);
	    System.exit(-1);
	}
	logger.info("Periodogram written (as QDP) to "+filename);
    }

    public static void writeAsQDP(Periodogram main, Periodogram[] psdArray, String filename) throws IOException {
	double[] freqs = main.getFreqs();
	double[] powers = main.getPowers();
  	PrintWriter printWriter = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	for ( int i=0; i < header1.length; i++ ) {
	    printWriter.println(header1[i]);
	}
	for ( int i=0; i < freqs.length; i++ ) {
	    printWriter.println(freqs[i]+"\t"+powers[i]);
	}
	printWriter.println("NO NO");
	for ( int k=0; k < psdArray.length; k++ ) {
	    freqs = psdArray[k].getFreqs();
	    powers = psdArray[k].getPowers();
	    for ( int i=0; i < freqs.length; i++ ) {
		printWriter.println(freqs[i]+"\t"+powers[i]);
	    }
	    if ( k < psdArray.length ) {
		printWriter.println("NO NO");
	    }
	}
	printWriter.close();
	logger.info("Periodogram array written (as QDP) to "+filename);
    }

    public static void writeAsQDP(Periodogram main, double[] func1, double[] func2, String lab1, String lab2, String dataLab, String panel, String filename) {
	double[] freqs = main.getFreqs();
	double[] powers = main.getPowers();	
	String[] header = getHeaderWithDataPanelAndTwoFunctionLabels(dataLab, panel, lab1, lab2);
	AsciiDataFileWriter out = null;
	try {
	    out = new AsciiDataFileWriter(filename);
	    out.writeData(header, freqs, powers, func1, func2);
	}
	catch ( IOException ex ) {
	    logger.error("Cannot write to file "+filename+"\n"+ex);
	    System.exit(-1);
	}
	logger.info("Periodogram written (as QDP) to "+filename);
    }

    public static void writeAsQDP(Periodogram main, String filename) {
	if ( main.binWidthIsConstant() ) {
	    writeAsQDP(main, header1, filename);
	}
	else {
	    writeAsQDP(main, header2, filename);
	}
    }

    public static void writeAsQDP(Periodogram main, double[] func1, double[] func2, double[] func3, String lab1, String lab2, String lab3, String filename) {
	double[] freqs = main.getFreqs();
	double[] powers = main.getPowers();	
	String[] header = getHeaderForThreePanelPlot(lab1, lab2, lab3);
	AsciiDataFileWriter out = null;
	try {
	    out = new AsciiDataFileWriter(filename);
	    out.writeData(header, freqs, powers, func1, powers, func2, powers, func3);
	}
	catch ( IOException ex ) {
	    logger.error("Cannot write to file "+filename+"\n"+ex);
	    System.exit(-1);
	}
	logger.info("Periodogram written (as QDP) to "+filename);
    }


    //  Headers
    private static String[] header1 = 
	new String[] {
	"DEV /XS",
	"SKIP SINGLE",
	"CS 1.5","LW 4",
	"LAB T", "LAB F",
	"TIME OFF",
	"LINE ON",
	"LINE STEP ON 1",
	"VIEW 0.2 0.1 0.8 0.9",
	"LOG X ON",
	"LAB Y Power",
	"LAB X Frequency (Hz)",
	"!"
    };
    private static String[] header2 = 
	new String[] {
	"DEV /XS",
	"SKIP SINGLE",
	"READ SERR 1",
	"CS 1.5","LW 4",
	"LAB T", "LAB F",
	"TIME OFF",
	"LINE ON",
	"LINE STEP ON 1",
	"VIEW 0.2 0.1 0.8 0.9",
	"LOG X ON",
	"LAB Y Power",
	"LAB X Frequency (Hz)",
	"!"
    };
    private static String[] getHeaderWithDataPanelAndTwoFunctionLabels(String dataLab, String panel, String lab1, String lab2) {
	String[] header = new String[] {
	    "DEV /XS",
	    "CS 1.5", "LW 4", 
	    "TIME OFF",
	    "LAB T", "LAB F",
	    "LAB 1 VPOS 0.721 0.805 \"("+panel+")\" JUST RIGHT CS 1.5",
	    "LAB 2 VPOS 0.23 0.805 \""+dataLab+"\" JUST LEFT MA 45 CO 1 CS 1.5",
	    "LAB 3 VPOS 0.23 0.755 \""+lab1+"\" JUST LEFT MA 45 CO 2 CS 1.5",
	    "LAB 4 VPOS 0.23 0.705 \""+lab2+"\" JUST LEFT MA 45 CO 4 CS 1.5",
	    "LAB X Frequency (Hz)",
	    "LAB Y Power",
	    "LINE ON",
	    "LINE STEP ON 1",
	    "LOG X ON",
	    "CO 4 ON 3",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "!"
	};
	return header;
    }
    private static String[] getHeaderForThreePanelPlot(String lab1, String lab2, String lab3) {
	String[] header = new String[] {
	    "DEV /XS",
	    "CS 1.5", "LW 4", 
	    "PLOT VERT",
	    "TIME OFF",
	    "LOG X ON",
	    "LAB T", "LAB F",
	    "LAB 1 VPOS 0.23 0.85 \""+lab1+"\" JUST LEFT CO 1 CS 1.5",
	    "LAB 2 VPOS 0.23 0.584 \""+lab2+"\" JUST LEFT CO 1 CS 1.5",
	    "LAB 3 VPOS 0.23 0.317 \""+lab3+"\" JUST LEFT CO 1 CS 1.5",
	    "LAB 11 VPOS 0.77 0.85 \"(a)\" JUST RIGHT CS 1.5",
	    "LAB 12 VPOS 0.77 0.584 \"(b)\" JUST RIGHT CS 1.5",
	    "LAB 13 VPOS 0.77 0.317 \"(c)\" JUST RIGHT CS 1.5",
	    "CO 2 ON 2 4 6",
	    "CO 1 ON 3 5 7",
	    "LINE STEP ON 2 4 6",
	    "LINE ON 3 5 7",
	    "WIN 1", // top
	    "YPL 6 7",
	    "LOC 0 0.0667 1 0.4",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "LAB Y Power",
	    "WIN 2", // middle
	    "YPL 4 5",
	    "LOC 0 0.333 1 0.667",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "LAB Y Power",
	    "WIN 3", // top
	    "YPL 2 3",
	    "LOC 0 0.6 1 0.933", 
	    "VIEW 0.2 0.1 0.8 0.9",
	    "LAB Y Power",
	    "WIN ALL",
	    "LAB X Frequency (Hz)", 
	    "!"
	};
	String[] headerWithSpaceBtwPanels = new String[] {
	    "DEV /XS",
	    "CS 1.5", "LW 4", 
	    "TIME OFF",
	    "PLOT VERT",
	    "LOG X ON",
	    "LAB T", "LAB F",
	    "LAB 1 VPOS 0.23 0.83 \""+lab1+"\" JUST LEFT CO 1 CS 1.5",
	    "LAB 2 VPOS 0.23 0.58 \""+lab2+"\" JUST LEFT CO 1 CS 1.5",
	    "LAB 3 VPOS 0.23 0.33 \""+lab3+"\" JUST LEFT CO 1 CS 1.5",
	    "LAB 11 VPOS 0.77 0.83 \"(a)\" JUST RIGHT CS 1.5",
	    "LAB 12 VPOS 0.77 0.58 \"(b)\" JUST RIGHT CS 1.5",
	    "LAB 13 VPOS 0.77 0.33 \"(c)\" JUST RIGHT CS 1.5",
	    "CO 2 ON 2 4 6",
	    "CO 1 ON 3 5 7",
	    "LINE STEP ON 2 4 6",
	    "LINE ON 3 5 7",
	    "VIEW 0.2 0.1 0.8 0.9",
	    "WIN 1", // bottom
	    "LAB Y Power",
	    "YPL 6 7",
	    "LOC 0 0.1 1 0.395", 
	    "WIN 2", // middle
	    "LAB Y Power",
	    "YPLT 4 5",
	    "LOC 0 0.345 1 0.655",
	    "WIN 3", // top
	    "LAB Y Power",
	    "LOC 0 0.605 1 0.9",
	    "YPL 2 3",
	    "WIN ALL",
	    "LAB X Frequency (Hz)", 
	    "!"
	};
	return header;
    }


}