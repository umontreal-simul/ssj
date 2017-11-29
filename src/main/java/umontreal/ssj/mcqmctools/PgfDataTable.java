package umontreal.ssj.mcqmctools;

import java.util.Arrays;
import java.util.ArrayList;

/**
 * Represents a data table that could be used to produce LaTeX code to draw two-dimensional
 * plots with the pgfplot package.  Each table has a name, a number of observations (rows),
 * a number of fields (columns), an array that contains the names (identifiers) of the fields, 
 * and a two-dimensional array that contains the data.  Methods are available to format the table,
 * or selected columns of the table, or to produce LaTeX code to draw a plot of one field against
 * another field in pgfplot.
 */

public class PgfDataTable {

	String tableName;     // The name of the table, as a character string.
	String tableLabel;    // A label (short name) that might be used to identify the curve in plots.
	String[] fields;      // The names of the fields.  Should be short.
	double[][] data;      // The data points. data[s][j] for observation s, field j.
	int numFields;        // Number of fields for each data point.
	int numObservations;  // Number of observations (data points).
	
	/**
	 * Constructs a table with the given table name, array of field names, and data observations.
	 * 
	 * @param tableName    Name (descriptor) of the table.
	 * @param tableLabel   A label (very short name) for the table, to be used in plots to identify curves.
	 * @param fields       The names of the fields (the variables).
	 * @param data         The data points; the array @ref data[s] contains observation s.
	 * @return the table as a string.
	 */
	public PgfDataTable (String tableName, String tableLabel, String[] fields, double[][] data) {
		super();
		this.tableName = tableName;
		this.tableLabel = tableLabel;
		this.fields = fields;
		this.data = data;
		numFields = fields.length;
		numObservations = data[0].length;
	}

	/**
	 * Formats the full table as a String, with the field names in the first row, 
	 * and the data points (observations) in the following rows.
     * This format can be used directly by the pgfplot package.
	 * 
	 * @return the table as a string.
	 */
	public String formatTable() {
		StringBuffer sb = new StringBuffer("");
		sb.append("%      " + tableName + "\n      ");
		for (int j = 0; j < numFields; j++)  // For each field j
	       sb.append(fields[j] + "  ");
		sb.append("\n");
		for (int s = 0; s < numObservations; s++) { // For each cardinality n
			for (int j = 0; j < numFields; j++)  // For each field j
			   sb.append("      " + data[s][j]);
		    sb.append("\n");
		}
		return sb.toString();
	}

	/**
	 * Similar to @ref formatTable, but retains only two selected columns of the table 
	 * (i.e., two selected fields), specified by j1 and j2.  Note that the fields are numbered from 0.
	 * 
	 * @param j1 first selected field.
	 * @param j2 second selected field.
	 * @return the table as a string.
	 */
	public String formatTableTwoFields (int j1, int j2) {
		StringBuffer sb = new StringBuffer("");
		// sb.append("%     " + tableName + "\n");
	    sb.append("      " + fields[j1] + "  " + fields[j2] + " \n");
		for (int s = 0; s < numObservations; s++)  // For each cardinality n
		    sb.append("      " + data[s][j1] + "  " + data[s][j2] + " \n");
		return sb.toString();
	}

	/**
	 * Similar to @ref formatTableTwoFields, but outputs complete LaTeX code in an appropriate format
	 * that adds a curve of the field j2 against field j1 with the pgfplot package.
	 * 
	 * @param j1 first selected field.
	 * @param j2 second selected field.
	 * @param  plotoptions  is used to specify the options of addplot: <tt>addplot[plotoptions]</tt>
	 * @return LaTeX code as a string.
	 */
	public String formatPgfCurveAddPlot (int j1, int j2, String plotoptions) {
		StringBuffer sb = new StringBuffer("");
		sb.append("      \\addplot "   //  [" + plotoptions + "] 
			 + "table[x=" + fields[j1] + ",y=" + fields[j2] + "] { \n");
		sb.append( formatTableTwoFields (j1, j2) + " }; \n");
		sb.append("      \\addlegendentry{" + tableLabel + "}\n"); 
		sb.append("      \\label{" + tableLabel + "}\n"); 
		sb.append("%   \n");
		return sb.toString();
	}


	public String formatPgfCurveAddPlot (String xaxis, String yaxis, String plotoptions) {
		int j1 = Arrays.asList(fields).indexOf(xaxis);
		int j2 = Arrays.asList(fields).indexOf(yaxis);
		return formatPgfCurveAddPlot (j1, j2, plotoptions);
	}


	/**
	 * Returns a string that contains a complete tikzpicture for the pgfplot package,
	 * showing the field j2 against field j1.  
	 * 
	 * @param  title  is the title of the plot.
	 * @param j1 first selected field.
	 * @param j2 second selected field.
	 * @param  loglog  says if we want a log-log plot (true) or not (false).
	 * @param  plotoptions  is used to specify the options of addplot: <tt>addplot[plotoptions]</tt>
	 * @return LaTeX code as a string.
	 */
	public String drawPgfPlotSingleCurve (String title, String axistype, String xaxis, String yaxis, 
			 int logbasis, String axisoptions, String plotoptions) {
		int j1 = Arrays.asList(fields).indexOf(xaxis);
		int j2 = Arrays.asList(fields).indexOf(yaxis);
        return drawPgfPlotSingleCurve (title, axistype, j1, j2, 
 			   logbasis, axisoptions, plotoptions);
	}

	public String drawPgfPlotSingleCurve (String title, String axistype, int j1, int j2, 
			   int logbasis, String axisoptions, String plotoptions) {
		StringBuffer sb = new StringBuffer("");
		sb.append("  \\begin{tikzpicture} \n");
		sb.append("    \\begin{" + axistype + "}[ \n");
		sb.append("      title ={" + title  + "},\n");
		sb.append("      xlabel=" + fields[j1] + ",\n");
		sb.append("      ylabel=" + fields[j2] + ",\n");
		if (axistype == "loglogaxis") 
			sb.append("     log basis x=" + logbasis + ", log basis y=" + logbasis + ",\n");
		sb.append("       grid,\n");
		sb.append( axisoptions + ",\n");
		sb.append("      ] \n");
		sb.append( formatPgfCurveAddPlot (j1, j2, plotoptions));
		sb.append("    \\end{" + axistype + "}\n"); 
		sb.append("  \\end{tikzpicture}\n"); 
		sb.append("  "); 
		return sb.toString();
	}

	
	/**
	 * Returns a string that contains a complete tikzpicture for the pgfplot package,
	 * showing the field j2 against field j1, for all the curves that belong to listCurves.
	 * 
	 * 
	 * @param  title  is the title of the plot.
	 * @param j1 first selected field.
	 * @param j2 second selected field.
	 * @param  listCurves there will be one curve for each table in this list.
	 * @param  loglog  says if we want a log-log plot (true) or not (false).
	 * @param  plotoptions  is used to specify the options of addplot: <tt>addplot[plotoptions]</tt>
	 * @return LaTeX code as a string.
	 */
	public static String drawPgfPlotManyCurves (String title, String axistype, int j1, int j2, 
			ArrayList<PgfDataTable> listCurves, int logbasis, String axisoptions, String plotoptions) {
		StringBuffer sb = new StringBuffer("");
		sb.append("  \\begin{tikzpicture} \n");
		sb.append("    \\begin{" + axistype + "}[ \n");
		sb.append("      title ={" + title  + "},\n");
		sb.append("      xlabel=" + listCurves.get(0).fields[j1] + ",\n");
		sb.append("      ylabel=" + listCurves.get(0).fields[j2] + ",\n");
		if (axistype == "loglogaxis") 
			sb.append("      log basis x=" + logbasis + ", log basis y=" + logbasis + ",\n");
		sb.append("       grid,\n");
		sb.append( axisoptions + ",\n");
		// sb.append("      legend style={xshift=-3em,yshift=-2em}\n");
		sb.append("      ] \n");
		for (PgfDataTable curve : listCurves) {
			sb.append(curve.formatPgfCurveAddPlot (j1, j2, plotoptions));
		}
		sb.append("    \\end{" + axistype + "}\n"); 
		sb.append("  \\end{tikzpicture}\n"); 
		sb.append("  "); 
		return sb.toString();
	}
	
	
	public static String pgfplotFileHeader () {
		StringBuffer sb = new StringBuffer("");
		sb.append("\\documentclass[12pt]{article}\n" + 
         "\\usepackage{pgfplots}\n" + 
         "\\pgfplotsset{compat=1.12}\n" + 
         "\\usepackage{xcolor}\n" + 
         "\\usetikzlibrary{shapes,decorations,arrows,automata,plotmarks,patterns}\n" + 
         "\\usepackage{tikz-inet}\n" + 
         "\\usepgfplotslibrary{groupplots}\n" + 
         "%\n" + 
         "\\pgfplotscreateplotcyclelist{defaultcolorlist}{%\n" + 
         "  {blue!85!black,line width=0.9pt,mark=*,solid},\n" + 
         "  {red!85!black,line width=0.9pt,mark=square*,dashed},\n" + 
         "  {lime!60!black,line width=0.9pt,mark=triangle*,dotted},\n" + 
         "  {orange!90!black,line width=0.9pt,mark=diamond*,loosely dashed},\n" + 
         "  {black,mark=star,loosely dotted},\n" + 
         "  {violet!85!black,mark=o,solid},\n" + 
         "  {brown!85!black,mark=square,dashed},\n" + 
         "  {purple!85!black,mark=triangle,dotted}}\n" + 
         "%\n" + 
         "\\pgfplotsset{\n" + 
         "  every axis/.append style={\n" + 
         "    font=\\footnotesize,\n" + 
         "    width=.8\\columnwidth,\n" + 
         "    height=.6\\columnwidth,\n" + 
         "    %!line width=0.6pt,\n" + 
         "    legend style={at={(1.02,1)}, anchor={north west}},\n" + 
         "    cycle list name=defaultcolorlist,\n" + 
         "  },\n" + 
         "  every axis title/.style={\n" + 
         "    at={(0.5,1)},\n" + 
         "    below\n" + 
         "  },\n }\n" + 
         "%\n" + 
         "\\begin{document}\n\n");
        return sb.toString();
	}
	
	public static String pgfplotEndDocument () {
		return "\\end{document}\n";
	}
}
