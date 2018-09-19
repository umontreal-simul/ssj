package umontreal.ssj.stat;

//import   org.jfree.chart.ChartFactory;
//import   org.jfree.chart.ChartPanel;
//import   org.jfree.chart.axis.NumberAxis;
//import   org.jfree.chart.plot.XYPlot;
//import   org.jfree.chart.plot.PlotOrientation;
//import   org.jfree.data.statistics.HistogramBin;
//import org.omg.CORBA.portable.OutputStream;
//import   cern.colt.list.DoubleArrayList;
//import   java.util.ListIterator;
import java.util.Locale;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Formatter;

public class HistogramChartToLatex {

	/*
	 * return the min of the array of height
	 * 
	 */
	public double getMinY(double[] height) {
		double min = height[0];
		for (int i = 1; i < height.length; i++) {
			if (height[i] < min)
				min = height[i];
		}
		return min;
	}

	/*
	 * 
	 * Return the max of array
	 */
	public double getMaxY(double[] height) {
		double max = height[0];
		for (int i = 1; i < height.length; i++) {
			if (height[i] > max)
				max = height[i];
		}
		return max;
	}

	/*
	 * 
	 * return a array that contain the bound of the histogram with an TallyHistogram object
	 * 
	 * @param scaledH Scaledhistogram
	 */
	public double[] getHistogramBound(ScaledHistogram scaledH) {
		double a = scaledH.getA();
		double b = scaledH.getB();
		int nBin = scaledH.getNumBins();
		double h = (b - a) / nBin;
		double bound[] = new double[nBin + 1];
		bound[0] = a;
		for (int i = 1; i <= nBin; i++)
			bound[i] = bound[i - 1] + h;
		return bound;
	}

	/*
	 * 
	 * return tex file for an histogram by using a Scaledhistogram object and can add the polygonal
	 * histogram if the flag 'poly' is false.
	 * 
	 * @param scaledH a ScaledHistogram object
	 * 
	 * @param poly Plots a polygonal (piecewise linear) interpolation of the histogram if
	 * 'poly=true'
	 * 
	 * @param hist Plot the histogram if 'hist=true'
	 * 
	 * If 'poly=true' and 'hist=false', the method plot only the Polygonal interpolation of the
	 * histogram If 'poly=false' and 'hist=true', the method plot only the Histogram 'poly=true' and
	 * 'hist=true', the plot the histogram and the the polygonal We need to plot some things. If
	 * 'poly=false' and 'hist=false', the method plot the histogram only.
	 * 
	 * 
	 */

	public String toLatex(ScaledHistogram scaledH, boolean poly, boolean hist) {

		double height[] = scaledH.getHeights();
		double bound[] = getHistogramBound(scaledH);
		return toLatex(bound, height, poly, hist);

	}

	/*
	 * 
	 * 
	 * return tex file for an histogram by using the array of bound and the array of height of the
	 * histogram
	 * 
	 * @param bound The array that contain the bound of the histogram
	 * 
	 * @param height The array of rescaled counters: height[j] is the height of bin j.
	 * 
	 * @param poly Plots a polygonal (piecewise linear) interpolation of the histogram if
	 * 'poly=true'
	 * 
	 * @param hist Plot the histogram if 'hist=true'
	 * 
	 * 
	 * If 'poly=true' and 'hist=false', the method plot only the Polygonal interpolation of the
	 * histogram If 'poly=false' and 'hist=true', the method plot only the Histogram 'poly=true' and
	 * 'hist=true', the method plot the histogram and the the polygonal We need to plot some things.
	 * If 'poly=false' and 'hist=false', the method plot the histogram only.
	 * 
	 * 
	 */
	public String toLatex(double[] bound, double[] height, boolean poly, boolean hist) {
		if (poly == false && hist == false) {
			hist = true;
		}

		double h = (bound[1] - bound[0]) / 2;
		double yMin = getMinY(height);
		double yMax = getMaxY(height) + 0.2 * getMaxY(height);
		Formatter formatter = new Formatter(Locale.US);
		formatter.format("\\documentclass[border=3mm, %n");
		formatter.format("           tikz, %n");
		formatter.format("           preview %n");
		formatter.format("           ]{standalone} %n%n");
		formatter.format("\\usepackage{pgfplots}");
		formatter.format("%n%n");
		formatter.format("\\begin{document}%n%n");
		formatter.format("%%---------------------------------------------------------------%%%n");
		formatter.format("\\begin{tikzpicture} %n%n");
		formatter.format("\\begin{axis}[ %n");
		formatter.format("        ymin=%s, ymax=%s,%n", yMin, yMax);
		formatter.format("        minor y tick num = 3, %n");
		formatter.format("        %%area style, %n");
		formatter.format("        ] %n");
		if (hist) {
			formatter.format("\\addplot+[ybar interval,mark=no] plot coordinates { ");
			for (int i = 0; i < height.length; i++)
				formatter.format("\n (%s,%s) ", bound[i], height[i]);
			formatter.format("};%n%n");
		}
		if (poly) {
			formatter.format("\\addplot+[sharp plot] plot coordinates { ");
			formatter.format("(%s,%s) ", bound[0], height[0]);
			for (int i = 0; i < height.length; i++)
				formatter.format("\n (%s,%s) ", bound[i] + h, height[i]);
			formatter.format("(%s,%s) ", bound[bound.length - 1], height[height.length - 1]);
			formatter.format("};%n%n");
		}
		formatter.format("\\end{axis} %n%n");
		formatter.format("\\end{tikzpicture} %n%n");
		formatter.format("\\end{document}%n");

		String ch = formatter.toString();
		formatter.close();
		return ch;
	}

	/*
	 * 
	 * 
	 * write a string in a file
	 * 
	 * @param name The name of the generate tex file. You can find it in your current working
	 * Directory you must change this path
	 * 
	 * @param chaine is the String to write
	 */
	public void writeStringTofile(String name, String chaine) throws IOException {
		String currentDir = System.getProperty("user.dir");
		File file = new File(currentDir + "/" + name + ".tex");
		if (file.exists()) {
			file.delete(); // you might want to check if delete was successfull
		}
		file.createNewFile();
		final FileOutputStream stream = new FileOutputStream(file, true);
		final PrintWriter out = new PrintWriter(new OutputStreamWriter(stream));
		out.println(chaine);                                                                                                                                                                                                        // Ajouter
		out.close();
	}

}
