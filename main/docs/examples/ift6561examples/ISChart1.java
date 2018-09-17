package ift6561examples;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import umontreal.ssj.probdist.*;
import umontreal.ssj.util.Num;
import umontreal.ssj.charts.*;

// Produces plots for is-expon examples in simulation book.

public class ISChart1
{

   // First exponential example.
   private static double[][] isExpon1 (double lambda, double T, 
                    double xmin, double xmax, int numPoints)
   {
      double x, z, v;
      double[][] points = new double[2][numPoints];
      double p = Math.exp (-lambda * T);
      double delta = (xmax-xmin) / (numPoints-1);
      double min = lambda + 1.0/T - Math.sqrt(lambda * lambda + 1/(T*T));
      System.out.println("Minimal lambda0 = " + min);
      for (int i = 0; i < numPoints; i++) 
      {
         x =  xmin + i * delta;
         z = (2.0 * lambda - x);
         v = (Math.exp (-z * T)) * lambda * lambda / (x * z);
         points[0][i] = x;
         points[1][i] = (v - p*p) / (p * (1.0-p));
      }
      return points;
   }

   // Third exponential example.
   private static double[][] isExpon3(double lambda, double T,
				  double xmin, double xmax, int numPoints)
   {
	   double x, z, v;
	   double[][] points = new double[2][numPoints];
	   double p = 1.0 - Math.exp(-lambda * T);
	   double delta = (xmax - xmin) / (numPoints - 1);
	   for (int i = 0; i < numPoints; i++)
	   {
		   x = xmin + i * delta;
		   z = (2.0 * lambda - x);
		   if (Math.abs(z) < 0.00001)
			   v = lambda * lambda * T / x;
		   else
			   v = (1.0 - Math.exp(-z * T)) * lambda * lambda / (x * z);
		   points[0][i] = x;
		   points[1][i] = (v - p * p) / (p * (1.0 - p));
	   }
	   // z = (2.0 * lambda - 3.5);
	   // v = (1.0 - Math.exp (-z * T));
	   // v *= lambda * lambda / (3.5 * z);
	   // v = (v - p*p) / (p * (1-p));
	   // System.out.println ("v = " + v);
	   return points;
   }

   // Binomial example in rare07 paper.
   private static double[][] isBinomial (int n, int na, double p,
                 double xmin, double xmax, int numPoints)
   {
      double q, z, v, w, sum;
      BinomialDist bin = new BinomialDist (n, p);
      double gamma = bin.barF (na);         // Exact probability.
      double[][] points = new double[2][numPoints];
      double delta = (xmax-xmin) / (numPoints-1);
      for (int j = 0; j < numPoints; j++) 
      {
         q =  xmin + j * delta;   // q in the paper.
         sum = 0.0;
         v = p * p / q;
         w = (1.0-p) * (1.0-p) / (1.0-q);
         z = Math.exp (na * Math.log (v) + (n-na) * Math.log (w));
         for (int i = na; i <= n; i++)
         {
            sum += Num.combination (n, i) * z;
            z *= (v / w);
         }
         points[0][j] = q;
         points[1][j] = (sum - gamma * gamma) / (gamma * (1-gamma));
         // points[0][j] = -Math.log (q);
         // points[1][j] = -Math.log (sum - gamma * gamma);
      }
      return points;
   }

   // Constant function.
   private static double[][] constant (double xmin, double xmax, double y)
   {
      double[][] points = new double[2][2];
      points[0][0] = xmin;
      points[0][1] = xmax;
      points[1][0] = y;
      points[1][1] = y;
      return points;
   }


   public static void main (String[] args) throws IOException
   {

	  double[][] points = isExpon3(1.0, 0.15, 0.25, 45.0, 181);
	  double[][] flat = constant(0.0, 45.0, 1.0);
	  XYLineChart chart = new XYLineChart(null, "lambda0", "ratio", points, flat);
	  Writer file = new FileWriter("is-expon3-raw.tex");
	  file.write(chart.toLatex(8, 5));
	  file.close();

      points = isExpon1 (1.0, 4.0, 0.005, 1.20, 241);
      flat = constant (0.0, 1.20, 1.0);
      chart = new XYLineChart (null, "lambda0", "ratio", points, flat);
	  file = new FileWriter("is-expon1-raw.tex");
	  file.write(chart.toLatex(8, 5));
	  file.close();
//      System.out.println (chart.toLatex (8, 5));

	  /*     
			points = isBinomial (20, 15, 0.25, 0.25, 0.98, 150);
			flat = constant (0.20, 1.0, 1.0);
			chart = new XYLineChart (null, "q", "ratio", points, flat);
			System.out.println (chart.toTikz (8, 5));
	   */ 
   }

}

