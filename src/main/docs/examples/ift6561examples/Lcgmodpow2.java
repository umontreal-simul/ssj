package ift6561examples;

import umontreal.ssj.util.PrintfFormat;

public class Lcgmodpow2 {

	public void printValues(long m, long a, long c, long x0, int n) {
		long x = x0;
		for (int i=0; i < n; i++) {
			x = (a * x + c) % m;
			System.out.print("x_" + i + " ~=~ " + x + " &=& ");
			System.out.println(PrintfFormat.formatBase (2, x) + "_2 ");
			}
		}

	public static void main(String[] args) {
		Lcgmodpow2 rng = new Lcgmodpow2();
		rng.printValues(4096 * 4096, 1140671485, 12820163, 12345, 10);
	}
}
