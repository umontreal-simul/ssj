package umontreal.ssj.latnetbuilder;

import org.junit.Test;

import umontreal.ssj.latnetbuilder.OrdinaryLatticeSearch;
import umontreal.ssj.latnetbuilder.DigitalNetSearch;
import umontreal.ssj.latnetbuilder.PolynomialLatticeSearch;

import umontreal.ssj.hups.Rank1Lattice;
import umontreal.ssj.hups.DigitalNetBase2;

public class LatNetBuilderTest {

	@Test
	public void testOrdinaryLatticeSearch(){
		System.out.println("====================");
		OrdinaryLatticeSearch search = new OrdinaryLatticeSearch();
		search.setDimension(5);
		search.setSizeParameter("2^10");
		search.setFigureOfMerit("CU:P2");
		search.setNormType("2");
		search.addWeight("product:1");
		search.setExplorationMethod("random-CBC:10");
		System.out.println(search.toString());
		System.out.println("Result:");
		System.out.println("Merit: " + search.merit());
		System.out.println("Time: " + search.time());
		Rank1Lattice lattice = search.search();
		System.out.println(lattice.toString());
	}

	@Test
	public void testDigitalNetSobol(){
		System.out.println("====================");
		DigitalNetSearch search = new DigitalNetSearch("sobol");
		search.setDimension(5);
		search.setSizeParameter("2^10");
		search.setFigureOfMerit("CU:P2");
		search.setNormType("2");
		search.addWeight("product:1");
		search.setExplorationMethod("random-CBC:10");
		System.out.println(search.toString());
		DigitalNetBase2 net = search.search();
		System.out.println("Result:");
		System.out.println("Merit: " + search.merit());
		System.out.println("Time: " + search.time());
		System.out.println(net.toString());
	    net.printGeneratorMatrices(5);
	}

	@Test
	public void testDigitalNetExplicit(){
		System.out.println("====================");
		DigitalNetSearch search = new DigitalNetSearch("explicit");
		search.setDimension(5);
		search.setSizeParameter("2^10");
		search.setFigureOfMerit("CU:P2");
		search.setNormType("2");
		search.addWeight("product:1");
		search.setExplorationMethod("random-CBC:10");
		System.out.println(search.toString());
		DigitalNetBase2 net = search.search();
		System.out.println("Result:");
		System.out.println("Merit: " + search.merit());
		System.out.println("Time: " + search.time());
		System.out.println(net.toString());
	    net.printGeneratorMatrices(5);
	}

	@Test
	public void testDigitalNetPolynomial(){
		System.out.println("====================");
		DigitalNetSearch search = new PolynomialLatticeSearch("net");
		search.setDimension(5);
		search.setSizeParameter("2^10");
		search.setFigureOfMerit("CU:P2");
		search.setNormType("2");
		search.addWeight("product:1");
		search.setExplorationMethod("random-CBC:10");
		System.out.println(search.toString());
		DigitalNetBase2 net = search.search();
		System.out.println("Result:");
		System.out.println("Merit: " + search.merit());
		System.out.println("Time: " + search.time());
		System.out.println(net.toString());
	    net.printGeneratorMatrices(5);
	}

	@Test
	public void testDigitalNetPolynomialLattice(){
		System.out.println("====================");
		DigitalNetSearch search = new PolynomialLatticeSearch("lattice");
		search.setDimension(5);
		search.setSizeParameter("2^10");
		search.setFigureOfMerit("CU:P2");
		search.setNormType("2");
		search.addWeight("product:1");
		search.setExplorationMethod("fast-CBC");
		System.out.println(search.toString());
		DigitalNetBase2 net = search.search();
		System.out.println("Result:");
		System.out.println("Merit: " + search.merit());
		System.out.println("Time: " + search.time());
		System.out.println(net.toString());
	    net.printGeneratorMatrices(5);
	}

	@Test
	public void testInterlacedDigitalNet(){
		System.out.println("====================");
		DigitalNetSearch search = new DigitalNetSearch("explicit");
		search.setInterlacing(4);
		search.setDimension(5);
		search.setSizeParameter("2^10");
		search.setFigureOfMerit("CU:IB");
		search.setNormType("1");
		search.addWeight("product:1");
		search.setExplorationMethod("random-CBC:10");
		System.out.println(search.toString());
		DigitalNetBase2 net = search.search();
		System.out.println("Result:");
		System.out.println("Merit: " + search.merit());
		System.out.println("Time: " + search.time());
		System.out.println(net.toString());
	    net.printGeneratorMatrices(5);
	}
}
