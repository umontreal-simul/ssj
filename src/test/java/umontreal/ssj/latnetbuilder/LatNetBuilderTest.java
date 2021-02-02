package umontreal.ssj.latnetbuilder;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import umontreal.ssj.hups.Rank1Lattice;
import umontreal.ssj.hups.DigitalNetBase2;

public class LatNetBuilderTest {

	@Test
	public void testOrdinaryLatticeSearch(){
		System.out.println("====================");
		OrdinaryLatticeSearch search = new OrdinaryLatticeSearch();
		search.setDimension("5");
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
		search.setDimension("5");
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
		search.setDimension("5");
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
		search.setDimension("5");
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
		search.setDimension("5");
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
		search.setPathToLatNetBuilder("/Users/pierre/Documents/stage-2018/code/latnetsoft/bin/latnetbuilder");
		search.setInterlacing("4");
		search.setDimension("5");
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

	@Test
	public void testJSONLoader(){
		String testJsonPolynomial = "{"+
			"\"pathToLatNetBuilder\":\"/Users/pierre/Documents/stage-2018/code/latnetsoft/bin/latnetbuilder\","+
			"\"pathToOutputFolder\":\"latnetbuilder_results\","+
			"\"dimension\":\"15\","+
			"\"sizeParameter\":\"2^3\","+
			"\"figure\":\"CU:IC2\","+
			"\"weights\":[\"product:1.0\"],"+
			"\"normType\":\"1\","+
			"\"explorationMethod\":\"random:100\","+
			"\"interlacing\":\"2\","+
			"\"pointSetType\":\"net\","+
			"\"construction\":\"polynomial\","+
			"\"filters\":[]"+
		"}";

		Search search = Search.fromJSON(testJsonPolynomial);
		DigitalNetBase2 pointSet = (DigitalNetBase2) search.search();
		assertEquals(30, pointSet.getDimension());

		String testJsonOrdinary = "{"+
			"\"pathToLatNetBuilder\":\"/Users/pierre/Documents/stage-2018/code/latnetsoft/bin/latnetbuilder\","+
			"\"pathToOutputFolder\":\"latnetbuilder_results\","+
			"\"dimension\":\"15\","+
			"\"sizeParameter\":\"2^3\","+
			"\"figure\":\"CU:P2\","+
			"\"weights\":[\"product:1.0\"],"+
			"\"normType\":\"2\","+
			"\"explorationMethod\":\"random:100\","+
			"\"pointSetType\":\"lattice\","+
			"\"construction\":\"ordinary\","+
			"\"filters\":[]"+
		"}";

		search = Search.fromJSON(testJsonOrdinary);
		Rank1Lattice lattice = (Rank1Lattice) search.search();
		assertEquals(15, lattice.getDimension());
	}
}
