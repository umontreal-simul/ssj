package latnetbuilder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Optional;
import java.util.Set;

import com.google.gson.Gson;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.latnetbuilder.Search;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class Experiment {
    private String experimentName;
    private transient Search search;
    private String serializedSearch;
    private transient PointSet pointSet;
    private String serializedPointSet;
    private SobolTestFunc integrand;
    private String randomization;
    private double variance;
    private double varianceMC;
    private double latnetbuilderComputeTime;
    private double ssjComputeTime;
    public transient Gson gson;
    private int hash;

    private static ArrayList<String> listRandomizations = new ArrayList<String>(
        Arrays.asList("LMS - Shift", "NUS - Interlace", "LMS - Interlace - Shift"));

    public Experiment(String name, Search search, SobolTestFunc integrand, String randomization) {
        this.gson = new Gson();
        this.experimentName = name;
        setSearch(search);
        this.integrand = integrand;
        setRandomization(randomization);
    }

    public Experiment(String name, PointSet pointSet, SobolTestFunc integrand, String randomization) {
        this.gson = new Gson();
        this.experimentName = name;
        setPointSet(pointSet);
        this.integrand = integrand;
        setRandomization(randomization);
    }

    public void setRandomization(String randomization) {
        this.randomization = randomization;
        assert listRandomizations.contains(randomization);
    }

    public void setVarianceMC(double varianceMC) {
        this.varianceMC = varianceMC;
    }

    public void setVariance(double variance) {
        this.variance = variance;
    }

    public void setPointSet(PointSet pointSet) {
        this.pointSet = pointSet;
        this.serializedPointSet = gson.toJson(pointSet);
    }

    public void setSearch(Search search) {
        this.search = search;
        this.serializedSearch = gson.toJson(search);
    }

    public int getHash() {
        updateHash();
        return this.hash;
    }

    private void simulateLMSShift(int m, RandomStream noise, Tally statQMC) {
        Tally statValue = new Tally("stat on value simple test function");
        DigitalNetBase2 net = (DigitalNetBase2) this.pointSet; // This randomization only works for nets.
        PointSetIterator stream = net.iterator();
        for (int j = 0; j < m; j++) {
            net.leftMatrixScramble(noise);
            net.addRandomShift(0, net.getDimension(), noise);
            stream.resetStartStream();
            this.integrand.simulateRuns(net.getNumPoints(), stream, statValue);
            statQMC.add(statValue.average());
        }
    }

    private void simulateNUSInterlace(int m, RandomStream noise, Tally statQMC) {
        Tally statValue = new Tally("stat on value simple test function");
        DigitalNetBase2 net = (DigitalNetBase2) this.pointSet; // This randomization only works for nets.
        for (int j = 0; j < m; j++) {
            int[][] randomizedPoints = new int[net.getNumPoints()][net.getDimension()];
            double[][] interlacedPoints = new double[net.getNumPoints()][net.getDimension() / net.getInterlacing()];
            net.nestedUniformScramble(noise, randomizedPoints, 0);
            net.outputInterlace(randomizedPoints, interlacedPoints);
            this.integrand.simulateRuns(net.getNumPoints(), interlacedPoints, statValue);
            statQMC.add(statValue.average());
        }
    }

    private void simulateLMSInterlaceShift(int m, RandomStream noise, Tally statQMC) {
        Tally statValue = new Tally("stat on value simple test function");
        DigitalNetBase2 net = (DigitalNetBase2) this.pointSet; // This randomization only works for nets.
        for (int j = 0; j < m; j++) {
            net.leftMatrixScramble(noise);
            DigitalNetBase2 interlacedNet = net.matrixInterlace();
            interlacedNet.addRandomShift(0, interlacedNet.getDimension(), noise);
            PointSetIterator stream = interlacedNet.iterator();
            this.integrand.simulateRuns(interlacedNet.getNumPoints(), stream, statValue);
            statQMC.add(statValue.average());
        }
    }

    public void runExperiment(Path experimentDirPath, int experimentNumber) throws IOException {

        if (this.pointSet == null) {
            System.out.println("Computing pointSet " + experimentNumber);
            long startTime = System.nanoTime();
            this.pointSet = this.search.search();
            long endTime = System.nanoTime();
            this.latnetbuilderComputeTime = (endTime - startTime) / 1000000;
            this.serializedPointSet = gson.toJson(pointSet);
        }

        if (this.variance == 0){
            System.out.println("Computing variance " + experimentNumber);
            int m = 100; // Number of QMC randomizations.
            Tally statValue = new Tally("RQMC method");
            long startTime = System.nanoTime();
            if (this.randomization.equals("LMS - Shift")){
                simulateLMSShift(m, new MRG32k3a(), statValue);
            }
            else if (this.randomization.equals("NUS - Interlace")){
                simulateNUSInterlace(m, new MRG32k3a(), statValue);
            }
            else if (this.randomization.equals("LMS - Interlace - Shift")){
                simulateLMSInterlaceShift(m, new MRG32k3a(), statValue);
            }
            else {
                throw new RuntimeException("randomization not understood.");
            }
            long endTime = System.nanoTime();
            this.ssjComputeTime = (endTime - startTime) / 1000000;
            this.variance = this.pointSet.getNumPoints() * statValue.variance();
        }

        saveToFile(experimentDirPath, experimentNumber);
    }

    private void saveToFile(Path experimentDirPath, int experimentNumber) throws IOException {
        String json = gson.toJson(this);
        String fileName = experimentDirPath.resolve(experimentNumber + ".json").toString();
        try (PrintStream out = new PrintStream(new FileOutputStream(fileName))) {
            out.print(json);
        }
    }

    public void toJson(Path experimentDirPath, Set<Integer> experimentHashes) throws IOException {
        if (!experimentHashes.contains(getHash())) {
            saveToFile(experimentDirPath, experimentHashes.size());
            experimentHashes.add(getHash());
        }
    }

    public static Experiment[] fromJSONs(Path experimentDirPath) throws IOException {
        Gson gson = new Gson();
        File[] listFiles = experimentDirPath.toFile().listFiles();
        int length = listFiles.length;
        ArrayList<Experiment> experimentList = new ArrayList<Experiment>(length);
        for (int i = 0; i < length; i++) {
            String serializedJson = readFromInputStream(new FileInputStream(listFiles[i]));
            if (getExtension(listFiles[i].getName()).get().equals("json")) {
                Experiment experiment = gson.fromJson(serializedJson, Experiment.class);
                experiment.gson = gson;
                if (experiment.serializedSearch != null) {
                    experiment.setSearch(Search.fromJSON(experiment.serializedSearch));
                }
                if (experiment.serializedPointSet != null) {
                    PointSet pointSet;
                    if (experiment.serializedPointSet.contains("minit_from_file")) {
                        pointSet = gson.fromJson(experiment.serializedPointSet, SobolSequence.class);
                    } else {
                        pointSet = gson.fromJson(experiment.serializedPointSet, DigitalNetBase2.class);
                    }
                    experiment.setPointSet(pointSet);
                }
                experimentList.add(experiment);
            }
        }
        Experiment[] result = new Experiment[experimentList.size()];
        return experimentList.toArray(result);
    }

    private static Optional<String> getExtension(String filename) {
        return Optional.ofNullable(filename).filter(f -> f.contains("."))
                .map(f -> f.substring(filename.lastIndexOf(".") + 1));
    }

    private static String readFromInputStream(InputStream inputStream) throws IOException {
        StringBuilder resultStringBuilder = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(inputStream))) {
            String line;
            while ((line = br.readLine()) != null) {
                resultStringBuilder.append(line).append("\n");
            }
        }
        return resultStringBuilder.toString();
    }

    private void updateHash() {
        this.hash = gson.toJson(integrand).hashCode() + randomization.hashCode();
        if (this.serializedSearch != null) {
            this.hash += serializedSearch.hashCode();
        } else if (this.serializedPointSet != null) {
            this.hash += serializedPointSet.hashCode();
        } else {
            throw new RuntimeException("serializedSearch or serializedPointSet should be set.");
        }
    }

}
