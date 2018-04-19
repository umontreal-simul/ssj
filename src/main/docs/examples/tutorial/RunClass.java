package tutorial;

import static org.junit.Assert.*;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Pattern;

public class RunClass {

    static class RunClassException extends Exception {
        public RunClassException(String msg) {
            super(msg);
        }
    }

    // This class provides threads that read the streams from subprocesses to
    // avoid blocking the main thread.
    private static class StreamSink extends Thread {
        private InputStream inputStream;
        private ByteArrayOutputStream baos = new ByteArrayOutputStream();
        private byte[] buffer = new byte[1024];

        private StreamSink(InputStream inputStream) {
            this.inputStream = inputStream;
        }

        @Override
        public void run() {
            int count = 0;
            try {
                while ((count = inputStream.read(buffer)) != -1) {
                    baos.write(buffer, 0, count);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }

        @Override
        public String toString() {
            try {
                return baos.toString("UTF-8");
            }
            catch (Exception e) {
                return e.toString();
            }
        }
    }

    public static String run(Class prog) throws RunClassException {
        return run(prog, null);
    }

    // We need to run programs in a distinct JVM because the initial states of
    // the generators need to be reset for every program.
    public static String run(Class prog, String[] args) throws RunClassException {

        if (args == null)
            args = new String[]{};

        String separator = System.getProperty("file.separator");
        String classpath = System.getProperty("java.class.path");
        String path = System.getProperty("java.home") + separator + "bin" + separator + "java";

        ArrayList<String> cmd = new ArrayList<String>(5 + args.length);
        cmd.add(path);
        cmd.add("-cp");
        cmd.add(classpath);
        cmd.add("-Duser.language=C");
        cmd.add(prog.getName());
        for (String arg : args)
            cmd.add(arg);

        ProcessBuilder processBuilder = new ProcessBuilder().command(cmd);

        try {

            Process process = processBuilder.start();
            StreamSink outSink = new StreamSink(process.getInputStream());
            StreamSink errSink = new StreamSink(process.getErrorStream());
            outSink.start();
            errSink.start();
            process.waitFor();
            outSink.join();
            errSink.join();

            String err = errSink.toString();
            if (err.length() > 0) {
                System.err.println("==================== " + prog.getName() + " ====================");
                System.err.print(err);
                System.err.println();
            }

            String out = outSink.toString();
            if (out.length() > 0) {
                System.out.println("==================== " + prog.getName() + " ====================");
                System.out.print(out);
                System.out.println();
            }
            return out;
        }
        catch (Exception e) {
            throw new RunClassException(e.toString());
        }
    }

    public static String readFile(String fileName) throws IOException {
        return readFile(new File(fileName));
    }

    public static String readFile(File file) throws IOException {
        return new Scanner(file).useDelimiter("\\A").next();
    }

    public static List<String> splitLines(String s) {
        return new ArrayList<String>(Arrays.asList(s.split("\r?\n\r?")));
    }

    public static void compareLineByLine(
            String label,
            String expected,
            String actual,
            Pattern ignore) {

        Iterator<String> it1 = splitLines(expected).iterator();
        Iterator<String> it2 = splitLines(actual).iterator();
        int lineNum = 0;
        // check that every output line matches expected output
        while (it1.hasNext()) {
            String s1 = it1.next().trim();
            // if it2 has finished just skip empty lines in expected output
            if (s1.length() == 0 && !it2.hasNext())
                continue;
            String s2 = it2.next().trim();
            lineNum += 1;
            // ignore selected pattern
            if (ignore.matcher(s1).matches())
                continue;
            assertEquals(label + ":" + lineNum, s1, s2);
        }
        // check that there are no non-empty trailing lines
        int trailingLines = 0;
        while (it2.hasNext()) {
            if (it2.next().trim().length() > 0)
                trailingLines += 1;
        }
        assertEquals(label + ":" + lineNum + " trailing lines", 0, trailingLines);
    }
}


