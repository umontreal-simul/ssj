package tutorial;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.*;
import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;
import java.util.Collection;
import java.util.Arrays;

@RunWith(Parameterized.class)
public class CompareOutputs {

    final static Pattern ignorePat = Pattern.compile(".*(\\bCPU time|\\bEfficiency ratio).*");

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][] {
                 { "tutorial", AsianGBM.class,        null },
                 { "tutorial", AsianGBMQMC.class,     null },
                 { "tutorial", BankEv.class,       null },
                 { "tutorial", CallCenter.class,   new String[]{"tutorial/CallCenter.dat"} },
             //  { "tutorial", Collision.class,    null },
                 { "tutorial", InventoryCRN.class, null },
                 { "tutorial", Inventory.class,    null },
            //   { "tutorial", Nonuniform.class,   null },
                 { "tutorial", PreyPred.class,     null },
                 { "tutorial", QueueEv.class,      null },
                 { "tutorial", QueueLindley.class, null },
                 { "tutorial", QueueObs.class,     null },
                 // { "probdistmulti/norta", ExampleNortaInitDisc.class, null }
        });
    }

    private String prefix;
    private Class targetClass;
    private String[] args;

    // Constructor.
    public CompareOutputs(String prefix, Class targetClass, String[] args) {
        this.prefix = prefix;
        this.targetClass = targetClass;
        this.args = args;
    }

    @Test
    public void runTutorialClass() throws RunClass.RunClassException, IOException {
        String expected = RunClass.readFile(new File(prefix, targetClass.getSimpleName() + ".res"));
        String actual = RunClass.run(targetClass, args);
        RunClass.compareLineByLine(targetClass.getName(), expected, actual, ignorePat);
    }
}
