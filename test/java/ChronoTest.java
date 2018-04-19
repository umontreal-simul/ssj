import org.junit.Test;
import umontreal.ssj.util.GlobalCPUTimeChrono;

public class ChronoTest {
   @Test
   public void testChrono() {
      GlobalCPUTimeChrono timer = new GlobalCPUTimeChrono();
      System.out.println(timer.format());
      // no assert; just test execution
   }
}
