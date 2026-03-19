/**
 * @package umontreal.ssj.rng
 *
 * Streams of independent uniform random numbers.
 *
 * This package offers basic facilities for generating multiple streams and substreams
 * of uniform random numbers over the interval @f$U(0,1)@f$ or over a range of integer values. 
 * The design is based on the interface @ref RandomStream and the package offers various 
 * implementations of this interface, with recurrence-based (sequential) random number generators. 
 * The interface specifies that each stream of random numbers is
 * partitioned into multiple substreams and that methods are available to
 * jump between the substreams, as discussed in @cite rLEC91a, @cite rLEC97d,
 * @cite rLEC02a, @cite rLEC17p. For examples of how to use these streams
 * properly, see  @cite rLEC15a and the `InventoryCRN.java` program in the tutorial examples. 
 *
 * Each implementation uses a specific backbone uniform random number
 * generator (RNG), whose period length is typically partitioned into very
 * long non-overlapping segments to provide the streams and substreams. A
 * stream can generate uniform variates (real numbers) over the interval
 * (0,1), uniform integers over a given range of values @f$\{i,...,j\}@f$, and
 * arrays of these.
 *
 * The generators provided are all recommendable.
 * They have been selected to be reasonably fast, to have a reasonably 
 * long period, good multivariate uniformity (based on theory), good statistical behavior,
 * and the capacity to implement the interface effectively.  
 * The @ref LFSR113 generator produces 
 * sequences of bits that obey a linear recurrence, so they eventually fail
 * statistical tests that measure the linear complexity of these binary
 * sequences. But this can affect only very special types of applications.
 *
 * The following tables give the approximate period length (period), 
 * the CPU time (in seconds) to generate @f$10^9@f$
 * @f$U(0,1)@f$ random numbers (gen. time), and the CPU time to jump ahead
 * @f$10^9@f$ times to the next substream (jump time). These timings
 * were on a (old) 2100 MHz 32-bit AMD Athlon XP 2800+ computer running Linux, with
 * the JDK 1.4.2.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">RNG</td>
 *   <td class="l bl">period</td>
 *   <td class="c">gen. time</td>
 *   <td class="r br">jump time</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">LFSR113<span style="height: 11.0pt;"></span></td>
 *   <td class="l bl">@f$2^{113}@f$</td>
 *   <td class="c">&ensp;51</td>
 *   <td class="r br">0.08</td>
 * </tr><tr>
 *   <td class="l bl br">WELL512</td>
 *   <td class="l bl">@f$2^{512}@f$</td>
 *   <td class="c">&ensp;55</td>
 *   <td class="r br">372</td>
 * </tr><tr>
 *   <td class="l bl br">WELL1024</td>
 *   <td class="l bl">@f$2^{1024}@f$</td>
 *   <td class="c">&ensp;55</td>
 *   <td class="r br">1450</td>
 * </tr><tr>
 *   <td class="l bl br">MT19937</td>
 *   <td class="l bl">@f$2^{19937}@f$</td>
 *   <td class="c">&ensp;56</td>
 *   <td class="r br">60</td>
 * </tr><tr>
 *   <td class="l bl br">WELL607</td>
 *   <td class="l bl">@f$2^{607}@f$</td>
 *   <td class="c">&ensp;61</td>
 *   <td class="r br">523</td>
 * </tr><tr>
 *   <td class="l bl br">MRG31k3p</td>
 *   <td class="l bl">@f$2^{185}@f$</td>
 *   <td class="c">&ensp;66</td>
 *   <td class="r br">1.8</td>
 * </tr><tr>
 *   <td class="l bl br">MRG32k3a</td>
 *   <td class="l bl">@f$2^{191}@f$</td>
 *   <td class="c">109</td>
 *   <td class="r br">2.3</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * The following timings were made on a 2400 MHz 64-bit AMD Athlon 64 Processor
 * 4000+ computer running Linux, with the JDK 1.5.0.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">RNG</td>
 *   <td class="l bl">period</td>
 *   <td class="c">gen. time</td>
 *   <td class="r br">jump time</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">LFSR113<span style="height: 11.0pt;"></span></td>
 *   <td class="l bl">@f$2^{113}@f$</td>
 *   <td class="c">&ensp;31</td>
 *   <td class="r br">0.08</td>
 * </tr><tr>
 *   <td class="l bl br">WELL512</td>
 *   <td class="l bl">@f$2^{512}@f$</td>
 *   <td class="c">&ensp;33</td>
 *   <td class="r br">234</td>
 * </tr><tr>
 *   <td class="l bl br">LFSR258</td>
 *   <td class="l bl">@f$2^{258}@f$</td>
 *   <td class="c">&ensp;35</td>
 *   <td class="r br">0.18</td>
 * </tr><tr>
 *   <td class="l bl br">MRG31k3p</td>
 *   <td class="l bl">@f$2^{185}@f$</td>
 *   <td class="c">&ensp;51</td>
 *   <td class="r br">0.89</td>
 * </tr><tr>
 *   <td class="l bl br">MRG32k3a</td>
 *   <td class="l bl">@f$2^{191}@f$</td>
 *   <td class="c">&ensp;70</td>
 *   <td class="r br">1.1</td>
 * </tr>
 * </table>
 * </center>
 *
 * Other tools offered in this package are as follows:
 * @ref RandomStreamManager give tools to manage and synchronize
 * several streams simultaneously,
 * @ref BasicRandomStreamFactory permits one to create random stream
 * factories for a given type of stream, and 
 * @ref AntitheticStream, @ref BakerTransformedStream, and @ref TruncatedRandomStream
 * permit one to apply automatic transformations to the output of a given stream.
 *
 * For further details about uniform RNGs, we refer the reader to
 * @cite rKNU98a, @cite rLEC06h, @cite rLEC12a, @cite rLEC15a, 
 * @cite rLEC17h, @cite rLEC21a.
 *
 */
 
 package umontreal.ssj.rng;
 