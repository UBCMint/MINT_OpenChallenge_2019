package ca.ubc.best.mint.museandroidapp.analysis;

import android.util.Log;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.io.*;

import ca.ubc.best.mint.museandroidapp.ParcelableResults;
import ca.ubc.best.mint.museandroidapp.vm.FlankerLiveRecorder;
import eeg.useit.today.eegtoolkit.Constants;
import eeg.useit.today.eegtoolkit.model.TimeSeriesSnapshot;

import static ca.ubc.best.mint.museandroidapp.Util.msToSamples;

/** Wraps results with post-processing logic. */
public class ResultsPostProcessing {
  // NOTE: These should be multiples of 100, as each sample covers 100ms.
  public static final int INITIAL_MS = -800; // Start recording 800ms *before* stimulus for baseline
  public static final int ALPHA_END_MS = 1200; // Alpha suppression ends around 1200ms post stimulus
  public static final int BETA_END_MS = 1300; // Beta suppression ends around 1300ms post stimulus

  private static final int BASELINE_SAMPLES = msToSamples(500); // -800ms to -300ms.
  private static final int TRIM_SAMPLES = msToSamples(Math.abs(INITIAL_MS));

  private static final int ALPHA_SUPPRESSION_START_SAMPLES = msToSamples(200);
  private static final int ALPHA_SUPPRESSION_END_SAMPLES = msToSamples(ALPHA_END_MS);

  private static final int BETA_SUPPRESSION_START_SAMPLES = msToSamples(800);
  private static final int BETA_SUPPRESSION_END_SAMPLES = msToSamples(BETA_END_MS);

  private static final int threshold = 0;

  private ResultsPostProcessing() { /* Don't create. */ }

  /**
   * Post-processing of results:
   *   1) For each epoch, subtract baseline and trim to simulus time.
   *   2) Calculate average suppression over the desired time periods.
   *   3) Calculate tap statistics (reaction time and accuracy).
   */
  public static ParcelableResults process(FlankerLiveRecorder recorder) {
    List<Map<String, TimeSeriesSnapshot<Double>>> alpha = processAllAlpha(recorder.getAlphaEpochs());
    List<Map<String, TimeSeriesSnapshot<Double>>> beta = processAllBeta(recorder.getBetaEpochs());
    Date timeOfExperiment = new Date();

    return new ParcelableResults(
        alpha, beta,
        calcAlphaSuppression(alpha),
        calcBetaSuppression(beta),
        timeOfExperiment,
        recorder.getAverageTapReactionTimeMs(),
        recorder.getTapAccuracyProportion()
    );
  }

  // Process all alpha epochs by processing each individial snapshot separately.
  public static List<Map<String, TimeSeriesSnapshot<Double>>> processAllAlpha(
      List<Map<String, TimeSeriesSnapshot<Double>>> alphas) {
    List<Map<String, TimeSeriesSnapshot<Double>>> results = new ArrayList<>();
    for (Map<String, TimeSeriesSnapshot<Double>> alpha : alphas) {
      Map<String, TimeSeriesSnapshot<Double>> snapshots = new HashMap<>();
      for (String snapshotKey : alpha.keySet()) {
        snapshots.put(snapshotKey, baselineCorrectAndTrim(alpha.get(snapshotKey)));
      }
      results.add(snapshots);
    }
    return results;
  }

  // Process all beta epochs by processing each individial snapshot separately.
  public static List<Map<String, TimeSeriesSnapshot<Double>>> processAllBeta(
          List<Map<String, TimeSeriesSnapshot<Double>>> betas) {
    List<Map<String, TimeSeriesSnapshot<Double>>> results = new ArrayList<>();
    for (Map<String, TimeSeriesSnapshot<Double>> beta : betas) {
      Map<String, TimeSeriesSnapshot<Double>> snapshots = new HashMap<>();
      for (String snapshotKey : beta.keySet()) {
        snapshots.put(snapshotKey, baselineCorrectAndTrim(beta.get(snapshotKey)));
      }
      results.add(snapshots);
    }

    return results;
  }

  /** @return Alpha suppression, calculated as the average drop in the desired time frame. */
  private static double calcAlphaSuppression(List<Map<String, TimeSeriesSnapshot<Double>>> alphas) {
    Double[] alphaSurpression = parseAlpha(alphas);
    return alphaSurpression[0];
  }

  // Saves the Alpha data and the corresponding time stamp into a text file
  private static Double[] parseAlpha(List<Map<String, TimeSeriesSnapshot<Double>>> alphas) {
    ArrayList<Double> alphaData = new ArrayList<>();
    //long[] timeStamps;
    for(int i = 0; i < alphas.size(); i++) {
      for (Map.Entry<String, TimeSeriesSnapshot<Double>> test : alphas.get(i).entrySet()) {
        for(int j = 0; j < test.getValue().values.length; j++) {
          alphaData.add(test.getValue().values[j]);
        }
        //System.arraycopy(test.getValue().values, 0, outputData, 0, test.getValue().values.length);
        //timeStamps = new long[test.getValue().timestamps.length];
        //System.arraycopy(test.getValue().timestamps, 0, timeStamps, 0, test.getValue().timestamps.length);
      }
    }
    Double[] alphaArray = new Double[alphaData.size()];
    for(int i = 0; i < alphaArray.length; i++) {
      alphaArray[i] = alphaData.get(i);
    }
    return alphaThreshold(alphaArray);

  }

  private static Double[] alphaThreshold(Double[] alphaData) {
    Double[] beforeAlpha = new Double[alphaData.length];
    Double[] afterAlpha = new Double[alphaData.length];
    Double[] filteredDips = new Double[alphaData.length];
    Double[] derivativeAlpha = getDerivative(alphaData);
    Double[] averageVoltages = new Double[2];
    Map<Integer, Double> alphaDips = new HashMap<>();

    double sumBefore = 0;
    double sumAfter = 0;
    int counter = 0;

    for(int i = 0; i < derivativeAlpha.length; i++) {
      if(derivativeAlpha[i] < threshold) {
        alphaDips.put(i, alphaData[i]);
      }
    }

    for (Map.Entry<Integer, Double> test : alphaDips.entrySet()) {
      alphaData[test.getKey() - 1] = beforeAlpha[counter];
      alphaData[test.getKey() + 1] = afterAlpha[counter];
      filteredDips[counter] = test.getValue();
      counter++;
    }

    for(int j = 0; j < counter; j++) {
      sumBefore += Math.abs(filteredDips[j] - beforeAlpha[j]);
    }

    for(int k = 0; k < counter; k++) {
      sumAfter += Math.abs(filteredDips[k] - afterAlpha[k]);
    }

    averageVoltages[0] = sumBefore/counter;
    averageVoltages[1] = sumAfter/counter;

    return averageVoltages;
  }

  private static Double[] getDerivative(Double[] alphaData) {
    Double[] output = {0.0};
    return output;
  }

  /** @return Beta suppression, calculated as the average drop in the desired time frame. */
  private static double calcBetaSuppression(List<Map<String, TimeSeriesSnapshot<Double>>> betas) {
    double sum = 0.0;
    int count = 0;
    for (Map<String, TimeSeriesSnapshot<Double>> epochs : betas) {
      for (TimeSeriesSnapshot<Double> snapshot : epochs.values()) {
        sum += averageInRange(
            snapshot.values,
            BETA_SUPPRESSION_START_SAMPLES,
            BETA_SUPPRESSION_END_SAMPLES
        ) * -1.0; // Negate, as higher suppression = lower values.
        count++;
      }
    }
    return sum / count;
  }

  /**
   * Given a snapshot, first calculate a baseline rate using the average of the first samples.
   * Apply baseline correction by subtracting this from everything.
   * Finally, skip the initial values as they are before the stimuli actually happened, and only
   * used for the baseline calculation
   */
  private static TimeSeriesSnapshot<Double> baselineCorrectAndTrim(
      TimeSeriesSnapshot<Double> snapshot
  ) {
    Log.i("MINT", "Base and Trim, " + snapshot.length + " values," +
        " base = " + BASELINE_SAMPLES + " trim = " + TRIM_SAMPLES);
    double baseline = averageInRange(snapshot.values, 0, BASELINE_SAMPLES);

    int newCount = snapshot.length - TRIM_SAMPLES;
    long[] newTimes = new long[newCount];
    Double[] newValues = new Double[newCount];
    for (int i = 0; i < newCount; i++) {
      newTimes[i] = snapshot.timestamps[i + TRIM_SAMPLES];
      newValues[i] = snapshot.values[i + TRIM_SAMPLES] - baseline;
    }
    return new TimeSeriesSnapshot<>(newTimes, newValues);
  }

  /** @return The mean of all values[start..end). */
  private static double averageInRange(Double[] values, int start, int end) {
    double sum = 0;
    for(int i = 0; i < values.length; i++) {
      sum += values[i];
    }
    return sum / (values.length);
  }
}


