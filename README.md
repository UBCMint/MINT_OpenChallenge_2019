# MINT OpenChallenge 2019 - Flank

![Logo](https://raw.githubusercontent.com/UBCMint/MuseAndroidApp/master/images/flank.png)

Flank is an Android App that can provide users with a collection of metrics regarding their attentional state, by connecting to a Muse EEG headset and analyzing signal during some behavioural testing.

_Disclaimer: The UBC MiNT App 'Flank' is a prototype tool for exploring attention data with a Muse headset.
It is not clinically approved, so the results may not be medically accurate and should not be interpreted as such._

### Background

Flank is named after the Erikson flanker task, a behavioural test where participants select a region according to a central stimulus, while it is flanked by distractor stimuli that may point in the same (congruent) or opposite (incongruent) direction, or even no direction at all (neutral).

In their paper "[Differential Oscillatory Electroencephalogram Between Attention-Deficit/Hyperactivity Disorder Subtypes and Typically Developing Adolescents](https://doi.org/10.1016/j.biopsych.2013.08.023)", Mazaher et al (2014) were able to pair the Flanker task with preceding cue stimuli, and find two frequency-based indicators of ADHD:
* _Reduced_ suppression of Alpha frequencies between 200ms and 1200ms after cue presented. During this period, a  participant is determining whether the cue actually corresponds to which side they should select.
* _Reduced_ suppression of Beta frequencies between 800ms and 1300ms after cue presented. During this period, a participant is sending motor commands to tap the correct side of the screen.

The Flank app implements the experiment protocol as closely as possible on an Android device using a Muse consumer EEG headset. Users can periodically test themselves using the phone-based Flanker tast, and get historic trends of their alpha- and beta-suppression, as well as overall reaction time and accuracy.
