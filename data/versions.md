### Version 1.0

| Branch | Type | Description |
| ------ | ---- | ----------- |
| calib_nChannel | `int` | number of hit channels |
| calib_channelId | `vector<int>` | channel id |
| calib_wf | `TClonesArray` of `TH1F*`s | waveform of the channel |

- Each waveform is stored as a TH1F, each has 9600 bins representing the 9600 tdcs.
- The unit of the waveform is "number of electrons" (or adc?).

### Version 2.0

- Added a TNamed Object `version` to describe the version of the data schema.
- Added branches regarding the "true" hits from sim:IDE.

| Branch | Type | Description |
| ------ | ---- | ----------- |
| simide_size | `int` | size of the stored sim:IDE steps |
| simide_channelIdY | `vector<int>` | channel id |
| simide_trackId | `vector<int>` | track id |
| simide_tdc | `vector<int>` | tdc number + TriggerOffsetTPC (3200) |
| simide_x | `vector<float>` | x position of the step |
| simide_y | `vector<float>` | y position of the step |
| simide_z | `vector<float>` | z position of the step |
| simide_numElectrons | `vector<float>` | number of deposited electrons of the step |

- The channel id runs though all wires on the three wire-planes, so the energy deposits are counted three times (on the other hand, there seems to be missing energy deposits if only run through one plane).
- A negative track id means this is a EM particle (for example, track id -1 means this is a EM daughter from track 1).
- There could be (and often are) many steps within one tdc.

### Version 3.0

- Added branches for the raw waveforms (before decovolution)
  - same formats as the calib waveforms
  - however size are much larger as all channels have the raw waveforms (with pedastal and noise)

| Branch | Type | Description |
| ------ | ---- | ----------- |
| raw_nChannel | `int` | number of raw channels |
| raw_channelId | `vector<int>` | channel id |
| raw_wf | `TClonesArray` of `TH1F*`s | waveform of the channel |

- Added branches regarding the MC truth
  - `MAX_TRACKS` is defined to be 30000 when extracting the MC tracks.

| Branch | Type | Description |
| ------ | ---- | ----------- |
| mc_Ntrack   | `int`              | number of tracks in MC |
| mc_id       | `int[MAX_TRACKS]`  | array of track id |
| mc_pdg      | `int[MAX_TRACKS]`  | array of track particle pdg |
| mc_process  | `int[MAX_TRACKS]`  | array of track generation process code |
| mc_mother   | `int[MAX_TRACKS]`  | array of mother id of each track |
| mc_daughters|  `vector<vector<int> >` | daughters id's of each track |
| mc_startXYZT    | `float[MAX_TRACKS][4]` | start position of each track [cm]|
| mc_endXYZT      | `float[MAX_TRACKS][4]` | end position of each track [cm]|
| mc_startMomentum| `float[MAX_TRACKS][4]` | start momentum of each track [GeV]|
| mc_endMomentum  | `float[MAX_TRACKS][4]` | end momentum of each track [GeV]|

- Added more branches regarding the neutrino interaction MC truth

| Branch | Type | Description |
| ------ | ---- | ----------- |
|mc_isnu | int |  is it neutrino interaction? |
|mc_nGeniePrimaries | int |  number of Genie primaries |
|mc_nu_pdg | int |  pdg code of neutrino |
|mc_nu_ccnc | int |  cc or nc |
|mc_nu_mode | int |  mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html |
|mc_nu_intType | int |  interaction type |
|mc_nu_target | int | interaction target  |
|mc_hitnuc | int |  hit nucleon? |
|mc_hitquark | int |  hit quark? |
| mc_nu_Q2 | double | Q^2 |
| mc_nu_W | double | W |
| mc_nu_X | double | X |
| mc_nu_Y | double | Y |
| mc_nu_Pt | double | Pt |
| mc_nu_Theta | double | angle relative to lepton |
| mc_nu_pos | float[4] | interaction position of nu [cm] |
| mc_nu_mom | float[4] | interaction momentum of nu [GeV] |
