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
| simide_tdc | `vector<int>` | tdc number |
| simide_x | `vector<float>` | x position of the step |
| simide_y | `vector<float>` | y position of the step |
| simide_z | `vector<float>` | z position of the step |
| simide_numElectrons | `vector<float>` | number of deposited electrons of the step |

- The channel id only run though the collection wires (only need to use one wire plane).
- A negative track id means this is a EM particle (for example, track id -1 means this is a EM daughter from track 1).
- There could be (and often are) many steps within one tdc.