### Version 2.0

- Added branches regarding the "true" hits from sim:IDE

| Branch | Type | Description |
| ------ | ---- | ----------- |
| simide_size | int | size of the stored sim:IDE steps |
| simide_channelIdY | vector<int> | channel id |
| simide_trackId | vector<int> | track id |
| simide_tdc | vector<int> | tdc number |
| simide_x | vector<float> | x position of the step |
| simide_y | vector<float> | y position of the step |
| simide_z | vector<float> | z position of the step |
| simide_numElectrons | vector<float> | number of deposited electrons of the step |

- The channel id only run though the collection wires (only need to use one wire plane).
- A negative track id means this is a EM particle (for example, track id -1 means this is a EM daughter from track 1).
- There could be (and often are) many steps within one tdc.