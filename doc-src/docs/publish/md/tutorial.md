# Startup


Edit this to point at the path to where you installed the `neuropixel-utils` distribution or repo.



```matlab
% Load the Neuropixel Utils library
homeDir = string(java.lang.System.getProperty('user.home'));
npixUtilsDistRoot = fullfile(homeDir, 'local', 'repos', 'npxutils-apj-wip-01');
addpath(fullfile(npixUtilsDistRoot, 'Mcode'));
fprintf('Loaded Neuropixel Utils %s\n', npxutils.globals.version)
```


```text
Loaded Neuropixel Utils 0.5.0-SNAPSHOT
```



Set the default channel map file you want to use with your data. Some standard channel maps are included in the `map_files` directory in the `neuropixel-utils` distribution.



```matlab
setenv('NEUROPIXEL_MAP_FILE', fullfile(npxutils.globals.distroot, 'map_files', 'neuropixPhase3A_kilosortChanMap.mat'))
```



Optional: Configure Neuropixel Utils to generate paths beneath your data root directory, based on its naming conventions.




This is entirely optional; you can arrange your data however you'd like.



```matlab
myExampleDataDir = fullfile(homeDir, 'work', 'npxutils', 'example-data');
setenv('NEUROPIXEL_DATAROOT', myExampleDataDir)
```

# Load Raw Imec dataset


Here we construct the path to the original, raw dataset as it was recorded by SpikeGLX. You don't need to use this folder nesting structure.



```matlab
subject = 'Vinnie';
dateStr = '2018-08-17';
rawBinFile = 'Vinnie_20180817_All.imec.ap.bin';
imecFile = npxutils.generatePath(subject, 'raw_datasets', dateStr, rawBinFile)
```


```text
imecFile = '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.ap.bin'
```



We then use this path to construct an ImecDataset object which will facilitate access to the raw data and load the metadata from disk. You can pass the path of the Imec dataset in several ways:



   -  Path to raw `.imec.ap.bin` file or `imec.lf.bin` file. 
   -  Path to the containing directory, if there is only one set of `.bin` files within. 



To construct the ImecDataset object, pass the identifying path directly to the constructor.



```matlab
imec = npxutils.ImecDataset(imecFile)
```


```text
imec = 
  ImecDataset with properties:

          bytesPerSample: 2
                pathRoot: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17'
                fileStem: 'Vinnie_20180817_All'
          fileImecNumber: NaN
            creationTime: 7.3729e+05
               nChannels: 385
              fileTypeAP: 'ap'
              fileTypeLF: 'lf'
              nSamplesAP: 36177457
              nSamplesLF: 0
                    fsAP: 30000
                    fsLF: NaN
                  fsSync: NaN
        highPassFilterHz: 300
                  apGain: 500
                 apRange: [-0.6000 0.6000]
                  lfGain: 250
                 lfRange: []
                 adcBits: 10
             snsShankMap: "(1,2,480)(0:0:0:1)(0:1:0:1)(0:0:1:1)(0:1:1:1)(0:0:2:1)(0:1:2:1)(0:0:3:1)(0:1:3:1)(0:0:4:1)(0:1:4:1)(0:0:5:1)(0:1:5:1)(0:0:6:1)(0:1:6:1)(0:0:7:1)(0:1:7:1)(0:0:8:1)(0:1:8:1)(0:0:9:1)(0:1:9:1)(0:0:10:1)(0:1:10:1)(0:0:11:1)(0:1:11:1)(0:0:12:1)(0:1:12:1)(0:0:13:1)(0:1:13:1)(0:0:14:1)(0:1:14:1)(0:0:15:1)(0:1:15:1)(0:0:16:1)(0:1:16:1)(0:0:17:1)(0:1:17:1)(0:0:18:0)(0:1:18:1)(0:0:19:1)(0:1:19:1)(0:0:20:1)(0:1:20:1)(0:0:21:1)(0:1:21:1)(0:0:22:1)(0:1:22:1)(0:0:23:1)(0:1:23:1)(0:0:24:1)(0:1:24:1)(0:0:25:1)(0:1:25:1)(0:0:26:1)(0:1:26:1)(0:0:27:1)(0:1:27:1)(0:0:28:1)(0:1:28:1)(0:0:29:1)(0:1:29:1)(0:0:30:1)(0:1:30:1)(0:0:31:1)(0:1:31:1)(0:0:32:1)(0:1:32:1)(0:0:33:1)(0:1:33:1)(0:0:34:1)(0:1:34:1)(0:0:35:1)(0:1:35:1)(0:0:36:1)(0:1:36:1)(0:0:37:1)(0:1:37:0)(0:0:38:1)(0:1:38:1)(0:0:39:1)(0:1:39:1)(0:0:40:1)(0:1:40:1)(0:0:41:1)(0:1:41:1)(0:0:42:1)(0:1:42:1)(0:0:43:1)(0:1:43:1)(0:0:44:1)(0:1:44:1)(0:0:45:1)(0:1:45:1)(0:0:46:1)(0:1:46:1)(0:0:47:1)(0:1:47:1)(0:0:48:1)(0:1:48:1)(0:0:49:1)(0:1:49:1)(0:0:50:1)(0:1:50:1)(0:0:51:1)(0:1:51:1)(0:0:52:1)(0:1:52:1)(0:0:53:1)(0:1:53:1)(0:0:54:1)(0:1:54:1)(0:0:55:1)(0:1:55:1)(0:0:56:0)(0:1:56:1)(0:0:57:1)(0:1:57:1)(0:0:58:1)(0:1:58:1)(0:0:59:1)(0:1:59:1)(0:0:60:1)(0:1:60:1)(0:0:61:1)(0:1:61:1)(0:0:62:1)(0:1:62:1)(0:0:63:1)(0:1:63:1)(0:0:64:1)(0:1:64:1)(0:0:65:1)(0:1:65:1)(0:0:66:1)(0:1:66:1)(0:0:67:1)(0:1:67:1)(0:0:68:1)(0:1:68:1)(0:0:69:1)(0:1:69:1)(0:0:70:1)(0:1:70:1)(0:0:71:1)(0:1:71:1)(0:0:72:1)(0:1:72:1)(0:0:73:1)(0:1:73:1)(0:0:74:1)(0:1:74:1)(0:0:75:1)(0:1:75:0)(0:0:76:1)(0:1:76:1)(0:0:77:1)(0:1:77:1)(0:0:78:1)(0:1:78:1)(0:0:79:1)(0:1:79:1)(0:0:80:1)(0:1:80:1)(0:0:81:1)(0:1:81:1)(0:0:82:1)(0:1:82:1)(0:0:83:1)(0:1:83:1)(0:0:84:1)(0:1:84:1)(0:0:85:1)(0:1:85:1)(0:0:86:1)(0:1:86:1)(0:0:87:1)(0:1:87:1)(0:0:88:1)(0:1:88:1)(0:0:89:1)(0:1:89:1)(0:0:90:1)(0:1:90:1)(0:0:91:1)(0:1:91:1)(0:0:92:1)(0:1:92:1)(0:0:93:1)(0:1:93:1)(0:0:94:0)(0:1:94:1)(0:0:95:1)(0:1:95:1)(0:0:96:1)(0:1:96:1)(0:0:97:1)(0:1:97:1)(0:0:98:1)(0:1:98:1)(0:0:99:1)(0:1:99:1)(0:0:100:1)(0:1:100:1)(0:0:101:1)(0:1:101:1)(0:0:102:1)(0:1:102:1)(0:0:103:1)(0:1:103:1)(0:0:104:1)(0:1:104:1)(0:0:105:1)(0:1:105:1)(0:0:106:1)(0:1:106:1)(0:0:107:1)(0:1:107:1)(0:0:108:1)(0:1:108:1)(0:0:109:1)(0:1:109:1)(0:0:110:1)(0:1:110:1)(0:0:111:1)(0:1:111:1)(0:0:112:1)(0:1:112:1)(0:0:113:1)(0:1:113:0)(0:0:114:1)(0:1:114:1)(0:0:115:1)(0:1:115:1)(0:0:116:1)(0:1:116:1)(0:0:117:1)(0:1:117:1)(0:0:118:1)(0:1:118:1)(0:0:119:1)(0:1:119:1)(0:0:120:1)(0:1:120:1)(0:0:121:1)(0:1:121:1)(0:0:122:1)(0:1:122:1)(0:0:123:1)(0:1:123:1)(0:0:124:1)(0:1:124:1)(0:0:125:1)(0:1:125:1)(0:0:126:1)(0:1:126:1)(0:0:127:1)(0:1:127:1)(0:0:128:1)(0:1:128:1)(0:0:129:1)(0:1:129:1)(0:0:130:1)(0:1:130:1)(0:0:131:1)(0:1:131:1)(0:0:132:0)(0:1:132:1)(0:0:133:1)(0:1:133:1)(0:0:134:1)(0:1:134:1)(0:0:135:1)(0:1:135:1)(0:0:136:1)(0:1:136:1)(0:0:137:1)(0:1:137:1)(0:0:138:1)(0:1:138:1)(0:0:139:1)(0:1:139:1)(0:0:140:1)(0:1:140:1)(0:0:141:1)(0:1:141:1)(0:0:142:1)(0:1:142:1)(0:0:143:1)(0:1:143:1)(0:0:144:1)(0:1:144:1)(0:0:145:1)(0:1:145:1)(0:0:146:1)(0:1:146:1)(0:0:147:1)(0:1:147:1)(0:0:148:1)(0:1:148:1)(0:0:149:1)(0:1:149:1)(0:0:150:1)(0:1:150:1)(0:0:151:1)(0:1:151:0)(0:0:152:1)(0:1:152:1)(0:0:153:1)(0:1:153:1)(0:0:154:1)(0:1:154:1)(0:0:155:1)(0:1:155:1)(0:0:156:1)(0:1:156:1)(0:0:157:1)(0:1:157:1)(0:0:158:1)(0:1:158:1)(0:0:159:1)(0:1:159:1)(0:0:160:1)(0:1:160:1)(0:0:161:1)(0:1:161:1)(0:0:162:1)(0:1:162:1)(0:0:163:1)(0:1:163:1)(0:0:164:1)(0:1:164:1)(0:0:165:1)(0:1:165:1)(0:0:166:1)(0:1:166:1)(0:0:167:1)(0:1:167:1)(0:0:168:1)(0:1:168:1)(0:0:169:1)(0:1:169:1)(0:0:170:0)(0:1:170:1)(0:0:171:1)(0:1:171:1)(0:0:172:1)(0:1:172:1)(0:0:173:1)(0:1:173:1)(0:0:174:1)(0:1:174:1)(0:0:175:1)(0:1:175:1)(0:0:176:1)(0:1:176:1)(0:0:177:1)(0:1:177:1)(0:0:178:1)(0:1:178:1)(0:0:179:1)(0:1:179:1)(0:0:180:1)(0:1:180:1)(0:0:181:1)(0:1:181:1)(0:0:182:1)(0:1:182:1)(0:0:183:1)(0:1:183:1)(0:0:184:1)(0:1:184:1)(0:0:185:1)(0:1:185:1)(0:0:186:1)(0:1:186:1)(0:0:187:1)(0:1:187:1)(0:0:188:1)(0:1:188:1)(0:0:189:1)(0:1:189:0)(0:0:190:1)(0:1:190:1)(0:0:191:1)(0:1:191:1)"
              channelMap: [1x1 npxutils.ChannelMap]
             badChannels: [0x1 uint32]
            syncBitNames: [16x1 string]
     concatenationInfoAP: [1x1 npxutils.ConcatenationInfo]
     concatenationInfoLF: []
          sourceDatasets: []
                 syncRaw: []
                   hasAP: 1
                   hasLF: 0
       hasSourceDatasets: 0
             hasSourceAP: 0
             hasSourceLF: 0
           hasSourceSync: 0
          channelMapFile: "/home/janke/local/repos/npxutils-apj-wip-01/map_files/neuropixPhase3A_kilosortChanMap.mat"
          mappedChannels: [384x1 uint32]
       mappedChannelInds: [384x1 double]
         nChannelsMapped: 384
       connectedChannels: [374x1 uint32]
    connectedChannelInds: [374x1 double]
      nChannelsConnected: 374
            goodChannels: [374x1 uint32]
         goodChannelInds: [374x1 double]
           nGoodChannels: 374
              channelIds: [385x1 uint32]
            channelNames: [385x1 string]
      channelNamesPadded: [385x1 string]
               nSyncBits: 16
           syncBitsNamed: [11x1 double]
                  fileAP: 'Vinnie_20180817_All.imec.ap.bin'
                  pathAP: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.ap.bin'
              fileAPMeta: 'Vinnie_20180817_All.imec.ap.meta'
              pathAPMeta: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.ap.meta'
                  fileLF: 'Vinnie_20180817_All.imec.lf.bin'
                  pathLF: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.lf.bin'
              fileLFMeta: 'Vinnie_20180817_All.imec.lf.meta'
              pathLFMeta: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.lf.meta'
                fileSync: 'Vinnie_20180817_All.imec.ap.bin'
                pathSync: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.ap.bin'
          fileSyncCached: 'Vinnie_20180817_All.sync.mat'
          pathSyncCached: '/home/janke/work/npxutils/example-data/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.sync.mat'
         creationTimeStr: '17-Aug-2018 13:50:34'
             apScaleToUv: 2.3438
             lfScaleToUv: NaN
           syncChannelId: 385
        syncChannelIndex: 385
            syncInAPFile: 1
            syncInLFFile: 0
```

## Specifying the channel map manually


By default, the channel map will be loaded from whatever the environment variable NEUROPIXEL_MAP_FILE is set to, which we set above. If you want to manaully specify a channel map, pass it in to the ImecDataset constructor using the `channelMap` option:




imec = npxutils.ImecDataset(imecFile, 'channelMap', '/path/to/channelMap.mat');


## Explore the loaded metadata


Here we see that this option 3A probe has 384 total channels (saved in the `.imec.ap.bin` file), and that 374 of these channels are connected, namely those listed in `connectedChannels`.



```matlab
imec.channelMap
```


```text
ans = 
  ChannelMap with properties:

                 file: "/home/janke/local/repos/npxutils-apj-wip-01/map_files/neuropixPhase3A_kilosortChanMap.mat"
                 name: "neuropixPhase3A_kilosortChanMap"
     channelIdsMapped: [384x1 uint32]
            connected: [384x1 logical]
             shankInd: [384x1 double]
         nSpatialDims: 2
              xcoords: [384x1 double]
              ycoords: [384x1 double]
              zcoords: [384x1 double]
     syncChannelIndex: 385
        syncChannelId: 385
               coords: [384x2 double]
         syncInAPFile: 1
         syncInLFFile: 1
           channelIds: [385x1 uint32]
            nChannels: 385
      nChannelsMapped: 384
    connectedChannels: [374x1 uint32]
    referenceChannels: [10x1 uint32]
      invertChannelsY: 1
             yspacing: 20
             xspacing: 16
             zspacing: [0x1 double]
                 xlim: [-9 79]
                 ylim: [0 3860]
```



We can determine the duration of the recording using the number of samples and the frequency in the Imec data:



```matlab
durationMinutes = imec.nSamplesAP / imec.fsAP / 60;
fprintf('Duration of recording %s is %.1f minutes.\n', imec.fileStem, durationMinutes)
```


```text
Duration of recording Vinnie_20180817_All is 20.1 minutes.
```

# Accessing the raw AP data


You can access the raw AP data in several ways. (See the full list of methods with `methods(imec)` or `doc npxutils.ImecDataset`). The most straightforward is to memory-map the full binary data:




*This bit is broken because it references code that doesn't exist in this library. Probably left over from an older version of the library and never updated? Or maybe it has an external dependency?*



```matlab
% This is broken because NeuropixelExpt doesn't exist?

% NeuropixelExpt.DataLoad.setImecSyncBitNames(imec);

% What does "p" refer to here?

% if ~isempty(p.Results.badChannels)
%    imec.markBadChannels(p.Results.badChannels);
% end
```



*This stuff works:*



```matlab
rmsBadChannels = imec.markBadChannelsByRMS('rmsRange', [3 100]);
fprintf('Marked %d channels bad based on RMS\n', numel(rmsBadChannels));
```


```text
Marked 0 channels bad based on RMS
```


```matlab
% Save the bad channels and sync bit names back to the meta file
imec.writeModifiedAPMeta();

% CAR the file and flush the unused sync bits to the cleaned_datasets folder
cleanedBinFile = 'Vinnie_20180817_All_cleaned.imec.ap.bin'
```


```text
cleanedBinFile = 'Vinnie_20180817_All_cleaned.imec.ap.bin'
```


```matlab
cleanedPath = npxutils.generatePath(subject, 'cleaned_datasets', dateStr, cleanedBinFile);

fprintf('Writing CAR version at %s\n', cleanedPath);
```


```text
Writing CAR version at /home/janke/work/npxutils/example-data/Vinnie/cleaned_datasets/2018-08-17/Vinnie_20180817_All_cleaned.imec.ap.bin
```


```matlab
extraMeta = struct;
extraMeta.run_clearUnusedSyncBits = true;
extraMeta.run_detectAndMarkStimArtifactWindows = true;
fnList = {
    % This is broken because NeuropixelExpt doesn't exist?
    % @NeuropixelExpt.DataClean.clearUnusedSyncBits
    @npxutils.dataprocess.commonAverageReference
    }';
cleanedDataDir = fullfile(myExampleDataDir, 'Vinnie', 'cleaned_datasets');
npxutils.io.rmrf(cleanedDataDir);
imec = imec.saveTransformedDataset(cleanedPath, ...
    'goodChannelsOnly', false, 'writeSyncSeparate', false, ...
    'transformAP', fnList, 'extraMeta', extraMeta);
```


```text
Writing AP meta file Vinnie_20180817_All_cleaned.imec.ap.meta
Writing AP bin file /home/janke/work/npxutils/example-data/Vinnie/cleaned_datasets/2018-08-17/Vinnie_20180817_All_cleaned.imec.ap.bin
Writing contents of Vinnie_20180817_All
Copying ap file 1 / 1: Vinnie_20180817_All :    100% [_______________]   
Warning: AP bin file size is not an integral number of samples, file data may not be fully copied, truncating nSamplesAP
Error using npxutils.ConcatenationInfo (line 56)
Mismatch between time shifts and total number of samples
Error in npxutils.ImecDataset/readInfo (line 313)
                this.concatenationInfoAP = npxutils.ConcatenationInfo(this, 'ap', metaAP);
Error in npxutils.ImecDataset (line 193)
                    this.readInfo();
Error in npxutils.ImecDataset.writeConcatenatedFileMatchGains (line 2763)
            imecOut = npxutils.ImecDataset(outFile, 'channelMap', imecList{1}.channelMapFile, 'sourceDatasets', cat(1, imecList{:}));
Error in npxutils.ImecDataset/saveTransformedDataset (line 1943)
            imecOut = npxutils.ImecDataset.writeConcatenatedFileMatchGains({this}, outPath, p.Results);
```


```matlab
% Symlink into ks directory for Kilosort's use
ksPath = npxutils.generatePath(subject, 'ks', dateStr, cleanedBinFile);
fprintf('Symlinking to %s\n', ksPath);
% And run Kilosort on the symlinked data set
imec = imec.symLinkAPIntoDirectory(ksPath);
fprintf('Running Kilosort2\n');
npxutils.runKilosort2(imec);

% Here's a little summary output
info.rawPath = imecFile;
info.cleanedPath = cleanedPath;
info.ksPath = ksPath;
info
```

