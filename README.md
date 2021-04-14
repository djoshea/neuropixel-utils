# neuropixel-utils

[![View Neuropixel Utils on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/81238-neuropixel-utils)
[![View Documentation](https://img.shields.io/badge/docs-latest-blue)](https://djoshea.github.io/neuropixel-utils/)

Neuropixel Utils is a toolkit written in Matlab for manipulating datasets collected by [SpikeGLX](https://github.com/billkarsh/SpikeGLX) (e.g. `imec.ap.bin` files) and the results produced by [Kilosort](https://github.com/cortex-lab/KiloSort) / [Kilosort 2](https://github.com/MouseLand/Kilosort2/). Please note that some of this functionality is redundant with the tools found in the Cortex Lab's [spikes repository](https://github.com/cortex-lab/spikes), authored By Nick Steinmetz, Mush Okun, and others. Here, we prioritize an organized, easy to use, object-oriented approach to accessing, manipulating, and visualizing the data. This reduces the need to worry about metadata.

See full documentation at <https://djoshea.github.io/neuropixel-utils>.

Neuropixel Utils facilitates the following data processing steps:

- [Load](https://djoshea.github.io/neuropixel-utils/imec_dataset#constructing-a-neuropixelimecdataset) and [visualize](https://djoshea.github.io/neuropixel-utils/imec_dataset#plotting-specific-time-windows) raw neuropixel data from `imec.ap.bin` and `imec.lf.bin` files in Matlab
- [Write custom pre-processing functions](https://djoshea.github.io/neuropixel-utils/imec_dataset#building-a-preprocessing-pipeline) to apply to raw data either by writing a copy of the raw file or modifying it in place, optionally removing specific problematic time windows in the file
- [Concatenate multiple Imec data files](https://djoshea.github.io/neuropixel-utils/imec_dataset#concatenating-multiple-files-together) together while matching the amplifier gains
- [Run Kilosort/Kilosort2](https://djoshea.github.io/neuropixel-utils/kilosort#running-kilosort), and [load the results](https://djoshea.github.io/neuropixel-utils/kilosort#loading-kilosort-results) back into Matlab after manual inspection in Phy
- [Plot drift maps](https://djoshea.github.io/neuropixel-utils/analysis#plotting-drift-maps) using code adapted from the [spikes repository](https://github.com/cortex-lab/spikes)
- [Extract waveforms](https://djoshea.github.io/neuropixel-utils/waveforms#extracting-waveforms-via-kilosortdataset) for each cluster from the raw data, optionally cleaning the snippets by subtracting templates for other clusters spiking during the same time window
- [Visualize cluster electrical spiking images](https://djoshea.github.io/neuropixel-utils/analysis#plotting-electrical-images) in space and [cluster locations on the probe](https://djoshea.github.io/neuropixel-utils/analysis#plotting-cluster-centers-of-mass)
- [Determine trial boundaries](https://djoshea.github.io/neuropixel-utils/kilosort#segmenting-a-kilosort-dataset-into-trials) in the file, and efficiently [segment Kilosort results into individual trials](https://djoshea.github.io/neuropixel-utils/kilosort#kilosorttrialsegmenteddataset)

Neuropixel Utils was authored by [Daniel J O'Shea](http://djoshea.com) ([@djoshea](https://twitter.com/djoshea)) to facilitate precision artifact removal and careful inspection of raw data traces before running Kilosort, as well as post-hoc verification that the artifacts were removed successfully.

## Download and install

To get started, clone the repo:

```bash
git clone https://github.com/djoshea/neuropixel-utils.git
```

And add it to your path in Matlab:

```matlab
>> addpath('/path/to/neuropixel-utils')
```

## Requirements

Neuropixel Utils requires Matlab R2019b or later.
