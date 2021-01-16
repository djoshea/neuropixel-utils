# neuropixel-utils

[![View Neuropixel Utils on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/81238-neuropixel-utils)
[![View Documentation](https://img.shields.io/badge/docs-latest-blue)](https://djoshea.github.io/neuropixel-utils/)

Neuropixel Utils is a toolkit written in Matlab for manipulating datasets collected by [SpikeGLX](https://github.com/billkarsh/SpikeGLX) (e.g. `imec.ap.bin` files) and the results produced by [Kilosort](https://github.com/cortex-lab/KiloSort) / [Kilosort 2](https://github.com/MouseLand/Kilosort2/). Please note that some of this functionality is redundant with the tools found in the Cortex Lab's [spikes repository](https://github.com/cortex-lab/spikes), authored By Nick Steinmetz, Mush Okun, and others. Here, we prioritize an organized, easy to use, object-oriented approach to accessing, manipulating, and visualizing the data. This reduces the need to worry about metadata.

See full documentation at <https://djoshea.github.io/neuropixel-utils>.

Neuropixel Utils facilitates the following data processing steps:

- [Load](imec_dataset.md#constructing-a-neuropixelimecdataset) and [visualize](imec_dataset.md#plotting-specific-time-windows) raw neuropixel data from `imec.ap.bin` and `imec.lf.bin` files in Matlab
- [Write custom pre-processing functions](imec_dataset.md#building-a-preprocessing-pipeline) to apply to raw data either by writing a copy of the raw file or modifying it in place, optionally removing specific problematic time windows in the file
- [Concatenate multiple Imec data files](imec_dataset.md#concatenating-multiple-files-together) together while matching the amplifier gains
- [Run Kilosort/Kilosort2](kilosort.md#running-kilosort), and [load the results](kilosort.md#loading-kilosort-results) back into Matlab after manual inspection in Phy
- [Plot drift maps](analysis.md#plotting-drift-maps) using code adapted from the [spikes repository](https://github.com/cortex-lab/spikes)
- [Extract waveforms](waveforms.md#extracting-waveforms-via-kilosortdataset) for each cluster from the raw data, optionally cleaning the snippets by subtracting templates for other clusters spiking during the same time window
- [Visualize cluster electrical spiking images](analysis.md#plotting-electrical-images) in space and [cluster locations on the probe](analysis.md#plotting-cluster-centers-of-mass)
- [Determine trial boundaries](kilosort.md#segmenting-a-kilosort-dataset-into-trials) in the file, and efficiently [segment Kilosort results into individual trials](kilosort.md#kilosorttrialsegmenteddataset)

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
