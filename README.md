# brana

Brana is a collection of simple tools that improve the MCNP workflow. The name comes from a mountain in the Kamnik-Savinja Alps range. It is also a Slovenian word for a tool used to work on (farm) fields. 

There exist better MCNP tool-kits than brana.
One notable example is [F4Enix](https://github.com/Fusion4Energy/F4Enix) - it is much more versitile, has many features and is further along in development.
However, in a few aspects this package can be advantageous:
* Fewer dependencies (only h5py).
* Parsing of tallies and meshtallies as numpy arrays rather than pandas DataFrames can be easier to understand and may result in more direct mathematical operaitons (For example: direct 3D field gradient calculation). 
* Brana suports exporting of parsed results into HDF5 format. These can be imported back to python objects.
* Unlike in F4Enix pyvista is not integrated. However, it can be the most efficent way to plot 3D data, even from a NDarray.

## Installation

`pip install brana`

## Input

Flexible and tolerant parsing of MCNP input files. Enables extraction of geometry and limited editing.

## Output

Multi-dimensional numpy array based processing of tally results:
* Easily distinguish individual cell bins based on array axis.

Example:
```python
from brana import OutputFile

of = OutputFile('mcnp_out')

of.GetResults()
of.results # --> dictionary with results

of.ReadLostPart()
of.lost_particles # --> list of lost particles
```

## Meshtal

Multi-dimensional numpy array based processing of meshtally results:
* Simple numpy slice based plotting.
* Easy to plot in 3D with pyvista.

Features:
* Export `MeshtallyFile` object to HDF5 format and importing from HDF5 format -- does not yet support HDF5 format introduced in MCNP 6.3.
* Tacebility: Internal log keeps track of where files were created, imported, exported, modified, ...
* Merging: Merge multiple meshtally files with optimal merging which consideres individual voxel uncertinaties.

Limitations:
* Cylindircal geometry.

Example:
```python
from brana import MeshtallyFile as mtf

# read MCNP meshtally output file 
mf = mtf.FromMCNP('my_meshtally')

# check file
print(mf) # --> prints structure, header information, tally ids, array shapes, tally name

# access resuts 
mf.results # --> python dictionary with tally IDs as keys

# export to HDF5
mf_a = mtf.FromMCNP('meshtal_a')
mf_a.ToHFD5('meshtal_a.h5')

# import from HDF5
mf_a_import = mtf.FromHDF5('meshtal_a.h5')

# merging
mf_b = mtf.FromMCNP('meshtal_b')
merged = mtf.Merge([mf_a_import, mf_b, 'meshtal_c'])

```
