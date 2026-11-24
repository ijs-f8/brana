"""
Module for reading MCNP meshtally files and converting them into the HDF5 format.

Reading from MCNP tallies output (ASCII):
	`mf = MeshtallyFile.FromMCNP('my_file.t')`

Results are available as:
	`mf.results` -- python dictionary

Export to HDF5:
	`mf.ToHDF5('my_file.h5')`

Import from HDF5:
	`mf = MeshtallyFile.FromHDF5('my_file.h5')`

Warning:
	The last index of 0-th axis (energy spectra) is sum over energy bins.
"""

# SINCE PREVIOUS VERSION:
#	- Merging for multiple files.
#  TO DO:
#	- Fix tally indexing in MeshtallyFile.results: now the keys are 'Tally <tally_id>' for easier
# 		compatibility with HDF5.
#		Ideally we would have integer keys in python objects and str keys in HDF5 !
#	- Update help text at begining of module (begining of this file)
#	- Time bins
#	- Add function to find tallies by their name
#	- Add non-cartesian geometries
#	- Add function for plotting


import numpy as np
import h5py
import re
import os
import datetime
from . import __version__


def _split_list(list, splitpoints):
	"""
	Used by MeshtallyFile.SplitMeshtallies().
	"""
	sublists = []
	last_sp = 0
	for sp in splitpoints:
		sublists.append(list[last_sp:sp])
		last_sp = sp
	sublists.append(list[last_sp:])
	return sublists

def opt_merge_arrays(A:np.ndarray, B:np.ndarray,
					 sigmaA:np.ndarray, sigmaB:np.ndarray):
	"""
	Optimal merging is a weighted average which consideres
	uncertainties when asigning the weights.

	Arrays `A`, `B`, `sigmaA` and `sigmaB` must have the
	same shape.
	
	Returns tuple `(M, sigmaM)` where M is optimally merged
	array and sigmaM are standard deviations.
	"""

	u = np.square(sigmaA) + np.square(sigmaB)
	u_nonzero = np.not_equal(u, 0)

	_div = np.zeros(shape=u.shape)
	np.divide(
		np.multiply(np.square(sigmaA), np.square(sigmaB)),
		u,
		where=u_nonzero,
		out=_div)
	opt_sigma = np.sqrt(_div)

	opt_value = np.zeros(shape=u.shape)
	np.divide(
			np.multiply(A, np.square(sigmaB)) +
			np.multiply(B, np.square(sigmaA)),
		u,
		where=u_nonzero,
		out=opt_value)

	# We need to determine where these arrays have values zero.
	# If one of the two arrays has value zero for some voxel,
	# we cannot optimaly merge them.
	# We should instead take the non zero value of the two.
	# If both arrays have a zero value at some voxel,
	# we shoud just return (0,0) as the optimaly merged value.

	# if element i is:
	# A_i == 0 and B_i == 0 --> M_i = 0
	# A_i != 0 and B_i == 0 --> M_i = A_i
	# A_i == 0 and B_i != 0 --> M_i = B_i
	# A_i != 0 and B_i != 0 --> M_i = merge(A_i, B_i)

	zeroA = np.equal(A, 0)
	zeroB = np.equal(B, 0)
	AnonzeroBzero 	= np.logical_and( np.logical_not(zeroA), zeroB)
	BnonzeroAzero 	= np.logical_and( np.logical_not(zeroB), zeroA)
	ABnonzero		= np.logical_and( np.logical_not(zeroA), np.logical_not(zeroB))

	shape = A.shape
	zeros = np.zeros(shape)	
	
	M = ( np.where(ABnonzero, opt_value, zeros)
		+ np.where(AnonzeroBzero, A, zeros)
		+ np.where(BnonzeroAzero, B, zeros)
		+ zeros)
	sigmaM = ( np.where(ABnonzero, opt_sigma, zeros)
		+ np.where(AnonzeroBzero, sigmaA, zeros)
		+ np.where(BnonzeroAzero, sigmaB, zeros)
		+ zeros)
	
	return (M, sigmaM)

def _recursive_print_dict(_dict:dict, indent_level:int=0, indent_str:str='| '):
	_s = ''
	for k,v in _dict.items():
		if isinstance(v, dict):
			_s += f'{indent_level*indent_str}{k} (dict):\n'
			if not v: # empty dict
				_s += f'{(indent_level+1)*indent_str}empty dict\n'
			else:
				_s += _recursive_print_dict(v, indent_level=indent_level+1, indent_str=indent_str)
		elif isinstance(v, list):
			_s += f'{indent_level*indent_str}{k} (list)=\n'
			_s += f'{(indent_level+1)*indent_str}'
			_s += f'\n{(indent_level+1)*indent_str}'.join(v) + '\n'
		elif isinstance(v, tuple):
			_s += f'{indent_level*indent_str}{k} (tuple)=\n'
			_s += f'{(indent_level+1)*indent_str}'
			_s += f'\n{(indent_level+1)*indent_str}'.join(v) + '\n'
		else:
			_s += f'{indent_level*indent_str}{k}={v}\n'
	return _s

def py_to_hdf5(h5_root, name, value):
	"""
	Converts python objects to hdf5 objects.

	Python	 HDF5
	dict --> group
	str ---> attribute
	int ---> attribute
	float -> attribute
	bool --> attribute
	np.ndarray -> dataset
	list -------> dataset
	tuple ------> dataset

	Exceptions:
	list of str --> dataset of utf-8 string (*)
	dict with name results and int keys --> f'Tally {int}' keys (**)
	"""
	if isinstance(value, (str, int, float, bool)):		# attribute
		h5_root.attrs.create(name=name, data=value)
		return h5_root
	elif isinstance(value, dict):						# group
		h5_group = h5_root.create_group(name=name)
		for k,v in value.items():
			if name == 'results': # (**) convert key <number> (int) to 'Tally <number>' (str)
				py_to_hdf5(h5_group, f'Tally {k}', v)
			else:
				py_to_hdf5(h5_group, k, v)
		return h5_root
	elif isinstance(value, (np.ndarray, list, tuple)):	# dataset
		if name == 'log': # this should be a list of str (*)
			h5_root.create_dataset(name=name, data=value, dtype=h5py.string_dtype('utf-8'))
		else:
			h5_root.create_dataset(name=name, data=value)
		return h5_root
	else:
		raise TypeError(f"Unsupported type for variable {name=} {type(value)=}")

def hdf5_to_py(h5group):
	"""
	Converts HDF5 group to python dictionary.

	HDF5	 Python
	group -> dict
	attirbute -> python default (int, float or str)
	dataset ---> np.array

	Exceptions:
	dataset with ObjectDType -> list of str (*)
	group with name r'Tally <number>' --> dict with int keys (**)
	"""
	py_dict = {}
	if h5group.attrs.items():							# attributes
		for k,v in h5group.attrs.items():
			py_dict[k] = v
	for k,v in h5group.items():
		if isinstance(v, h5py.Dataset):					# datasets
			if isinstance(v.dtype, np.dtypes.ObjectDType):		# (*) this is a list of str
				py_dict[k] = np.asarray(v, dtype=str).tolist()
			else:
				py_dict[k] = np.array(v)
		elif isinstance(v, h5py.Group):					# groups
			m = re.match(r'Tally (\d+)', k)
			if m: # (**) convert key 'Tally <number>' (str) to <numer> (int)
				tid = int(m.group(1))
				py_dict[tid] = hdf5_to_py(v)
			else:
				py_dict[k] = hdf5_to_py(v)
	return py_dict


class MeshtallyFile:
	""" 
	`MeshtallyFile object` is used to parse, store, export and merge meshtally data.
		
	Initialization from MCNP meshtally output (ASCII):
		`meshtally = MeshtallyFile.FromMCNP('my_meshtally.t')`
		Parsing results:
		`meshtally.GetResults()`
		Results are stored in python dictionary`meshtally.results`.

	Exporting to HDF5 file:
		`meshtally.ToHDF5('my_meshtally.h5')`
	
	Loading from HDF5 file
		`meshtally = MeshtallyFile.FromHDF5('my_meshtally.h5')`
		
	A merged meshtally file can be created from two HDF5 files:
		`merged = MeshtallyFile.Merge('my_mt_A.h5', 'my_mt_B.h5')`
	"""

	def __init__(self):
		
		# temporary attributes
		self.file_text = None # this will be removed once GetResults() is run
		self.text_blocks = None # this will be removed once GetResults() is run

		# permanent attributes
		self.results = None
		self.metadata = {
			'name'		: '',	# used to be stored in self.name
			'header'	: [],	# used to be stored in self.header
			'log' : [], # --> HDF5 dataset
			# <date_time> <file_operation> <file_path>
			#   <date_time> = datetime.datetime.now()
			# eg.: ['2026-12-03 13:02:14.23531 read filename=/home/iter/111/a/ITER.t brana_version=8',
			#	    '2026-12-03 13:02:17.23531 export to filename=/home/iter/111/a/ITER.h5 brana_version=8',
			#	    '2026-12-03 17:02:17.23531 import from C:\my_files\ITER.h5',
			#	    '2026-12-03 20:02:17.23531 merged from C:\my_files\ITER_A.h5 and C:\my_files\ITER_A.h5']
			'parents' : {},		# --> HDF5 groups
			#	parent A : parentA.metadata, parent B : parentB.metadata
			'comment' : '',		# can be anything specified by the user - should be UTF-8 compatible
			
		}

		# state attributes
		self.is_converted_to_hdf5 = False
		self.is_read = False
		self.is_split_meshtallies = False
		self.parsed_results = False
		self.is_imported_from_hdf5 = False

	def __str__(self):
		indent_str = '| '
		s = '`MeshtallyFile` object with attributes:\n'
		s += _recursive_print_dict(
			{'metadata' : self.metadata},
			indent_level=1, indent_str=indent_str)
		s += f'{indent_str}results (dict): [tally_id, shape, tally_name]'
		s += f'\n{2*indent_str}'
		s += f'\n{2*indent_str}'.join([f'{k}\t{v['mean'].shape}\t{v['name']}' for k,v in self.results.items()])
		return s

	@staticmethod
	def FromMCNP(filename:str, print_progress:bool=False, skip_tallies:list=[]): # (***)
		"""
		Returns a `MeshtallyFile` object generated by reading a MCNP output file (ASCII).
		eg.:
		`mf = MeshtallyFile.FromMCNP('my_mcnp.t')`
		
		Results are available as
		`mf.results`

		Tally IDs provided in `skip_tallies` will be skipped when parsing with _get_results().
		This is ment for debuging large meshtally files.
		"""
		mt = MeshtallyFile()
		mt.metadata['name'] = '.'.join(filename.split('.')[:-1]) # remove suffix from filename
		file_handle = open(filename, 'r')
		mt.file_text = file_handle.read()
		file_handle.close()

		mt.is_read = True
		_abs_path = os.path.join(os.getcwd(), filename)
		mt.metadata['log'].append(f'{(datetime.datetime.now())} read file={_abs_path} brana_version={__version__}')

		mt._split_meshtallies()

		mt._get_results(print_progress=print_progress, skip_tallies=skip_tallies) # (***)

		return mt
		

	def _split_meshtallies(self):
		"""
		Function for method `MeshtallyFile.FromMCNP()` to parse MCNP output.
		"""
		tally_block_matches = re.finditer(r'\n\s*Mesh Tally Number\s+(\d+)', self.file_text)
		block_splitpoints = [tbm.start() for tbm in tally_block_matches]
		self.text_blocks = _split_list(self.file_text, block_splitpoints)
		self.metadata['header'] = re.split(r'\s*\n\s*', self.text_blocks[0])
		# self.is_split_meshtallies = True

	def _get_results(self, print_progress:bool=False, skip_tallies:list=[]): # (***)
		"""
		Function for method `MeshtallyFile.FromMCNP()` to parse MCNP output.
		
		Results are stored in `MeshtallyFile.results` attribute.
		"""
		#	Warning: should only be used for small files. -- is this true ??
		#	`save_asattr` if true, the results will be saved as -- not implemented ??
		#	  `MeshtallyFile.results`. Otherwise as above

		results = {}
		# results['info'] = self.header

		for tb in self.text_blocks[1:]:

			mt = MeshTally(tb, print_progress=print_progress, skip_tallies=skip_tallies) # (***)
			if mt.skip: # (***)
				print(f'Warning! Skipping tally {mt.id}, as it was specified in {skip_tallies=}.') # (***)
				continue  # (***)

			mt.ReshapeResults()

			results[mt.id] = { # see (**) # older: key is sotred as string to prevent issues in HDF5 -- try to fix later
				'name' : mt.name,
				'x bins' : mt.x_bins,
				'y bins' : mt.y_bins,
				'z bins' : mt.z_bins,
				'en bins' : mt.en_bins,
				'mean' : mt.mean_4d,
				'rel unc' : mt.rel_unc_4d,
			}

		self.file_text = 'removed by MeshtallyFile.GetResults()'
		self.text_blocks = 'removed by MeshtallyFile.GetResults()'

		self.results = results
		self.parsed_results = True
		
		# return results

	def ToHDF5(self, savename:str=None, print_progress=True): #  print_structure=False
		"""
		Export `MeshtallyFile` object to HDF5 file.
		"""
		
		if not self.results:
			raise ValueError(
			"""Nothing to export to HDF5 since `MeshtallyFile.results` are not available.
			Has the meshtally been read from MCNP output yet? Use either:
				`MEshtallyFile.FromMCNP(filename)`
				OR
				`MeshtallyFile.FromHDF5(filename)`""")

		if not savename:
			savename = f'{self.metadata['name']}.h5'

		_abs_path = os.path.join(os.getcwd(), savename)
		self.metadata['log'].append(f'{(datetime.datetime.now())} export file={_abs_path} brana_version={__version__}')

		h5file = h5py.File(savename, 'w')

		export_dict = {
			'metadata'				: self.metadata,
			'results'				: self.results,
			# 'is_converted_to_hdf5'	: self.is_converted_to_hdf5,
			# 'is_read'				: self.is_read,
			# 'is_split_meshtallies'	: self.is_split_meshtallies,
			# 'is_imported_from_hdf5'	: self.is_imported_from_hdf5,
		}

		py_to_hdf5(h5file, name='MeshtallyFile', value=export_dict)

		h5file.close()

		self.is_converted_to_hdf5 = True
		
	@staticmethod
	def FromHDF5(filename:str):
		"""
		Initialize `MeshtallyFile` object from HDF5 file.
		"""

		h5file = h5py.File(filename, mode='r')

		import_dict = hdf5_to_py(h5file)

		h5file.close()

		meshtally_file = MeshtallyFile()

		meshtally_file.metadata = import_dict['MeshtallyFile']['metadata']
		meshtally_file.results = import_dict['MeshtallyFile']['results']

		# meshtally_file.is_read = True
		meshtally_file.is_imported_from_hdf5 = True
		_abs_path = os.path.join(os.getcwd(), filename)
		meshtally_file.metadata['log'].append(
			f'{(datetime.datetime.now())} import file={_abs_path} brana_version={__version__}')

		return meshtally_file

	@staticmethod
	def Merge(files:list, _level_of_recursion:int = 0):
		"""
		List of meshtally files, h5 files or list of `MeshtallyFile` object instances.

		Returns merged `MeshtallyFile` object.
		"""

		print('| '*_level_of_recursion + f'{files=}')

		if len(files) == 0:
			raise ValueError("MeshtalFile.Merge() function takes a list of multiple meshtally files as input.")

		elif len(files) == 1:
			return MeshtallyFile.init_from_type(files[0])

		elif len(files) > 2: # solve recursion tree
			N = len(files)
			print('| '*_level_of_recursion + f'first half: {files[:N//2]}')
			print('| '*_level_of_recursion + f'second half: {files[N//2:]}')
			first_half = MeshtallyFile.Merge(files[:N//2], _level_of_recursion=_level_of_recursion+1)
			second_half = MeshtallyFile.Merge(files[N//2:], _level_of_recursion=_level_of_recursion+1)
			return MeshtallyFile.Merge(files=[first_half, second_half], _level_of_recursion=_level_of_recursion+1)

		else:
			pass
		
		fileA = MeshtallyFile.init_from_type(files[0])
		fileB = MeshtallyFile.init_from_type(files[1])

		meshtally_file = MeshtallyFile()

		meshtally_file.metadata['name'] = f'({fileA.metadata['name']} merge {fileB.metadata['name']})'
		# meshtally_file.metadata['info'] = f'merged from {fileA.metadata['name']} and {fileB.metadata['name']}'
		meshtally_file.metadata['header'] = ['see header of parents']
		# _abs_path_A = os.path.join(os.getcwd(), 
		# _abs_path_B = os.path.join(os.getcwd(), fileB._get_merge_name())
		meshtally_file.metadata['log'].append(
			f'{(datetime.datetime.now())} merge file_A={fileA._get_merge_name()} file_B={fileA._get_merge_name()} brana_version={__version__}')
		meshtally_file.metadata['parents'] = {
			'parent A' : fileA.metadata,
			'parent B' : fileB.metadata
		}

		results = {}
		for tally_id in fileA.results.keys():

			results[tally_id] = {}
			
			# copy all tally quantities except 'mean' and 'rel unc' from fileA
			for k in fileA.results[tally_id].keys():
				if k not in ['mean', 'rel unc']:
					try:
						results[tally_id][k] = fileA.results[tally_id][k]
					except KeyError:
						print(f'Tally {tally_id} is not abailable') # TODO fix

			# calculate optimally merged arrays
			merged, merged_std_dev = opt_merge_arrays(
				# mean values
				A=fileA.results[tally_id]['mean'],
				B=fileB.results[tally_id]['mean'],
				# standad deviation -- calcualted from rel uncertainty
				sigmaA=np.multiply(
					fileA.results[tally_id]['mean'],
					fileA.results[tally_id]['rel unc']),
				sigmaB=np.multiply(
					fileB.results[tally_id]['mean'],
					fileB.results[tally_id]['rel unc'])
			)
			
			results[tally_id]['mean'] = merged
			# where_non_zero = np.not_equal(merged, 0)
			results[tally_id]['rel unc'] = np.zeros(shape=merged.shape) # (****)
			np.divide( # (****)
				merged_std_dev, merged, # calculate rel unc form standard deviation # (****)
				out=results[tally_id]['rel unc'], # (****)
				where=merged != 0, # (****)
			) # (****)
			# results[tally_id]['rel unc'] = # CALCIULATE AND KEEP MCNP NOTATION FOR 0 PARTICLE RESULTS

		meshtally_file.results = results

		return meshtally_file

	@staticmethod
	def init_from_type(meshtally):
		if isinstance(meshtally, MeshtallyFile):
			return meshtally
		
		if isinstance(meshtally, (os.PathLike, str)):
			pass
		else:
			raise ValueError(f'Meshtally object type not supported. You provided {type(meshtally)=}.')

		if meshtally.endswith(('.h5', '.H5', '.hdf5', '.HDF5')):
			mt = MeshtallyFile.FromHDF5(meshtally)
			return mt
		else:
			mt = MeshtallyFile.FromMCNP(meshtally)
			return mt

	
	def _get_merge_name(self):
		for l in self.metadata['log']:
			m = re.match(r'read file=(.*) brana_version=', l)
			if m:
				return os.path.join(os.getcwd(), m.group(1))
		# else -- FromMCNP name not found. MeshtallyFile object was was made as merger.
		return self.metadata['name']
		# print('Warning! No birth name found in self.metadata["log"]. Searched for "read file=(.*) brana_version="')

	def SetComment(self, comment:str):
		"""
		User specified comment which will be stored in `MeshtallyFile.metadata['comment']`.

		Should be UTF-8 compatible.
		"""
		self.metadata['comment'] = comment

	def GetComment(self):
		"""
		Retrive user specified comment which is stored in `MeshtallyFile.metadata['comment']`.
		"""
		return self.metadata['comment']

	def Plot(self, tally_id, plt_axes, origin, ):
		raise NotImplementedError



class MeshTally:
	"""
	This is ment for parsing only. Not for storing tally data by the user.

	Used by `MeshtallyFile._get_results()`.
	"""

	def __init__(self, whole_tally_text:str, print_progress:bool=True, skip_tallies:list=[]): # (***)

		match_split = re.search(r'\s+Result\s+Rel Error\s*\n', whole_tally_text)
		self.info_text = whole_tally_text[:match_split.end()]
		
		self.id = int(re.search(r'Mesh Tally Number\s+(\d+)\s*\n', self.info_text).group(1))

		self.skip = False # (***)
		if self.id in skip_tallies: # (***)
			self.skip = True # (***)
			return None # (***)

		_name = re.search(r'Mesh Tally Number\s+\d+\s*\n\s*([^\n]+)', self.info_text).group(1)
		_empty_spaces = re.search(r'\s*$', _name)
		self.name = _name[:_empty_spaces.start()]

		if print_progress: print(f'\nmesh tally {self.id} : {self.name}')

		self.x_bins = np.array([float(s) for s in re.split(r'\s+', re.search(r'X direction:\s+(.*)\n', whole_tally_text).group(1))])
		self.y_bins = np.array([float(s) for s in re.split(r'\s+', re.search(r'Y direction:\s+(.*)\n', whole_tally_text).group(1))])
		self.z_bins = np.array([float(s) for s in re.split(r'\s+', re.search(r'Z direction:\s+(.*)\n', whole_tally_text).group(1))])
		self.en_bins = np.array([float(s) for s in re.split(r'\s+', re.search(r'Energy bin boundaries:\s+(.*)\n', whole_tally_text).group(1))])

		self.x_locs = np.average(np.array([self.x_bins[:-1], self.x_bins[1:]]), axis=0)
		self.y_locs = np.average(np.array([self.y_bins[:-1], self.y_bins[1:]]), axis=0)
		self.z_locs = np.average(np.array([self.z_bins[:-1], self.z_bins[1:]]), axis=0)

		rows = []
		results_lines = re.split(r'\n', whole_tally_text[match_split.end():])
		N = len(results_lines)
		for n,line in enumerate(results_lines):
			try:
				if n % int(N/100) == 0:
					if print_progress: print(f'  reading results ... {n/N*100:.0f} %', end='\r')
			except ZeroDivisionError: pass # this occurs for 1-voxel meshtallies
			# SKIP BLANK LINE
			if re.search(r'\S', line):
				# KEEP ONLY RESULT AND RELATIVE ERROR
				rows.append(  [ float(s) for s in re.split(r'\s+', line)[-2:] ]  )
			else:
				continue
		if print_progress: print(f'  reading results ... 100 %', end='\r')
		if print_progress: print('\n')

		self.results = np.array(rows)
		
	def ReshapeResults(self):
		"""
		Used by `MeshtallyFile._get_results()`.
		"""
		
		if len(self.en_bins) <= 2:
			en_lines = 1
		else:
			en_lines = len(self.en_bins)

		mesh_shape = (en_lines, len(self.x_bins)-1, len(self.y_bins)-1, len(self.z_bins)-1) # not considering energy bins
		meshtally_shape = (en_lines, mesh_shape[1], mesh_shape[2], mesh_shape[3], -1)  # not considering energy bins

		results_4d = self.results.reshape(meshtally_shape)

		self.mean_4d = results_4d[:,:,:,:,-2]
		self.rel_unc_4d =   results_4d[:,:,:,:,-1]


