"""
Module for handling MCNP (ASCII) output files.

Example:
	`of = OutputFile('my_file.o')`
	`of.GetResults()`
	`of.results` --> dictionary with tally_ids as keys and numpy arrays as values
"""


# TO DO:
#	- when saving to numpy, save to NPZ numpy dict file
#		* tallies will be saved as "tally_<tally_id>"=...
#		* perturbations as         "pert_<pert_id>_tally_<pert_tally_id>"=...
#	- a more efficient way to store tally text is needed:
#		* now each tally object stores whole tally ASCII text
#		* this means we store tally info twice:
#			OutputFile.mo_text and Tally.result_blocks_text
#		* for large output files this may be probelmatic
#		* for lost particles this shouldn't be done!
#			For large NPS lost particles can increase output file by > 10 GB !
#		* Maybe try using "mmap" ??
#	- is there a feature in MCNP to specify the hirearchy of tally bins?
#		* If this is the case this should also be considered
#		* the MCNP input part of the output file can also be read ;)
#		* this is decribed in section 6.3.3.7 (page 629 of MCNP 6.3 manual):
#			FDUSMET (Fred Died Under Some Mysterious Cicumstances Editing Tallies)
#	- we could create a global variable OUTPUT_FILES = {}
#		* a dictionary with ID keys and OutputFile object values
#		* a tally object would have a `parent` atttribute which would reference
#			the ID of its parent OutputFile
#		* we wouldn't need to pass the tally_text explicitly,
#			we would only pass the startpoint and endpoints and read the Outputfile.mo_text


import numpy as np
import re
import time
import os


# TALLY_PATTERN = r'1tally\s+(\d+)([\s|\S]+?)\n\n\n'
TALLY_PATTERN = r'1tally\s+(\d+)\s+nps\s*=\s*(\d+)([\s|\S]+?)tally\s+\d+\n' # groups: 1: tally_id, 2: nps, 3: tally_text
PERT_PATTERN = r'perturbation\s+(\d+)\.'
LOST_PATTERN = r'1\s+lost particle no.\s+(\d+)([\s|\S]+?)\n\n\n'

REAL_NUMBER_PATERN = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?" # REAL_NUMBER_PATERN = r'[-|+]?(\d*\.*\d*e?[\-|\+]?\d*(\((.*?)\))?)'

class Tally:
	"""
	Tally object is initiated wtith tally_id:
		tally = Tally(tally_id)
	MCNP out is read with ReadTally method:
		tally.ReadTally(MCNP_out_text)
	An array with tally results is retrived with the GetResults method:
		tally.GetResults()
	"""
	def __init__(self, tally_id:int, startpoint:int, endpoint:int):
		self.id = tally_id
		self.startpoint = startpoint
		self.endpoint = endpoint
		self.block_bins = {}	# dictionary by bin_type. Contains list of bins as they correspond to each block
		self.bins = {}	# 'bin_type' : [bins]	- dict of all tally bins. Keys are bin_type
		self.n_bins = {}	# number of tally bins
		self.results = None
		self.blocks_info_text = None
		# self.results = []
		self.bin_pattern = {
			'surface'	: r'surface\s+(\S+)',
			'cell'		: r'cell\s+(\S+)',
			'cosine'	: r'cosine bin\:\s+(\S+)\s+to\s+(\S+)'
		}
	
	def __str__(self):
		return f'Tally: id={self.id}'

	def ReadTally(self, tally_text:str):
		""" Splits `tally_text` into `blocks_text` and other info,
			which can be parsed by other functions such as `GetResults()` and
			`GetBins()`. """
		self.blocks_text = re.split(r'\n\s\n', tally_text) # maybe this is not such a good idea for large tally files ?
		# for lost particles this shouldn't be done! For large NPS lost particles can increase output file by > 10 GB !
		self.info_text = self.blocks_text[0]
		match_nps = re.search(r'nps =\s+(\d+)', self.info_text)
		self.nps = int(match_nps.group(1))
		self.result_blocks_text = self.blocks_text[1:]
		match_name = re.search(r'\+\s+([^\n]+)', self.info_text) # match everything after '+' exepct newline
		if match_name: self.name = match_name.group(1)
		else: self.name = None
		self.GetMultiplier()
		self.Flagging()

	def GetMultiplier(self):
		match_multiplier = re.search(r'this tally is all multiplied by\s+([\d|\.|E|e|\+|\-]+)', self.info_text)
		if match_multiplier:
			self.multiplier = float(match_multiplier.group(1))
		else:
			self.multiplier = 1.0

	def Flagging(self):
		match_flagging = re.search(r'flagging\s\S+\:\s+(\S+)', self.info_text)
		if match_flagging:
			self.flagging = match_flagging.group(1)
			return self.flagging
		else:
			self.flagging = False
			return self.flagging
	
	def GetBins(self, bin_type:str):
		block_bins = []
		bins = []
		#
		#	block bins
		for block in self.result_blocks_text:
			match = re.search(self.bin_pattern[bin_type], block)
			# print(match.__dir__())
			# print(match.groups())
			# print(match.group(1))
			if match:
				block_bins.append(list(match.groups()))
				bins.append(list(match.groups()))
		if len(block_bins) == 0:
			print(f'Warrning: No bins of type "{bin_type}" were found in tally {self.id}.')
			return None
		self.block_bins[bin_type] = block_bins
		#
		#	tally bins
		if bin_type not in list(self.bins.keys()):
			self.bins[bin_type] = []	# make new dictionary entry if it doesn't exist yet
		for b in bins:
			if not b in self.bins[bin_type]: self.bins[bin_type].append(b)
			else: pass
		#	number of bins by bin type
		for k,v in self.bins.items():
			self.n_bins[k] = len(v)
		#
		return self.bins
	
	def PrintBins(self, bin_type=None):
		if bin_type:
			self.GetBins(bin_type)
		for k, v in self.bins.items():
			print(f'{k} bins = {v}')

	def PrintBlocks(self, max_start_lines=8, max_end_lines=3):
		for i,block in enumerate(self.blocks_text):
			lines = block.split('\n')
			lines_short_space = [re.sub(r'\s{3}\s+', '    ', line) for line in lines]
			n = len(lines_short_space)
			print(f'\nBlock {i} has {n} lines:')
			max_print_lines = max_start_lines + max_end_lines
			if n > max_print_lines:
				# print shorthand
				for l in lines_short_space[:max_start_lines]:
					print(l)
				print('...')
				for l in lines_short_space[-max_end_lines:]:
					print(l)
			else:
				for l in lines_short_space: print(l)
			print('\n')

	def GetResults(self):	#, find_only_flagged=False):
		""" Parses `Tally.result_blocks_text` and
			returns the results as a numpy array. """
		tally_results = []
		blocks_info_text = []
		for block in self.result_blocks_text:
			"""if find_only_flagged == 'flagged':
				if re.search('flagged', block):
					pass
				else:
					continue
			elif find_only_flagged == 'unflagged':
				if re.search('flagged', block):
					continue
				else:
					pass
			else: pass"""
			lines = re.split(r'\n', block)
			block_results = []
			block_info_text = []
			for line in lines:
				if re.search(r'^\s{3}\s+[\d|\-]', line):
					_all = re.findall(r'\s+(\S+)', line)
					block_results.append([float(s) for s in _all])
					# print(f'result in line {line}')
				else:
					block_info_text.append(line)
			
			if len(block_results) > 0:
				tally_results.append(block_results)
			
			blocks_info_text.append(block_info_text)
		self.blocks_info_text = [''.join(ls) for ls in blocks_info_text]
		"""if find_only_flagged == 'flagged':
			self.flagged_results = np.array(tally_results)
		elif find_only_flagged == 'unflagged':
			self.unflagged_results = np.array(tally_results)
		else:"""
		# print(tally_results)
		self.results = np.array(tally_results)
		return np.array(tally_results)
	
	"""def GetFlaggedResults(self):
		return self.GetResults(find_only_flagged='flagged')"""

	"""def GetUnglaggedResults(self):
		return self.GetResults(find_only_flagged='unflagged')"""
	
	def ReshapeResults(self, n_primary_bins, n_secondary_bins, n_energy_bins):
		# primary bins:		cells / surfaces / ...
		# secondary bins:	cosine / ...
		#
		# returns shape (n_primary_bins, n_flaging_bins, n_secondary_bins, n_energy_bins, n_values)
		#	where:
		#	n_flaging_bins = 1 if no falgging or 2 if flagging
		#	n_values = 2 (mean and std_dev)
		if len(self.results) == 0:
			print('Error: Method Reshape() was called without results being read. Call method GetResults() first.')
			quit()
		#
		n_values = 2
		if self.flagging:
			n_flagging_bins = 2
		else:
			n_flagging_bins = 1
		#
		self.results = np.reshape(self.results, (n_primary_bins, n_flagging_bins, n_secondary_bins, n_energy_bins, n_values))
		#
		#	energy bins could be taken directly from self.results
		#
		return self.results

class LostParticle:
	def __init__(self, particle_number:int, startpoint:int, endpoint:int):
		""" `startpoint` and `endpoint` are positions of lost particle dedscription
		within the MCNP output text. """
		self.particle_number = particle_number
		self.startpoint = startpoint
		self.endpoint = endpoint
		self.xyz = None
		self.uvw_cos = None
		self.surface = None
		self.cell = None
		self.erg = None
		self.wgt = None
		self.tme = None
		self.n_cols = None
		self.dist = None


	
	def ReadLostPart(self, lost_part_text:str):
		m = re.search(rf'x,y,z coordinates:\s+({REAL_NUMBER_PATERN})\s+({REAL_NUMBER_PATERN})\s+({REAL_NUMBER_PATERN})', lost_part_text)
		self.xyz = np.array([float(m.group(i)) for i in range(1,6,2)])
		m = re.search(rf'u,v,w direction cosines:\s+({REAL_NUMBER_PATERN})\s+({REAL_NUMBER_PATERN})\s+({REAL_NUMBER_PATERN})', lost_part_text)
		self.uvw_cos = np.array([float(m.group(i)) for i in range(1,6,2)])
		# for i in range(10):
		# 	print(m.group(i))

	def __str__(self):
		return f'Lost particle number {self.particle_number}.\n  xyz={self.xyz}\n  uvw_cos={self.uvw_cos}'
		


class OutputFile:
	def __init__(self, filename:str):
		self.filename = filename
		self.name = '.'.join(filename.split('.')[:-1])
		self.time_stamps = {}
		self.comp_time = {}

		self.time_stamps['file_read_start'] = time.time()
		file_handle = open(self.filename, 'r')
		self.mo_text = file_handle.read()
		file_handle.close()
		self.time_stamps['file_read_end'] = time.time()
		self.comp_time['file_read'] = self.time_stamps['file_read_end'] - self.time_stamps['file_read_start']

		self.time_stamps['tallies_init_start'] = time.time()
		
		# perturbations will appear before regular tallies (not after tallies?)
		# first deal with perturbations, in order to avoid them when interfering with normal tallies
		
		#
		#	PERTURBATIONS
		#
		self.perturbations = {}
		found_perts = re.search(PERT_PATTERN, self.mo_text)
		if found_perts:
			# tally_end = list(pert_matches)[0].start()
			pert_matches = re.finditer(PERT_PATTERN, self.mo_text)
			pert_spans = []
			pert_ids = []
			for pm in pert_matches:
				pert_ids.append(int(pm.group(1)))
				pert_spans.append(pm.start())
			pert_spans.append(-1)
			pert_spans = [(p1, p2) for p1,p2 in zip(pert_spans[:-1], pert_spans[1:])]

			for ps, pert_id in zip(pert_spans, pert_ids):
				self.perturbations[pert_id] = {}
				for tmatch in re.finditer(TALLY_PATTERN, self.mo_text[ps[0]:ps[1]]):
					pert_tally_id = int(tmatch.group(1))
					self.perturbations[pert_id][pert_tally_id] = Tally(pert_tally_id,
										startpoint=tmatch.start() + ps[0], # these should be the global startpoint
										endpoint=tmatch.end() + ps[0])
			tally_end = pert_spans[0][0]
		else:
			tally_end = -1



		#
		#	LOST PARTICLES
		#
		self.lost_particles = []
		for match in re.finditer(LOST_PATTERN, self.mo_text):
			part_num = int(match.group(1))
			
			lp = LostParticle(part_num, match.start(), match.end())

			lp.ReadLostPart(self.mo_text[match.start():match.end()])
			
			self.lost_particles.append(lp)

			# self.lost_particles.append(LostParticle(part_num, match.start(), match.end()))


		#
		#	TALLIES
		#
		self.tallies = {}
		for match in re.finditer(TALLY_PATTERN, self.mo_text[:tally_end]): # only until perturbations start
			tally_id = int(match.group(1))
			self.tallies[tally_id] = Tally(tally_id,
								  startpoint=match.start(),
								  endpoint=match.end())
		self.time_stamps['tallies_init_end'] = time.time()
		self.comp_time['tallies_init'] = self.time_stamps['tallies_init_end'] - self.time_stamps['tallies_init_start']

		self.results = {}

	def GetTalliesText(self):
		""" This function is the most robust.
			It does not parse the tallies, just returns a dictionary
			with tally ids as keys and tally texts as values. """

		dict_of_tallies = {}
		for _id, t in self.tallies.items():
			s, e = t.startpoint, t.endpoint
			dict_of_tallies[_id] = self.mo_text[s:e]
		return dict_of_tallies

	def ReadTally(self, tally_id:int):
		""" Reads a specific tally and returns a Tally object. """

		self.time_stamps[f'tally_{tally_id}_read_start'] = time.time()
		startpoint = self.tallies[tally_id].startpoint
		endpoint = self.tallies[tally_id].endpoint - 3 # skip the triple newline at the end of the tally

		self.tallies[tally_id].ReadTally(tally_text = self.mo_text[startpoint : endpoint] )
		
		self.time_stamps[f'tally_{tally_id}_read_end'] = time.time()
		self.comp_time[f'tally_{tally_id}_read'] = self.time_stamps[f'tally_{tally_id}_read_end'] - self.time_stamps[f'tally_{tally_id}_read_start']

		return self.tallies[tally_id]
	
	def ReadPertTally(self, pert_id:int, tally_id:int):
		self.time_stamps[f'pert_{pert_id}_tally_{tally_id}_read_start'] = time.time()
		startpoint = self.perturbations[pert_id][tally_id].startpoint
		endpoint = self.perturbations[pert_id][tally_id].endpoint - 3 # skip the triple newline at the end of the tally

		self.perturbations[pert_id][tally_id].ReadTally(tally_text = self.mo_text[startpoint : endpoint] )
		
		self.time_stamps[f'pert_{pert_id}_tally_{tally_id}_read_end'] = time.time()
		self.comp_time[f'pert_{pert_id}_tally_{tally_id}_read'] = self.time_stamps[f'pert_{pert_id}_tally_{tally_id}_read_end'] - self.time_stamps[f'pert_{pert_id}_tally_{tally_id}_read_start']

		return self.perturbations[pert_id][tally_id]
		
	def GetResults(self):
		""" Reads all tally instances in OutputFile.
			Creates results array for each tally by calling Tally.ReadTally().
			Eg.: `of = OutputFile()` and `of.GetResults()`,
			then results can be read as `of.results` which is a dictionary with
			tally ids as keys and result arrays as values. """

		self.time_stamps['to_arrays_start'] = time.time()
		
		tally_ids = list(self.tallies.keys()) 
		for tid in tally_ids:
			self.ReadTally(tid)
			self.results[tid] = self.tallies[tid].GetResults()

		pert_ids = list(self.perturbations.keys())
		for pid in pert_ids:
			pert_tally_ids = list(self.perturbations[pid].keys())
			# print(pert_tally_ids)
			self.results[f'pert {pid}'] = {}
			for ptid in pert_tally_ids:
				self.ReadPertTally(pid, ptid)
				self.results[f'pert {pid}'][ptid] = self.perturbations[pid][ptid].GetResults()
		
		self.time_stamps['to_arrays_end'] = time.time()
		self.comp_time['to_arrays'] = self.time_stamps['to_arrays_end'] - self.time_stamps['to_arrays_start']

	def GetLostPosition(self):
		""" Returns (N,3) array for N lost particles.
			Second dimension is 3 coordinates (xyz). """

		xyz = np.array([lp.xyz  for lp in self.lost_particles])

		return xyz

	def GetLostCosine(self):
		""" Returns (N,3) array for N lost particles.
			Secpnd dimension is 3 direction cosines (u,v,w). """
		
		uvw_cos = np.array([lp.uvw_cos  for lp in self.lost_particles])
	
		return uvw_cos


	def ToHDF5(self, savename:str=None):
		self.time_stamps['to_hdf5_start'] = time.time()
		import h5py
		if not savename:
			savename = f'{self.name}.h5'
		h5file = h5py.File(savename, 'w')
		for k,v in self.results.items():
			tally_group = h5file.create_group(name=f'tally {k}')
			tally_group.create_dataset(name='results', data=v)
		h5file.close()
		self.time_stamps['to_hdf5_end'] = time.time()
		self.comp_time['to_hdf5'] = self.time_stamps['to_hdf5_end'] - self.time_stamps['to_hdf5_start']

	def ToNumpy(self, savedir:str=None): #  tally:int=None, savename:str=None, savedir:str=None):
		""" Converts tally to numpy array and saves as name_tallyID.npy
		If no tally_id is provided it will export all tallies to numpy arrays. """
		
		self.time_stamps['to_numpy_start'] = time.time()
		
		if savedir:
			# remove last slash in path if exists
			if savedir.endswith('/') or savedir.endswith('\\'):
				savedir = savedir[:-1]
			try: os.mkdir(savedir)
			except: FileExistsError
			for k,v in self.results.items():
				np.save(file=f'{savedir}/{k}.npy', arr=v)
		else:
			try: os.mkdir(self.name)
			except: FileExistsError
			for k,v in self.results.items():
				np.save(file=f'{self.name}/{k}.npy', arr=v)
		
		self.time_stamps['to_numpy_end'] = time.time()
		self.comp_time['to_numpy'] = self.time_stamps['to_numpy_end'] - self.time_stamps['to_numpy_start']

		"""

		if tally:
			if savename:
				if savedir:
					pass # export single tally to f'{savedir}/{savename}.npy'
				else:
					pass # export single tally to f'{savename}.npy'
			else:
				if savedir:
					pass # export single tally to f'{savedir}/{self.name}_{tally}.npy'
				else:
					pass # export single tally to f'{self.name}_{tally}.npy'
		else:
			if savename:
				if savedir:
					pass # export all tallies to f'{savedir}/{savename}.npy'
				else:
					pass # export all tallies to f'{savename}.npy'
			else:
				if savedir:
					pass # export all tallies to f'{savedir}/{self.name}_{tally}.npy'
				else:
					pass # export all tallies to f'{self.name}_{tally}.npy'"""
