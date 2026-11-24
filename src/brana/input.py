"""
Module for parsing and editing MCNP (ASCII) input files.

"""

# TO DO:
#	- Deck.ExportGeometry() needs to parse also material M and MT cards, maybe even MX ?
#	- Add parsing for periodic surfaces in Deck.ExportGeometry()
#	- Check get_id_and_data_type() for heavy ion particle designators that are defined with numeric values.
#	- Is Deck.SetSettings() redundant? -- see implementation of AddCell and AddSurface
#	- When exporting geometry, also loop over surfaces to obtain all transformations to be exported :)
#		* You should also loop over transforamtions refrenced in FILL parameters of cells!
#	- Making a new Deck from scratch is clunky. It could have a more openMC-like approach
#	- REAL_NUMBER_PATTERN has yet to be implemented for all instnces in MCNP input
#		* it has to be more universal:
#		* it should match all the expresions that python recognizes as FLOAT variables
#	- Implement Parse() methods for Surfaces, Materials, Transformations, ...
#		OR at least make it so that a user can run Deck.Parse() without having to specify scope=['cells']
#		* Maybe the Deck.Parse method can chech wether the parse methods are implemented ?
#	- Implement loop over complement_cells in ExtractGeometry. Even if universe recursion is off!
#	- Could I run MCNP from python notebook
#		* if I specify the MCNP executable and build the MCNP environment ??
#		* That would be super cool

import re
import numpy as np
import copy



#
#	SEE MCNP 6.3 manual pages 3-5 (table of contents: chapter 5) for an overview of CARDS
#
SETTINGS_DATA_CARDS = [
	# geometry focused settings cards (5.5)
	'vol', 'area', 'u', 'lat', 'fill', 'uran', 
	# material / nuclear data specific settings cards (5.6) -- others are in MATERIAL_DATA_CARDS
	'mt', 'mx', 'mpn', 'otfdb', 'totnu', 'nonu', 'awtab', 'xs', 'void', 'mgopt', 'drxs',
	# physics focused settings cards (5.7)
	'mode', 'phys', 'act',
	'cut', 'elpt',
	'tmp', 'thtme', 'dbrc',
	'mphys', 'lca', 'lcb', 'lcc', 'lea', 'leb',
	'fmult', 'tropt', 'unc',
	'cosyp', 'cosy', 'bfld', 'bflcl', 'field',
	# - source specification settings cards (5.8) are in SOURCE_DATA_CARDS
	# - tally specification settings cards (5.9, 5.10 and 5.11) are in TALLY_DATA_CARDS
	# variance reduction settings cards (5.12)
	'imp',
	'var',
	'wwe',
	'wwt',
	'wwn',
	'wwp',
	'wwg', 'wwge',
	'wwgt',
	'mesh',
	'esplt',
	'tsplt',
	'exp',
	'vect'
	'fcl',
	'dxt',
	'dd', 'pd', 
	'dxc',
	'bbrem', 'pikmt',
	'spabi',
	'pwt',
	# problem termination, output control and misc cards (5.13)
	'nps', 'ctme', 'stop',
	'rand', 'print', 'talnp', 'prdmp', 'ptrac', 'mplot', 'histp', 'pio',
	'read', 'dbcn', 'lost', 'idum', 'rdum', 'files', 'disable'
]
TALLY_DATA_CARDS = [
	# tally cards (5.9)
	'f', 'fc', 'e', 't', 'c', 'fq', 'fm', 'de', 'df',
	'em', 'tm', 'cm', 'cf', 'sf', 'fs', 'sd', 'fu',
	'ft', 'tf', 'notrn', 
	# tally perturbations (5.10)
	'pert', 'kpert', 'ksen',
	# superimposed mesh tallies (5.11)
	'rmesh', 'cmesh', 'smesh', 'cora', 'corb', 'corc', 'ergsh', 'mshmf', 
	'xmesh', 'xmesh',
	'fmesh', 'spdtl'
]
SOURCE_DATA_CARDS = ['sdef', 'si', 'sp', 'sb', 'ds', 'sc',
	'ssw', 'ssr', 'kcode', 'ksrc', 'kopts', 'hsrc', 'burn']
	# 'source', 'srcdx'] -- these are subrutines, not input cards
TRANSFORMATION_DATA_CARDS = ['tr', '*tr']
MATERIAL_DATA_CARDS = ['m'] # includes only cards which make the material object (not mt, mx, mpn, ...)

BASE_TALLY_CARDS = ['f']
BASE_SOURCE_CARDS = ['sdef', 'ssw', 'ssr', 'kcode', 'ksrc',
	'kopts', 'hsrc', 'burn']

REAL_NUMBER_PATERN = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?" # REAL_NUMBER_PATERN = r'[-|+]?(\d*\.*\d*e?[\-|\+]?\d*(\((.*?)\))?)'



def previous_card(line:str):
	# checks if line describes previous card
	#   Line is multiline description --> return TRUE
	#   Line is new entity --> return FALSE 
	return bool(re.match(r'\s{5}', line))


def comment(line:str):
	# is line a comment?
	# if there are more than 4 spaces it is not a comment
	return bool(re.match(r'\s{0,4}c(?!\S)', line, flags=re.IGNORECASE)) # (?!\S) --> 'c' should not be followed by non-whitespace character

# def seperate_id_and_definition(initial_line:str):
# 	match_id = re.search(r'\s{,4}(\**[a-z0-9]+)', initial_line, flags=re.IGNORECASE)
# 	_id = match_id.group(1).lower()
# 	endpoint = match_id.end()
# 	definition = initial_line[endpoint:]
# 	return _id, definition


def get_id_and_data_type(initial_line:str, data_card:bool):
	""" Returns tuple (id, data_type).
		`id` is intiger.
		`data_type` is `None` unless `data_card` is set to `True`.
		`data_type` will be returend as lowered string.
		Eg. m | tr | *tr | sdef | fmesh | fm | si | sp | ds | ... 
		"""
	match_id = re.search(r'\s{,4}([\*|\+]?([a-z]*)([0-9]*)\:*[a-z]*)', initial_line, flags=re.IGNORECASE)
	# print(f'g0: {match_id.group(0)} | g1: {match_id.group(1)} | g2: {match_id.group(2)} | g3: {match_id.group(3)}')
	if data_card:
		if match_id.group(3) == '':
			_id = ''
		else:
			_id = int(match_id.group(3))
		data_type = match_id.group(2).lower() # this will be used to sort data cards
		data_key = match_id.group(1).lower() # this will be used to set a specific key in dictionary
		return (_id, data_type, data_key)
	else:
		_id = int(match_id.group(3))
		return (_id, None, None)


def seperate_definition_and_comment(line:str):
	m = re.search(r'\$', line)
	if m:
		definition = line[:m.start()]
		comment = line[m.end():]
		return definition, comment
	else: return line, ''

# def seperate_line_comment(line:str):
# 	m = re.search(r'c', line, flags=re.IGNORECASE)
# 	comment = line[m.end():]
# 	return comment


def data_cards_types(deck_lines:list):
	""" Returns list of card type initializations. """

	types = []
	ids = []
	for i,line in enumerate(deck_lines):
		if comment(line):
			types.append('comment')
			ids.append(0)
			continue
		if previous_card(line):
			types.append('prev_crd')
			ids.append(0)
			continue
		

		m_card = re.search(r'^\s{0,4}(\**([a-z]+))(\d*)', line, re.IGNORECASE)
		if not m_card:
			raise ValueError(f'Unrecognized data entry on data line {i}: {line}')

		# print(f'{i}: {m_card.group(1)} | {m_card.group(2)} | {m_card.group(3)}')

		types.append(m_card.group(2).lower())
		_id = m_card.group(3)
		if _id:
			ids.append(int(_id))
		else:
			ids.append(0)

		# NOW WE CAN ASSIGN EACH COMMENT AND PREV. CARD TO AN ID
		# 

	return np.array(types), np.array(ids)


def is_splitpoint_char(char):
	""" Function for function split_into_lines(),
	which is used by funciton MCNP_format_region(). """

	if char in [' ', '(', ')']: # removed ':' as splitpoint char, since it is also used in cell parameters. Eg.: IMP:N = 1.0
		return True
	else:
		return False


def split_into_lines(string:str, max_len:int=75):
	""" Returns a list of lines, split at max_len. """

	s = string
	list_of_lines = []
	pattern=r'[\+|\-|\s][\d|\w]+' # The format for card definition entities is \+, \- or \s before multiple digits or word characters
	# the string will be split if there is a splitpoint character between entities

	# check if entity lengths are too long
	margin = 5
	all_entities = re.findall(pattern, string)
	lengths = [len(entity) for entity in all_entities]
	if max(lengths) >= max_len - margin:
		print("Error: Length of entities within margin of max_len.")
		return None

	while len(s) > max_len:
		for i,c in enumerate(s):
			if i+1 > max_len:
				list_of_lines.append(s[:last_possible_splitpoint])
				s = s[last_possible_splitpoint:]
				break
			elif is_splitpoint_char(c):
				last_possible_splitpoint = i
				
	if len(s) != 0:
		# write also the left-overs
		list_of_lines.append(s)
	
	return list_of_lines


def split_into_lines_MCNP(string:str, max_len:int=75, n_spaces:int=5):
	""" Returns a list of lines, split at max_len, with appended 5 spaces. """
	
	list_of_lines = split_into_lines(string[n_spaces:], max_len)

	list_of_lines[0] = string[:n_spaces] + list_of_lines[0]

	list_of_lines = [' '*n_spaces + l if n != 0 else l for n,l in enumerate(list_of_lines)]

	return list_of_lines



class Card:
	def __init__(self, _id:str, _type:str, definition0:str='', comment0:str='', header:list=[], subtype:str=''):
		self.type = _type # cell | surface | material | transformation | tally | source | settings
		self.id = _id
		self.subtype = subtype
		# for data cards: sdef, ds, si, fc, fmesh, nps, rand, ...
		# for surface cards: px, s, so, gq, cz, ...
		self.definition = np.array([definition0])
		self.comments = np.array([comment0])
		self.header = header
		# self.delete_comments = delete_comments

	def AppendLines(self, definitions:list=[], comments:list=[]):
		if len(definitions) != len(comments):
			raise ValueError('definitions and comments should be lists of the same length.')
		self.definition = np.concat((self.definition, np.array(definitions)))
		self.comments = np.concat((self.comments, np.array(comments)))

	""" def __str__(self):
		if self.delete_comments:
			return f'{self.name}{self.definition}'
		elif self.prepended_comment:
			return f'{self.prepended_comment}{self.name}{self.definition}'
		else:
			return f'{self.name}{self.definition}' """

	def __str__(self):

		# _out = self.header # f'{self.definition[0]} {self.comments[0]}\n'
		_out = []

		for _d,_c in zip(self.definition, self.comments):

			# print(f'{d=}    {c=}')

			if _d and _c: # definition with inline comment
				_out.append( _d + '$' + _c )
			elif _d and not _c: # definition without comment
				_out.append( _d )
			elif _c and not _d: # comment line
				# if comment(_c): # this is line comment
				_out.append( _c )
				# else:          # this used to be inline comment - now it is a full line comment
				# 	_out.append('c' + _c)
			else:
				raise ValueError('Empty entry in definitions and comments')
				# _out.append('c')
		return '\n'.join(list(self.header) + list(_out)) #  + '\n' + '\n'.join(_out)

	def __repr__(self):

		def _max_width(str_list:list):
			if len(str_list) == 0:
				return 0
			else:
				return max([len(_s) for _s in str_list])	

		wd = _max_width(self.definition)
		wc = _max_width(self.comments)
		wh = _max_width(self.header)
		# print(wd, wc)
		_out = f'Card object, type={self.type}, id={self.id}, subtype={self.subtype}\nheader=\n{"-"*(wh+1)}\n|{"\n|".join(self.header)}\n{"-"*(wh+1)}\ndefinitions= | comments=\n{"-"*(wd+wc+2)}\n'
		def_com_lines = []
		for _d, _c in zip(self.definition, self.comments):
			def_com_lines.append(f'|{_d:{"_"}<{wd}}|{_c:{"_"}<{wc}}')
		return _out + '\n'.join(def_com_lines) + f'\n{"-"*(wd+wc+2)}'

	def IntoComment(self):
		""" Moves entire definition into comment section. """
		self.comments = [f'c {_d} {_c}' for _d,_c in zip(self.definition, self.comments)]
		self.definition = ['' for _c in self.comments]

	def Raw(self):
		""" Returns entire card definition as single text line, with single space separations. """
		return re.sub(r'\s+', ' ', ' '.join(self.definition))
	
	def NewDefinition(self, new_definition:str):
		""" Assigns new definition to card. Accepts one-line definition.
			This also overwrites old comments but leaves the header.
			5 spaces will be added before each line. """
		new_def_lines = split_into_lines(new_definition[5:], max_len=75) # update to use split_into_lines_MCNP
		new_def_lines[0] = new_definition[:5] + new_def_lines[0]
		for i,ndl in enumerate(new_def_lines):
			if i == 0: continue
			new_def_lines[i] = ' '*5 + ndl
		self.definition = np.array(new_def_lines)
		self.comments = np.full(shape=(len(new_def_lines)), fill_value='')

	def NewDefinitionLines(self, new_def_lines:list):
		""" Assigns new definition to card. Accepts one-line definition. This also overwrites old comments but leaves the header. """
		self.definition = np.array(new_def_lines)
		self.comments = np.full(shape=(len(new_def_lines)), fill_value='')


class Cell(Card):
	""" Cell object has the following additional attributes
		(in addition to Card attributes):
		`material` Int, reference to `Mn` card, is 0 if void
		`density` Float or None, weight density in g/cm3
		`parameters` Dict with Str keys and values
		`region` Str, reference to surfaces
		`complement_cells` List of Int, cells referenced with complement operator (`#`)
		`n_surfaces` Int, number of surfaces in region
		`surfaces` List of Int, unique surfaces referenced by Cell
		
		returns `universe` of cell (same as `int(Cell.parameters['u'])`)"""

	def Parse(self):
		""" This function reads all materials, surfaces, transformations, universes, ... referenecd by this cell.
			Returns universe ID. """
		
		#	THIS IS NOT YET COMPLETE!
		#	
		#	A lot of things are not yet acocunted for.
		#	See MCNP 6.3 MANUAL, pages 257-- 259
		#	Particularly the LUKE n BUT ... parameter is not considered yet
		#	Also the EXT parameter is not considered
		#	CAUTION: IS "IMP: N=1.0" A VALID KWARG ?? -- THIS IS NOT YET ACCOUNTED FOR
		
		_def = self.Raw()

		# Skip ID
		m_id = re.match(r'\s*(\d+)', _def)
		_def = _def[m_id.end():]

		# Material
		m_mat = re.match(r'\s*(\d+)', _def)
		self.material = int(m_mat.group(1))
		_def = _def[m_mat.end():]
		
		# Density
		if self.material == 0:
			# print('this cell is void')
			self.density = None
		else:
			m_dens = re.match(r'\s*(' + REAL_NUMBER_PATERN + r')', _def)
			self.density = float(m_dens.group(1))
			_def = _def[m_dens.end():]
		# print(f'Density = {cell_dens}')

		# Parameters
		self.parameters = {}
		# m_p = re.search(r'([a-z]+\:?[a-z]*)\s*\=\s*(\d*\.*\d*e?[\-|\+]?\d*)', _def, flags=re.IGNORECASE)
		# while m_p:
		# 	_key = m_p.group(1).lower()
		# 	_value = m_p.group(2).lower() # some are floats, some ints, some are strings
		# 	# if _key in ['trcl', 'u', 'fill', ]
		# 	self.parameters[_key] = _value
		# 	_def = _def[:m_p.start()] + _def[m_p.end():]
		# 	m_p = re.search(r'([a-z]+\:?[a-z]*)\s*\=\s*(\d*\.*\d*e?[\-|\+]?\d*)', _def, flags=re.IGNORECASE)
		m_p = re.search(r'(\*?[a-z]+\:?[a-z]*)\s*\=\s*(\d*\.*\d*e?[\-|\+]?\d*\s*(\((.*?)\))?)', _def, flags=re.IGNORECASE) # r'([a-z]+\:?[a-z]*)\s*\=\s*(\-?\d*\.*\d*[a-z]*)' # r'(\*?[a-z]+\:?[a-z]*)\s*\=\s*(\d*\.*\d*\s*(\((.*?)\))?)'
		while m_p:
			_key = m_p.group(1).lower()
			_value = m_p.group(2).lower() # some are floats, some ints, some are strings
			# if _key in ['trcl', 'u', 'fill', ]
			self.parameters[_key] = _value
			# print(f'Cell kwarg: {_key} = {_value}')
			_def = _def[:m_p.start()] + _def[m_p.end():]
			# print(f'{_def = }')

			# print(f'{m_p.group(1)}')
			# print(f'{m_p.group(2)}')
			# print(f'{m_p.group(3)}')
			# print(f'{m_p.group(4)}')
			# print(f'{m_p.group(1)} | {m_p.group(2)}')

			m_p = re.search(r'(\*?[a-z]+\:?[a-z]*)\s*\=\s*(\d*\.*\d*e?[\-|\+]?\d*\s*(\((.*?)\))?)', _def, flags=re.IGNORECASE)

		# Region
		self.region = re.sub(r'\s+', ' ', _def)

		# Cells in complement
		self.complement_cells = [] # list of cells in the complement (surfaces are not considered here! -- they are defined in region and will be saved in list of surfaces)

		m_c_c = re.search(r'(?<=\#)\s*(\d+)', _def)
		while m_c_c: 
			# print(f'Referenced cell {m_c_c.group(1)} (complement)')
			self.complement_cells.append(int(m_c_c.group(1)))
			_def = _def[:m_c_c.start()] + _def[m_c_c.end():]
			m_c_c = re.search(r'(?<=\#)\s*(\d+)', _def)

		ref_surfaces = [int(m.group(1)) for m in re.finditer(r'(\d+)\.?\d*', _def)]

		self.n_surfaces = len(ref_surfaces)
		self.surfaces = [int(s) for s in np.unique(ref_surfaces)]

		#
		#	RETURN UNIVERSE, BECAUSE THIS IS A PARENT STRUCTURE
		#

		try:
			universe = int(re.match(r'\d+', self.parameters['u']).group()) # this is done because universe value of universe parameter is string
		except KeyError:
			universe = 0

		return universe

	def Update(self, keep_header:bool=True):
		""" Update cell definition based on attributes. """
		_new_def = f'{self.id} {self.material} '
		if self.material != 0:
			_new_def += str(self.density) + ' '
		_new_def += self.region + ' '
		for k,v in self.parameters.items():
			_new_def += f'{k}={v} '
		_new_def = re.sub(r'\s+', ' ', _new_def)
		
		self.NewDefinition(_new_def)
		
		if not keep_header:
			self.header = []

	def GetFill(self):
		""" Returns `(fill_type, fill_id, fill_transformation)`.
		    `fill_type` is `fill` | `*fill`
			`fill_id` is intiger
			`fill_transformation` is all that remains """
		if 'fill' in self.parameters.keys():
			fill_type = 'fill'
		elif '*fill' in self.parameters.keys():
			fill_type = '*fill'
		else:
			return None
		
		fill_parameters = self.parameters[fill_type]

		mtch = re.match(r'\s*(\d+)\s*', fill_parameters)
		fill_id = int(mtch.group(1))
		fill_transformation = fill_parameters[mtch.end():]

		return (fill_type, fill_id, fill_transformation)

	def Clone(self, new_id:int=None):
		""" Returns a copy of the Card.
			`new_id` : int, If specified, the cloned card will be assigend
			a new ID. This will also trigger an update of the card definition.
			Warning: The card has the same ID as the source Card.
			Card is not assigned to a Deck.
			If card is assigned to the same Deck, the ID should be modified. """
		
		new_card = Cell(_id='none', _type='none')
		for attr_name, attr_value in self.__dict__.items():
			setattr(new_card, copy.deepcopy(attr_name), copy.deepcopy(attr_value))
		if new_id:
			new_card.id = new_id
			new_card.Update()
		return new_card

class Surface(Card):
	""" Surface object has the following additional attributes
		(in addition to Card attributes):
		`boundary_type` Str, either: 'vacuum' (default), 'white' (+ID) or 'reflective' (*ID)
		`transformation` Int, reference to `TRn` card
		`peiodic` Int, ID of corresponding periodic surface
		`subtype` Str, type of surface: `px`, `s`, `gq`, `cz`, ...
		`parameters` List of Str, parameters that define surface """

	def Parse(self):
		""" This method reads the surface type, its parameters and a potencial tranformation. """

		_def = self.Raw()

		# Parse ID and boundary_type
		m_id = re.match(r'\s*([\*|\+]?)(\d+)', _def)
		
		if m_id.group(1) == '*':	# a star (*) before the ID - reflective boundary surface
			self.boundary_type = 'reflective'
		elif m_id.group(1) == '+':	# a plus (+) before the ID - white boundary surface
			self.boundary_type = 'white'
		else:
			self.boundary_type = 'vacuum'

		_def = _def[m_id.end():]

		# Transformation or Periodic surface
		m_tr_pr = re.match(r'\s*(\-?)(\d+)', _def)

		if m_tr_pr:
			if m_tr_pr.group(1) == '-':
				self.periodic = int(m_tr_pr.group(2))
				self.transformation = None
			else:
				self.transformation = int(m_tr_pr.group(2))
				self.periodic = None
			_def = _def[m_tr_pr.end():]
		else:
			self.transformation = None
			self.periodic = None

		# Surface type
		m_st = re.match(r'\s*([a-z]+)', _def, flags=re.IGNORECASE)
		self.subtype = m_st.group(1)

		# Parameters
		self.parameters = [m.group() for m in re.finditer(REAL_NUMBER_PATERN, _def)]
		# --> this is a list of strings
		# otherwise we would need to specify some precission when running Surface.Update() 
			

	def Update(self, keep_header:bool=True):
		""" Update surface definition based on attributes: `id`, `transformation`, `subtype`, `parameters`. """

		if self.boundary_type == 'vacuum':
			_id = f'{self.id}'
		elif self.boundary_type == 'reflective':
			_id = f'*{self.id}'
		elif self.boundary_type == 'white':
			_id = f'+{self.id}'
		else:
			raise ValueError(f'Invalid surface `boundary_type`. Must be either `vacuum`, `reflective` or `white`. Provided value was `{self.boundary_type}`.')

		if self.transformation:
			_trans_period_and_type = f'{self.transformation} {self.subtype}'
		elif self.periodic:
			_trans_period_and_type = f'-{self.periodic} {self.subtype}'
		else:
			_trans_period_and_type = f'{self.subtype}'

		_new_def = f'{_id} {_trans_period_and_type} '
		
		_new_def += ' '.join(self.parameters)
		
		self.NewDefinition(_new_def)
		
		if not keep_header:
			self.header = []


# The Data Card should never be used, instead Material, Transforamtion, Source, ... should be used
# class Data(Card):
# 	""" This encompases Material, Transformation, Settings, Tally and Source cards."""
# 	pass

class Transformation(Card):
	pass

class Material(Card):
	pass

class Source(Card):
	pass

class Tally(Card):
	pass

class Settings(Card):
	pass



class Deck:

	def new(title:str=None):
		""" Returns a new blank Deck. """

		d = Deck(filename='', is_blank=True)

		if title: d.title = title

		return d

	def __init__(self, filename:str, is_blank=False, enconding='ansi'):
		
		is_clone = False

		if is_blank:

			self.filename = filename
			self.title = ''
			self.cells = {}
			self.cell_comment_tail = []
			self.surfaces = {}
			self.surf_comment_tail = []
			self.materials = {}
			self.transformations = {}
			self.source = {}
			self.tallies = {}
			self.settings = {}
			self.data_comment_tail = []
			self.universes = {}
			self.parsed = []
		
		else:

			self.filename = filename
			
			fh = open(filename, 'r', encoding=enconding)
			try: everything = fh.read()
			except UnicodeDecodeError:
				raise UnicodeDecodeError('Try changing enocding.')
			
			cells, surfaces, other, *_ = re.split(r'\n\s*\n', everything)

			self.title = cells.split('\n')[0]
			
			cell_lines = cells.split('\n')[1:]  # skip first line - title
			surf_lines = surfaces.split('\n')
			other_lines = other.split('\n')

			
			
			# self.cell_cards = cell_cards
			# self.surf_cards = surf_cards


			cell_cards, self.cell_comment_tail = self.ReadCards(cell_lines, 'cells')
			self.cells = {c.id : c for c in cell_cards}

			surf_cards, self.surf_comment_tail = self.ReadCards(surf_lines, 'surfaces')
			self.surfaces = {s.id : s for s in surf_cards}

			data_cards, self.data_comment_tail = self.ReadCards(other_lines, 'data')

			self.materials = {_d.id : _d for _d in data_cards if _d.type == 'material'}
			
			self.transformations = {_d.id : _d for _d in data_cards if _d.type == 'transformation'}

			self.source = {f'{_d.subtype}' : _d for _d in data_cards if _d.type == 'source'}

			self.tallies = {f'{_d.subtype}' : _d for _d in data_cards if _d.type == 'tally'}

			self.settings = {f'{_d.subtype}' : _d for _d in data_cards if _d.type == 'settings'}


			self.universes = {} # dict whose keys are universe IDs, values are cell IDs
			# Eg.:	cells 1-100 are universe 0 (base universe)
			#		cells 1000-2000 fill universe 1, which fills cell 10
			#		cells 10,000-20,000 fill universe 111, which fills cell 1111
			#		:
			#		universes =	{	0 : [1-100],
			#						1 : [1000-2000],
			#						111 : [10,000-20,000]}
			#		... fill information is stored in individual cell
			#
			#		This is just a dictionary which will alow us to find the necessary cells in order to extract geometry

			fh.close()

			self.parsed = []

	def __str__(self):

		n_cells = len(list(self.cells.keys()))
		n_surfaces = len(list(self.surfaces.keys()))
		n_materials = len(list(self.materials.keys()))
		n_tallies = len(list(self.tallies.keys()))
		n_source = len(list(self.source.keys()))

		out = f'MCNP input deck with {n_cells} cells, {n_surfaces} surfaces and {n_materials} materials.\nThere are {n_tallies} tally cards and {n_source} source cards.\nThe settings cards are:'

		for k,sc in self.settings.items():
			if sc.definition:
				out += f'\n  {sc.Raw()}'

		return out + '\n'

	def ReadCards(self, lines:list, deck_part:str):
		""" This is the main part of the Deck parser.
			This method loops over a part of the deck,
			specified by `deck_part`: `cells`, `surfaces` or `data`.
			The type of data card is determined inside the loop,
			either `material`, `transformation`, `source`, `tally`
			or `settings`. """
		
		cards = []
		stored_comments = []

		card_object = {
			'cell' : Cell,
			'surface' : Surface,
			'material' : Material,
			'transformation' : Transformation,
			'source'	: Source,
			'tally'	: Tally,
			'settings'	: Settings,
		}

		if deck_part == 'data':
			data_card = True
		else:
			data_card = False

		for i,l in enumerate(lines):
			
			#
			#	IN THIS LOOP LINES OF THE MCNP INPUT ARE READ
			#
			
			# if first line starts with
			#   4 or less spaces AND
			#   doesn't start with a letter
			# it's a cell or surface definition
			# else: it's a comment or a definition for previous entity

			if previous_card(l): # append line to definition of previous card
				
				# first append comments to previous cards
				_d = ['' for sc in stored_comments]
				cards[-1].AppendLines(definitions=_d, comments=stored_comments)
				stored_comments = []

				# print(f'From previous card: {l}')
				d, c = seperate_definition_and_comment(l)
				cards[-1].AppendLines(definitions=[d], comments=[c])


			elif comment(l):
				c = l # seperate_line_comment(l) # 'c'.join(l.split('c')[1:]) # remove first c in comment
				stored_comments.append(c)


			else: # card is definition of new entity
				_id, data_card_type, data_card_key = get_id_and_data_type(l, data_card)
				d, c = seperate_definition_and_comment(l)

				if data_card_type:
					#
					#	SEE MCNP 6.3 manual pages 3-5 (chapter 5 table of contents) for an overview of CARDS
					#
					if data_card_type in MATERIAL_DATA_CARDS:
						_type = 'material'
					elif data_card_type in TRANSFORMATION_DATA_CARDS:
						_type = 'transformation'
					elif data_card_type in SOURCE_DATA_CARDS:
						_type = 'source'
					elif data_card_type in TALLY_DATA_CARDS:
						_type = 'tally'
					elif data_card_type in SETTINGS_DATA_CARDS:
						_type = 'settings'
					else:
					
						raise NotImplementedError(f'Unknown data_card_type: "{data_card_key=}". Please implement - sort into appropriate card type: material | transformation | source | tally | settings.')

				elif deck_part == 'cells':
					_type = 'cell'
				elif deck_part == 'surfaces':
					_type = 'surface'
				else:
					if re.match(r'\s*', l): # this is an end of file empty line
						pass
					else:
						raise NotImplementedError(f'Unknown deck_part "{deck_part}". {data_card_type=}. {data_card_key=}. Problem on line {i} in deck_part {deck_part}, line: {l}.')

				cards.append(
					card_object[_type](_id, _type, definition0=d, comment0=c, header=stored_comments, subtype=data_card_key)
					)
				stored_comments = []

		return cards, stored_comments
	

	def Parse(self, scope=['cells', 'surfaces', 'materials', 'transformations', 'source', 'tallies', 'settings']):
		""" Parses MCNP input deck. Assigns attributes to all deck objects. """

		for _type in scope:
			
			if _type == 'cells':
				cards_dict = self.GetCardGroup(_type)
				
				for c in cards_dict.values():
					_u = c.Parse() # for cells, Parse returns universe
					
					if _u in list(self.universes.keys()):
						self.universes[_u].append(c.id)
					else:
						self.universes[_u] = [c.id]

			else:
				cards_dict = self.GetCardGroup(_type)
				for c in cards_dict.values():
					c.Parse()
			
		self.parsed = scope

	def GetCardGroup(self, _type:str):
		if _type == 'cells':
			return self.cells
		if _type == 'surfaces':
			return self.surfaces
		if _type == 'materials':
			return self.materials
		if _type == 'transformations':
			return self.transformations
		if _type == 'settings':
			return self.settings
		if _type == 'source':
			return self.source
		if _type == 'tallies':
			return self.tallies
	
	def Save(self, savename:str, data_cards_order:list=['transformations', 'materials', 'settings', 'tallies', 'source']):
		
		fh = open(savename, 'w')

		out = self.title + '\n'

		""" for c in self.cells.values():
			out += str(c) + '\n' """
		out += '\n'.join([str(s) for s in self.cells.values()] + self.cell_comment_tail) #  + '\n'
		# out += '\n'.join(self.cell_comment_tail) + '\n'

		out += '\n\n'

		# for s in self.surfaces.values():
		out += '\n'.join([str(s) for s in self.surfaces.values()] + self.surf_comment_tail) # + '\n'
		# out += '\n'.join(self.surf_comment_tail) + '\n'

		out += '\n\n'

		# for ol in self.other_lines:
		# out += '\n'.join(self.other_lines)

		data_out = []

		for dt in data_cards_order:
			data_out += [str(s) for s in self.GetCardGroup(_type=dt).values()]
		data_out += self.data_comment_tail
		out += '\n'.join(data_out)
		# 	data_out.append('\n'.join([str(s) for s in self.GetCardGroup(_type=dt).values()]))
		# out += '\n'.join(_do for _do in data_out) + '\n'
		# out += '\n'.join(self.data_comment_tail) + '\n'
		
		fh.write(out)

		fh.close()

	def Condense(self):
		""" Function goes trough whole MCNP input file and removes all comments and unnecessary whitespace. """

		raise NotImplementedError()

	def ExtractGeometry(self, cell_ids:list, universe_recursion:bool=True):
		""" Extract a part of the geometry.
			Returns a Deck object with all necesarry cells, surfaces, materials, ... """
		
		if 'cells' not in self.parsed:
			#
			#	To extract geometry we need to parse at least all cells. This is needed to obtain whole universes!
			#	Surfaces, materials and transformations can be parsed as needed.
			#
			print('Cells not yet parsed. Parsing all cells  ...')
			self.Parse(scope=['cells'])

		if 'surfaces' not in self.parsed:
					#
					#	Surfaces need to be parsed as they can reference transformation cards
					# 	or other surfaces (periodic surfaces).
					#
					print('Surfaces not yet parsed. Parsing all surfaces  ...')
					self.Parse(scope=['surfaces'])

		_cells = []
		_surfaces = []
		_materials = []
		_transformations = []
		
		_universes = [] # This should be a dictionary not a list!

		subuniverses = True # do our cells reference other universes ('fill')

		_new_cells = cell_ids

		while subuniverses:

			subuniverses = False
			_new_fills = []

			for cid in _new_cells:
				
				_cells.append(cid)
				_surfaces += list(self.cells[cid].surfaces)
				if 0 in self.cells[cid].surfaces:
					print(f'\n0 in surfaces. Found in cell {cid}:')
					print(self.cells[cid].__repr__())
				
				if self.cells[cid].material != 0:
					_materials.append(self.cells[cid].material)

				try:
					_trcl = self.cells[cid].parameters['trcl']
					_transformations.append(int(re.search(r'\d+', _trcl).group()))
				except KeyError: pass

				if 'fill' in self.cells[cid].parameters.keys():
					_fill = self.cells[cid].parameters['fill']
					_new_fills.append(int(re.match(r'\d+', _fill).group()))
					subuniverses = True
				elif '*fill' in self.cells[cid].parameters.keys():
					_fill = self.cells[cid].parameters['*fill']
					_new_fills.append(int(re.match(r'\d+', _fill).group()))
					subuniverses = True
				# except KeyError: pass

			if not universe_recursion: break

			_universes += _new_fills

			_new_cells = []

			for _f in _new_fills: _new_cells += list(self.universes[_f])

			_new_cells = np.unique(_new_cells)
		
		_cells = np.unique(_cells)
		_surfaces = np.unique(_surfaces)

		for _s in _surfaces:
			if self.surfaces[_s].transformation:
				_transformations.append(self.surfaces[_s].transformation)

			#
			#	ADD FUNCTION TO PARSE ALL PERIODIC SURFACES
			#
			

		_materials = np.unique(_materials)
		_transformations = np.unique(_transformations)
		_universes = np.unique(_universes)

		# if 0 in _surfaces:
		#	print('0 in surfaces')

		#
		#	CREATE & RETURN NEW DECK
		#

		new_deck = Deck(filename='', is_blank=True)
		
		new_deck.is_clone = True
		new_deck.title = f' extracted from {self.filename}'
		new_deck.cells = {int(i) : self.cells[int(i)] for i in _cells}
		new_deck.cell_comment_tail = []
		new_deck.surfaces = {int(i) : self.surfaces[int(i)] for i in _surfaces}
		new_deck.surf_comment_tail = []
		new_deck.materials = {int(i) : self.materials[int(i)] for i in _materials}
		new_deck.transformations = {int(i) : self.transformations[int(i)] for i in _transformations}
		new_deck.data_comment_tail = []
		new_deck.settings = {}
		new_deck.tallies = {}
		new_deck.source = {}


		return new_deck #  _cells, _surfaces, _materials, _transformations, _universes

		# loop trough all subuniverses

		# loop trough all surfaces -- surface transformations # we need a surface parser

	def CommentCards(self, _type:str, ids:list):
		""" Changes cards of `_type`, with corresponding `ids` into comments.
			`_type` is cells | surfaces | materials | transformations | tallies | source
			`ids` is list of int
			Not for removing settings. Settings may be removed with:
			`deck.settings.pop(key)`. """
		
		_group = self.GetCardGroup(_type)
		# pop_keys = []
		for k,v in _group.items():
			if v.id in ids:
				v.IntoComment()		
		# raise NotImplementedError()
		
	def Clone(self):
		""" Returns a copy of the Deck. """
		new_deck = Deck(filename='',is_blank=True)
		for attr_name, attr_value in self.__dict__.items():
			setattr(new_deck, copy.deepcopy(attr_name), copy.deepcopy(attr_value)) # this is slow but unfortunately it maybe necessary
		new_deck.is_clone = True
		return new_deck

	def SetSettings(self, settings:dict={}):
		""" Overwrites Settings cards. Usage: `settings={'nps' : '1e9'}` """

		self.settings = {}
		for k,v in settings.items():
			definition = f'{k} {v}'
			self.settings[k.lower()] = Settings(_id='', _type='settings', definition0=definition, subtype=k.lower())
	
	def AddSettings(self, settings:dict={}):
		""" Adds Settings cards. Usage: `settings={'nps' : '1e9'}` """

		try:
			existing_settings = list(self.settings.keys())
		except ValueError:
			existing_settings = []
			self.settings = {}

		for _k,v in settings.items():
			k = _k.lower()
			if k in existing_settings:
				print(f'Warning! Overwriting setting {k}:')
				
			definition = f'{k} {v}'
			self.settings[k.lower()] = Settings(_id='', _type='settings', definition0=definition, subtype=k.lower())

	def AddCell(self, _id:int, region:str, material:int=None, density:float=None, name:str='', parameters:dict={}):
		""" Adds Cell object to Deck.cells dictionary."""

		if _id in list(self.cells.keys()):
			print(f'Warning! Cell with ID={_id} already exists. Overwriting existing cell ...')

		if material and density:
			mat_and_dens = f'{material} {density}'
		else:
			mat_and_dens = '0'
		def0 = f'{_id} {mat_and_dens} {region}'
		c = Cell(_id=_id, _type='cell', definition0=def0, comment0=name)
		c.Parse()
		for k,v in parameters.items():
			c.parameters[k.lower()] = v
		c.Update()
		
		self.cells[_id] = c

		return c

	def AddSurface(self, _id:int, _type:str, parameters:list, transformation:int=None, periodic:int=None, boundary_type:str='vacuum' ,name:str='', precission:int='none', exponential_notation:bool=False):
		""" Adds Surface object to Deck.surfaces dictionary.
			The surfaces `parameters` will be weitten as an ASCII string.
			The user should specify the `precission` and if the `exponential_notation` should be used.
			`parameters` may also be provided as a list of strings, and will not be reformated unless precission is specified.
			If `parameters` are provided as a list of floats and `precission` is `'none'`,
			then defult python f-string formating will be used. """

		_type = _type.lower()

		if _id in list(self.surfaces.keys()):
			print(f'Warning! Surface with ID={_id} already exists. Overwriting existing surface ...')

		if boundary_type == 'vacuum':
			_id_ = f'{_id}'
		elif boundary_type == 'reflective':
			_id_ = f'*{_id}'
		elif boundary_type == 'white':
			_id_ = f'+{_id}'
		else:
			raise ValueError(f'Invalid surface `boundary_type`. Must be either `vacuum`, `reflective` or `white`. Provided value was `{boundary_type}`.')

		if transformation:
			_trans_period_and_type = f'{transformation} {_type}'
		elif periodic:
			_trans_period_and_type = f'-{periodic} {_type}'
		else:
			_trans_period_and_type = str(_type)

		# IF parameters is list of string
		#	IF precission is 'none'	=> keep same string format -- done
		#	IF precission is Int 	=> convert to float and reformat -- exponential?
		# IF parameters is list of float
		#	IF precission is 'none'	=> use default f-strings precission -- exponential?
		#	IF precission is Int 	=> convert to float and reformat -- eponential?

		if type(parameters[0]) == str and precission == 'none':
			params = ' '.join(parameters)
		else:

			_parameters = [float(p) for p in parameters]

			if exponential_notation:
				_dt = 'e' # set data type
			else:
				_dt = 'f'

			if precission == 'none':
				params = ' '.join([f'{p:{_dt}}' for p in _parameters])
			else:
				params = ' '.join([f'{p:{precission}{_dt}}' for p in _parameters])

			""" parameters = [float(p) for ]

		if exponential_notation:
			params = ' '.join([f'{p:.{precission}e}' for p in parameters])
		else:
			params = ' '.join([f'{p:.{precission}f}' for p in parameters]) """

		def0 = f'{_id_} {_trans_period_and_type} {params}'
		s = Surface(_id=_id, _type='surface', definition0=def0, comment0=name)
		s.Parse()
		s.Update()
		
		self.surfaces[_id] = s

		return s
	
	def AddTally(self, _id:int, _type:str, parameters:str, particle:str='n'):
		""" Add tally card.
			`_type` corresonds to the letter before the tally ID.
			Eg.: `f`, `fc`, `e`, `c`, `t`, ...
			`parameters` define the tally card. Eg.: For `e` type these are energy bins,
			while for `fc` this is the string comment of the tally card.
			`particle` is considered only for base tally cards such as `f`.
			
			Often, to fully define a tally the `Deck.AddTally()` method
			has to be called multiple times. """
		
		_type = _type.lower()
		particle = particle.lower()

		if _type in BASE_TALLY_CARDS:
			_id_def = f'{_type}{_id}:{particle}' # this is used in the tally definition
		else:
			_id_def = f'{_type}{_id}'
		
		_id_key = f'{_type}{_id}' # this is used to identify the tally card in the deck
		
		if _id_key in list(self.tallies.keys()):
			print(f'Warning! Tally card "{_id_key}" already exists. Overwriting existing surface ...')

		def0 = f'{_id_def} {parameters}'
		
		### this is diferent compared to AddCell and AddSurface, #############
		### since Tally.Parse() and Update have not yet bee implemented ######

		t = Tally(_id=_id_key, _type='tally', definition0='', subtype=_type)
		t.NewDefinition(def0) 

		######################################################################

		self.tallies[_id_key] = t

		return t
		
	def AddSource(self, _type:str, parameters:str, _id:int=0):
		""" Add source card.
			`_type` is `sdef`, `si`, `sp`, `kcode`, `burn`, ...
			`_id` should be used with source distribution and bias types (eg. `si`)
			`parameters` define the source card. Eg.: they are keyword strings for the `sdef` type (`pos=0 0 1 erg=d1) or strings of values for distribution cards such as `si` (`0.1 1.0 10.0 ...`).
			`_id` should be 0 for `sdef` card.
			
			Often, to fully define a source the `Deck.AddSource()` method
			has to be called multiple times. """
		
		_type = _type.lower()

		if _id == 0 or _type in BASE_SOURCE_CARDS:
			_id_key = _type
		else:
			_id_key = f'{_type}{_id}'
		
		if _id_key in list(self.tallies.keys()):
			print(f'Warning! Tally card "{_id_key}" already exists. Overwriting existing surface ...')

		def0 = f'{_id_key} {parameters}'
		
		### this is diferent compared to AddCell and AddSurface, #############
		### since Source.Parse() and Update have not yet bee implemented #####

		s = Source(_id=_id_key, _type='source', definition0='', subtype=_type)
		s.NewDefinition(def0) 

		######################################################################

		self.source[_id_key] = s

		return s


	def GetAvailableCellIDs(self):
		""" Returns an array of ranges of available cell IDs.
			Array has shape (N,2), where `N` is the number of available spots.
			The range is non-inclusive: e.g. [3,7] means that cell IDs 4,5 and 6
			are available.
			The first range entry (N=0) is [last_cell_ID, first_cell_ID]. """
		cids = np.sort(list(self.cells.keys()))
		diff_idxs = np.argwhere(np.greater(np.diff(cids), 1)).flatten()
		available_ids = np.array([[cids[di], cids[di+1]] for di in diff_idxs])
		return available_ids