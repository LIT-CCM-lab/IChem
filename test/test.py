#!/usr/bin/python3

import sys
import os
import re
import glob
import json
import subprocess
import pdb
import logging
# import difflib
from pathlib import Path

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


# Create stdout handler for logging to the console (logs all five levels)
stdout_handler = logging.StreamHandler()
stdout_handler.setLevel(logging.DEBUG)
stdout_handler.setFormatter("%(asctime)s | %(levelname)8s | %(message)s")

# Add both handlers to the logger
logger.addHandler(stdout_handler)

BASE_DIRECTORY = Path.cwd()
ICHEM_PROD_PATH = BASE_DIRECTORY / "../distrib/IChem_original_project"
ICHEM_TEST_PATH = BASE_DIRECTORY / "../build/IChem" 

if not os.path.exists(ICHEM_PROD_PATH):
    sys.exit(f"Path error for executable ICHEM_PROD_PATH at {ICHEM_PROD_PATH}")
if not os.path.exists(ICHEM_TEST_PATH):
    sys.exit(f"Path error for executable ICHEM_TEST_PATH at {ICHEM_TEST_PATH}")


ICHEM_PROD = ICHEM_PROD_PATH
ICHEM_TEST = ICHEM_TEST_PATH



SCRIPT_FILEPATH = os.path.realpath(__file__)
TEST_DIRPATH = Path(os.path.dirname(SCRIPT_FILEPATH))
DATASET_DIRPATH = TEST_DIRPATH / "dataset"
PROTLIG_DIRPATH = DATASET_DIRPATH / "protlig"
PROTLIG_ENTRY_DIRPATHS = [Path(x) for x in glob.glob(os.path.join(PROTLIG_DIRPATH, "[a-z0-9][a-z0-9][a-z0-9][a-z0-9]"))]


class Bcolors:
	"""
	https://stackoverflow.com/questions/287871/how-do-i-print-colored-text-to-the-terminal
	"""
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'


class IChemCommand:
	STDOUT_FILENAME = "log.out"
	STDERR_FILENAME = "log.err"

	def __init__(self, ichem_path, input_dirpath, root_dirpath, cmd_config):
		self.ichem_path = ichem_path
		self.cmd_config = cmd_config
		self.out_filenames = cmd_config["out_files"] + [self.STDOUT_FILENAME, self. STDERR_FILENAME]

		self.input_dirpath = input_dirpath
		self.root_dirpath = root_dirpath
		self.out_dirpath = self.root_dirpath / cmd_config["id"]

		self.is_input_checked = False
		self.missing_input_files = set()

		self.is_executed = False

	def format_input_paths(self):
		return {k: self.input_dirpath / v for k, v in self.cmd_config["in_files"].items()}

	@property
	def cmd_string(self):
		return self.cmd_config["cmd"].format(
			ichem=self.ichem_path,
			**self.format_input_paths()
		)

	def check_inputs(self):
		self.is_input_checked = True
		for filepath in self.format_input_paths().values():
			if not os.path.isfile(filepath):
				self.is_input_checked = False
				self.missing_input_files.add(filepath)

	def run(self):
		if not os.path.isdir(self.out_dirpath):
			raise FileNotFoundError(f"Output directory {self.out_dirpath} does not exist!")

		os.chdir(self.out_dirpath)

		# Write command
		with open(os.path.join(self.out_dirpath, "cmd.sh"), "w") as f:
			f.write(self.cmd_string+'\n')

		self.check_inputs()

		# TODO: Check if IChem runs properly (check=True catch core dumped only
		if self.is_input_checked:
			try:
				with open(self.STDOUT_FILENAME, "wb") as out, open(self. STDERR_FILENAME, "wb") as err:
					subprocess.run(self.cmd_string.split(), stdout=out, stderr=err, check=True)

				with open(self.STDERR_FILENAME) as f:
					for line in f:
						if 'TRACE' in line:
							self.is_executed = False
							break
					else:
						self.is_executed = True
			except subprocess.CalledProcessError:
				self.is_executed = False
		else:
			self.is_executed = False


class IChemBatch:
	def __init__(self, ichem, input_dirpath, out_dirpath, commands):
		self.ichem = ichem
		self.out_dirpath = out_dirpath

		self.cmds = [IChemCommand(ichem, input_dirpath, out_dirpath, cmd_config) for cmd_config in commands]

	def run(self):
		for cmd in self.cmds:
			cmd.run()


class TestContainer:
	def __init__(self, desc, func, *args):
		self.desc = desc
		self.func = func
		self.args = args

		self.tests = list(func(*args))


class Test:
	_IGNORED_FILES = set(['cmd.sh',])

	PRINT_COLORS = {
		True: Bcolors.OKGREEN,
		False: Bcolors.FAIL,
	}

	STATUS_NAMES = {
		True: "SUCCEED",
		False: "FAILED",
	}

	@classmethod
	def is_files(cls, filepaths):
		filepaths = set(filepaths) - cls._IGNORED_FILES
		for filepath in filepaths:
			is_file = os.path.isfile(filepath)

			if is_file:
				yield f"{filepath} EXISTS", is_file
			else:
				yield f"{filepath} DOES NOT EXIST", is_file

	@classmethod
	def executed_command(cls, command):
		if command.is_executed:
			yield f"command OK", command.is_executed
		else:
			yield f"command NOK", command.is_executed

	@classmethod
	def same_files(cls, ref_dirpath, comp_dirpath):
		is_same_files = len(set(os.listdir(ref_dirpath)) ^ set(os.listdir(comp_dirpath))) == 0

		if is_same_files:
			yield "SAME files", is_same_files
		else:
			yield "DIFFERENT files", is_same_files

	@classmethod
	def missing_files(cls, dirpath, expected_filenames):
		existing_filenames = set(os.listdir(dirpath))
		for expected_filename in expected_filenames:
			# is_existing = expected_filename in existing_filenames
			is_existing = False
			for existing_filename in existing_filenames:
				if re.match(expected_filename, existing_filename) is not None:
					is_existing = True
					break

			if is_existing:
				yield f"{expected_filename} EXISTS", is_existing
			else:
				yield f"{expected_filename} DOES NOT EXIST", is_existing

	@classmethod
	def extra_files(cls, dirpath, expected_filenames):
		existing_filenames = set(os.listdir(dirpath)) - cls._IGNORED_FILES
		for existing_filename in existing_filenames:
			# is_expected = existing_filename in expected_filenames
			is_expected = False
			for expected_filename in expected_filenames:
				if re.match(expected_filename, existing_filename) is not None:
					is_expected = True
					break

			if is_expected:
				yield f"{existing_filename} IS an expected file", is_expected
			else:
				yield f"{existing_filename} IS NOT an expected file", is_expected

	@classmethod
	def inter_files(cls, ref_dirpath, comp_dirpath):
		ref_filenames = set(os.listdir(ref_dirpath)) - cls._IGNORED_FILES
		comp_filenames = set(os.listdir(comp_dirpath)) - cls._IGNORED_FILES

		for ref_filename in ref_filenames:
			is_comp_in_ref = ref_filename in comp_filenames

			if is_comp_in_ref:
				yield f"{ref_filename} ref file IN comp dir", is_comp_in_ref
			else:
				yield f"{ref_filename} ref file NOT IN comp dir", is_comp_in_ref

	@classmethod
	def diff_files(cls, ref_dirpath, comp_dirpath):
		# TODO: add method to exclude some patterns
		inter_filenames = set(os.listdir(ref_dirpath)) & set(os.listdir(comp_dirpath)) - cls._IGNORED_FILES

		for inter_filename in inter_filenames:
			ref_filepath = ref_dirpath / inter_filename
			comp_filepath = comp_dirpath / inter_filename

			is_same_file = True

			# Check file size first
			if os.path.getsize(ref_filepath) != os.path.getsize(comp_filepath):
				is_same_file = False
			else:
				# Compare line by line
				with open(ref_filepath) as f_ref, open(comp_filepath) as f_comp:
					for ref_line, comp_line in zip(f_ref, f_comp):

						if ref_filepath == '/projects/shared/users/cjacquemard/Programmation/IChem/test/grim/ref/1eou/02/log.err':
							pdb.set_trace()

						if 'Modification time' in ref_line:
							continue

						if ref_line != comp_line:
							is_same_file = False
							break

			if is_same_file:
				yield f"{inter_filename} files are EQUAL", is_same_file
			else:
				yield f"{inter_filename} files are DIFFERENT", is_same_file

	@classmethod
	def format_message(cls, message, test):
		return "        +-> " + "{:<9}".format(cls.STATUS_NAMES[test]) + f"{message}"

	@classmethod
	def print_message(cls, message, test):
		print("        +-> " + cls.PRINT_COLORS[test] + "{:<9}".format(cls.STATUS_NAMES[test]) + Bcolors.ENDC + f"{message}")


class RunTest:
	def __init__(self):
		self.__ichem_module = None
		self.__ids = []
		self.__wk_dirpath = None
		self.__conf = None
		self.__is_prepared = False
		self.__is_executed = False
		self.__is_tested = False

		self.prod_batches = {}
		self.test_batches = {}

		self.__test_containers = {}

	@property
	def status(self):
		return self.__status

	@property
	def ref_dirpath(self):
		return self.__working_dirpath / "ref"

	@property
	def comp_dirpath(self):
		return self.__working_dirpath / "comp"

	def add_test(self, entry_name, cmd_id, desc, func, *args):
		if entry_name not in self.__test_containers:
			self.__test_containers[entry_name] = {}

		if cmd_id not in self.__test_containers[entry_name]:
			self.__test_containers[entry_name][cmd_id] = []

		test_container = TestContainer(desc, func, *args)
		self.__test_containers[entry_name][cmd_id].append(test_container)

	def is_succeed(self):
		tests = []
		for entry_name in self.__test_containers:
			for cmd_tests in self.__test_containers[entry_name].values():
				for test_container in cmd_tests:
					for message, test in test_container.tests:
						tests.append(test)

		return all(tests)

	def load_conf(self, conf_filepath):
		with open(conf_filepath) as f:
			conf = json.load(f)

		if "module" not in conf:
			raise KeyError("'module' is missing!")

		if "subdirs" not in conf:
			raise KeyError("'subdirs' is missing!")

		if not conf["subdirs"]:
			raise ValueError("No subdir has been defined")

		cmd_ids = set()

		# Check if one command is defined at least
		if not conf["commands"]:
			raise ValueError("Commands list is empty!")

		for i, ichem_conf in enumerate(conf["commands"]):
			if "id" not in ichem_conf:
				raise KeyError("'id' in entry f{i} is missing!")
			elif "cmd" not in ichem_conf:
				raise KeyError("'cmd' in entry f{i} is missing!")
			elif"in_files" not in ichem_conf:
				raise KeyError("'in_files' in entry f{i} is missing!")
			elif"out_files" not in ichem_conf:
				raise KeyError("'out_files' in entry f{i} is missing!")

			# Check if ids are unique
			if ichem_conf["id"] in cmd_ids:
				raise ValueError(f"Id {ichem_conf['id']} duplicate!")
			else:
				cmd_ids.add(ichem_conf["id"])

			# Check if input files list is empty
			if not ichem_conf["in_files"]:
				raise ValueError(f"No input file for command {ichem_conf['id']}!")

			# Check if output files list is empty
			if not ichem_conf["out_files"]:
				raise ValueError(f"No out files for command {ichem_conf['id']}!")

		self.__conf = conf
		self.__entry_names = self.__conf['subdirs']
		self.__working_dirpath = Path(os.path.dirname(os.path.realpath(conf_filepath)))
		self.__input_dirpath = Path(os.path.realpath(self.__working_dirpath / self.__conf['input_dir']))
		self.__status = 'INIT'

	def prepare(self):
		if self.__conf is None:
			raise ValueError("No configuration has been loaded!")

		if not os.path.isdir(self.ref_dirpath):
			os.mkdir(self.ref_dirpath)

		if not os.path.isdir(self.comp_dirpath):
			os.mkdir(self.comp_dirpath)

		for entry_name in self.__entry_names:
			if not os.path.isdir(self.ref_dirpath / entry_name):
				os.mkdir(self.ref_dirpath / entry_name)

			if not os.path.isdir(self.comp_dirpath / entry_name):
				os.mkdir(self.comp_dirpath / entry_name)

			for cmd in self.__conf["commands"]:
				if not os.path.isdir(self.ref_dirpath / entry_name / cmd["id"]):
					os.mkdir(self.ref_dirpath / entry_name / cmd["id"])
				else:
					for filename in cmd["out_files"]:
						if os.path.isfile(self.ref_dirpath / entry_name / cmd["id"] / filename):
							os.remove(self.ref_dirpath / entry_name / cmd["id"] / filename)

				if not os.path.isdir(self.comp_dirpath / entry_name / cmd["id"]):
					os.mkdir(self.comp_dirpath / entry_name / cmd["id"])
				else:
					for filename in cmd["out_files"]:
						if os.path.isfile(self.comp_dirpath / entry_name / cmd["id"] / filename):
							os.remove(self.comp_dirpath / entry_name / cmd["id"] / filename)

		self.__is_prepared = True

	def run(self):
		if not self.__is_prepared:
			raise ValueError("Test not prepared!")

		print(f"| Testing {self.__conf['module']}...")

		for entry_name in self.__conf['subdirs']:
			entry_dirpath = self.__input_dirpath / entry_name
			out_ref_entry_dirpath = self.ref_dirpath / entry_name
			out_comp_entry_dirpath = self.comp_dirpath / entry_name

			prod_batch = IChemBatch(ICHEM_PROD, entry_dirpath, out_ref_entry_dirpath, self.__conf["commands"])
			self.prod_batches[entry_name] = prod_batch
			prod_batch.run()

			test_batch = IChemBatch(ICHEM_TEST, entry_dirpath, out_comp_entry_dirpath, self.__conf["commands"])
			self.test_batches[entry_name] = test_batch
			test_batch.run()

		self.__is_executed = True

	def compare(self):
		if not self.__is_executed:
			raise ValueError("Command not executed!")

		for entry_name in self.__conf['subdirs']:
			entry_dirpath = self.__input_dirpath / entry_name

			# print(f"+-->| {entry_name}")
			prod_batch = self.prod_batches[entry_name]
			test_batch = self.test_batches[entry_name]

			for i, command in enumerate(self.__conf["commands"]):
				# print(f"    +-->| {command['id']}")
				# print("        +->| Check input files")
				prod_cmd = prod_batch.cmds[i]
				test_cmd = test_batch.cmds[i]

				ref_dirpath = self.__working_dirpath / "ref" / entry_name / command['id']
				comp_dirpath = self.__working_dirpath / "comp" / entry_name / command['id']

				self.add_test(entry_name, command['id'], "Input files", Test.is_files, prod_cmd.format_input_paths().values())
				self.add_test(entry_name, command['id'], "Production command", Test.executed_command, prod_cmd)
				self.add_test(entry_name, command['id'], "Test command", Test.executed_command, test_cmd)
				self.add_test(entry_name, command['id'], "Same number of files", Test.same_files, ref_dirpath, comp_dirpath)
				self.add_test(entry_name, command['id'], "Missing expected files", Test.missing_files, comp_dirpath, command["out_files"])
				self.add_test(entry_name, command['id'], "Extra files", Test.extra_files, comp_dirpath, command["out_files"])
				self.add_test(entry_name, command['id'], "Ref files in comp dir", Test.inter_files, ref_dirpath, comp_dirpath)
				self.add_test(entry_name, command['id'], "File difference", Test.diff_files, ref_dirpath, comp_dirpath)

		if self.is_succeed():
			print(Bcolors.OKGREEN + "SUCCEED" + Bcolors.ENDC)
		else:
			print(Bcolors.FAIL + "FAILED" + Bcolors.ENDC)

		self.__is_tested = True

	def write_results(self, filepath, only_failed=False):
		if not self.__is_tested:
			raise ValueError('No test!')

		with open(filepath, 'w') as f:
			f.write(self.__conf['module']+'\n')

			for entry_name in self.__test_containers:
				f.write(f"+-->| {entry_name}\n")

				for cmd_id, cmd_tests in self.__test_containers[entry_name].items():
					f.write(f"    +-->| {cmd_id}\n")
					for test_container in cmd_tests:
						f.write(f"        +->| {test_container.desc}\n")
						for message, test in test_container.tests:
							if not test and only_failed or not only_failed:
								f.write(Test.format_message(message, test)+'\n')



if __name__ == "__main__":
	tested_modules = ('bsa', 'frag', 'grim', 'ifp', 'ints', 'realign', 'pdbconv', 'sims', 'volsite')
	# tested_modules = ('grim',)

	for tested_module in tested_modules:
		conf_file = os.path.join(TEST_DIRPATH, tested_module, 'conf.json')
		T = RunTest()
		T.load_conf(conf_file)
		T.prepare()
		T.run()
		T.compare()
		T.write_results(os.path.join(TEST_DIRPATH, tested_module, 'test_results.txt'), only_failed=True)
