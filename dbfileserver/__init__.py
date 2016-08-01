
# SYSTEM IMPORTS
import os
import sys

# these variables are used to modify the "PATH" internal to python.
# "sys.path" is a list of directories (strings), and when you try
# to import a python module (call it <module>), python will
# enumerate sys.path and search for <element_in_sys.path>/<module>/__init__.py
# OR <element_in_sys.path>/<module>.py. The first one it comes across, python
# will attempt to import, raising errors if there are any.

# we need access to all subdirectories that have python scripts we wish
# to export. So I am getting the "absolute" path to these directories
# and adding them to sys.path. This allows the db module to be imported
# from any python code.
currentDir = os.path.dirname(os.path.realpath(__file__))
scriptDir = os.path.join(currentDir, "scripts")
sys.path.extend([currentDir, scriptDir])

# PYTHON PROJECT IMPORTS
from DBManager import DBManager
import stdUnpickler
from HTTPRequest import HTTPRequest, urljoin
from DataManager import DataManager

# delete variables so that outside users do not have access to them
del currentDir
del scriptDir