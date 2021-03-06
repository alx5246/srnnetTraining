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

# currentDir = os.path.dirname(os.path.realpath(__file__))
# scriptDir = os.path.join(currentDir, "cythonBuilt")
# sys.path.extend([currentDir, scriptDir])


print(os.path.abspath(".."))
newDir = os.path.join(os.path.abspath("../"),'networkGenerators')
if newDir not in sys.path:
    print('not there')
    sys.path.insert(0,newDir)

#if newDir not in sys.path:
#    print('still not there')