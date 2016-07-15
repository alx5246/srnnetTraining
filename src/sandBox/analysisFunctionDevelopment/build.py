import os
import subprocess

if __name__ == '__main__':
    wd = os.path.dirname(os.path.realpath(__file__))

    # fileName = raw_input("Enter file name to run:")
    fileName = "controller.py"
    visualStudioVcVarsAllPath = "C:\\Program Files (x86)\\Microsoft Visual Studio 10.0\\VC\\vcvarsall.bat"


    argsAndApp="(\"%s\" amd64 > nul) && python \"%s\"" % (visualStudioVcVarsAllPath,
                                                  os.path.join(wd, fileName))

    # argsAndApp="(\"%s\" amd64 > nul)" % visualStudioVcVarsAllPath
    print("RUNNING COMMAND: %s" % argsAndApp)
    childProcess = subprocess.Popen(argsAndApp, cwd=wd, shell=True, stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT, env=os.environ, bufsize=1)
    for line in iter(childProcess.stdout.readline, b''):
        print(str(line.rstrip()))
    childProcess.communicate()
    if childProcess.returncode != 0:
        print("\n\nFAILURE: return code %s" % childProcess.returncode)
    else:
        print("\n\nSUCCESS")