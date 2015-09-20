import subprocess
import os
import sys


def install():
    
    allan_path = os.getcwd() + "/" + "allan_5.c"
    print " A path: %s" % allan_path
    
    #Determine platform    
    if sys.platform == "darwin":
        cmd = "gcc -o Run %s" % allan_path
    elif sys.platform == "linux":
        pass #TODO: Add linux support
    else:
        print "OS not supported"

    print "Compiling required C program: allan_5.c"
    print cmd
    subprocess.call(cmd, shell=True)
    

    #Sys-Link to /usr/local/bin
    subprocess.call('cp allan_runner.py allan_runner', shell=True) #Make executable copy
    subprocess.call('chmod +x allan_runner', shell=True) #Make an executable
    print "Sys-linking to /usr/local/bin"
    subprocess.call('ln -s allan_runner /usr/local/bin', shell=True) #Sys link




if __name__ == "__main__":
    install()

