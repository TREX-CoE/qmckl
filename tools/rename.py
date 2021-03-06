#!/usr/bin/env python

"""
Changes the name of a function into all the org files.
This script should be run in the src directory.
"""

import sys
import os


def help():
    print("Syntax : {0}  OLD_FUNC_NAME  NEW_FUNC_NAME".format(sys.argv[0]))



def replace_in_file(filename, old_func_name, new_func_name):
    with open(filename,'r') as f:
        text = f.read()
        
    new_text = text.replace(old_func_name, new_func_name)
                 
    with open(filename,'w') as f:
        f.write(new_text)


def main():
    if len(sys.argv) != 3:
        help()
        sys.exit(-1)
    old_func_name = sys.argv[1]
    new_func_name = sys.argv[2]

    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".org"):
            replace_in_file(filename, old_func_name, new_func_name)

    print("Done. run git diff to check what has been changed.")
    



if __name__ == "__main__":
    main()
