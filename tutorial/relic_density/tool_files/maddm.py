#! /usr/bin/env python
import os
import sys
root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
exe_path = os.path.join(root_path,'bin','mg5_aMC')
sys.argv.pop(0)
os.system('%s  -O -W ignore::DeprecationWarning %s %s --mode=maddm' %(sys.executable, str(exe_path) , ' '.join(sys.argv) ))
