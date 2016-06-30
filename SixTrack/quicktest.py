#!/usr/bin/env python

import os, sys, shutil, time
from glob import glob

def get_all_tests():
    return  sorted(glob("quicktest/*/ref_*"))

def get_timestamp():
    return time.strftime("%Y%m%dT%H%M%S%Z")

def get_hash():
    return os.popen('git log -1 --format="%h"').read().strip()

class Ref(object):
    def __init__(self,ref_dir):
      print "Parsing %s"%(ref_dir)
      if not os.path.isdir(ref_dir):
          raise ValueError, "error: %s"%ref_dir
      self.ref_dir=ref_dir
      head,options=os.path.split(ref_dir)
      self.options=tuple(options.split('_')[1:])
      self.head=head
      self.input_dir=os.path.join(head,"input")
      if not os.path.isdir(self.input_dir):
          raise ValueError, "error: no %s"%self.input_dir
      self.run_script=os.path.join(self.input_dir,"run_test")
      if not os.path.isfile(self.run_script):
          self.mk_default_script()
    def run_prepare(self,exe_path):
      run_dir="run_%s_%s_%s"%("_".join(options),get_hash(),get_timestamp())
      run_dir=os.path.join(self.head,run_dir)
      print run_dir
      os.mkdir(run_dir)
      os.system("cp %s %s"%(exe_path,run_dir))
      os.system("cp %s/* %s"%(self.input_dir,run_dir))
      os.system("cp -a %s %s/ref"%(self.ref_dir,run_dir))
      return run_dir
    def mk_default_script(self):
        tmp="""#!/bin/bash
        ./sixtrack >fort.6
        for iii in ref/*;
        do
          diff $iii ${iii#ref/} >>diff
        done
        """
        open(self.run_script,'w').write(tmp)
        os.chmod(self.run_script,0700)

def expand_tests(tests):
    compiles={}
    for ref_dir in tests:
        ref=Ref(ref_dir)
        compiles.setdefault(ref.options,[]).append(ref)
    return compiles

def compile_exe(options):
    cmd="make OPTIONS='%s'"%(" ".join(options))
    ret=os.popen(cmd).readlines()
    for rl in ret:
        if rl.startswith("Sixtrack build dir:"):
            exe=rl.split(':')[1].strip()+'/sixtrack'
    if not os.path.isfile(exe):
        raise ValueError, "error: no %s"%exe
    print exe
    return exe

def run_test(run_dir):
    cmd="(cd %s; ./run_test >run_test.out 2>run_test.err)"%run_dir
    result=os.system(cmd)
    return result

if __name__=='__main__':
    alltest=sys.argv[1:]
    if len(alltest)==0:
        tests=get_all_tests()
    compiles=expand_tests(tests)
    for options,refs in compiles.items():
        exe_path=compile_exe(options)
        for ref in refs:
           run_dir=ref.run_prepare(exe_path)
           result=run_test(run_dir)
           print result


