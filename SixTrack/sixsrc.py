#!/usr/bin/ipython


""" Sixtrack code analyzer

Source:
  asts: containts templates (.ast) for generating fortran files (.f)
  sources: source files containing the code in decks and common decks
  blocks: common decks that compose the decks or other common decks
  flags: flags used to compose the decks

  Usage:
  s=Source('./')
"""



#from objdebug import ObjDebug
import os
import re

#object=ObjDebug

def parse_comma_list(s):
  lst=s.split(',')
  lst=[f.strip() for f in lst]
  lst=[f.lower() for f in lst if f!='']
  return lst

def parse_astuce(filenames):
  blocks={}
  stats={}
  for filename in filenames:
    fh=open(filename)
    code=[]
    def store():
      if len(code)>0:
        block.code[os.path.split(filename)[1]]=''.join(code)
        blocks[block.name]=block
    for iline,line in enumerate(fh):
      if line.startswith('+cd') or line.startswith('+dk'):
        store()
        code=[]
        btype,name=line.split()
        name=name.lower()
        block=blocks.setdefault(name,Block(name,btype))
      else:
        code.append(line)
      if line.startswith('+ca'):
        blockname=line.split()[1]
        block.block_used.add(blockname)
      elif line.startswith('+if'):
        statement=line.split()[1]
        flags=[]
        for token  in statement.split('.'):
          if token not in ['and','or','not','']:
            block.flag_used.add(token)
      stats.setdefault(line,[]).append((iline,name))
    store()
  return blocks,stats

class Source(object):
  def get_ast_names(self):
    astdn=os.path.join(self._basedir,'ast_mask')
    lst=[f for f in os.listdir(astdn) if f.endswith('.ast')]
    return lst
  def parse_asts(self):
    self.asts=Lookup()
    for k in self.get_ast_names():
      name=os.path.splitext(k)[0]
      fname=os.path.join(self._basedir,'ast_mask',k)
      try:
        setattr(self.asts,name,Ast(name,fname))
      except ValueError,msg:
        print msg
  def __init__(self,basedir):
    self._basedir=basedir
    self.process()
  def process(self):
    self.parse_asts()
    self.parse_sources()
    self.fill_deps()
    self.find_unused_blocks()
    self.find_multiple_defined_blocks()
  def get_fortran_files(self):
    return [ast.outfn for ast in self.asts]
  def get_source_names(self):
    return sorted(set([ast.src for ast in self.asts]))
  def parse_sources(self):
    self.blocks=Lookup()
    fnames=[os.path.join(self._basedir,fn) for fn in self.get_source_names()]
    blocks,stats=parse_astuce(fnames)
    self.blocks.__dict__.update(blocks)
    self.stats=stats
  def fill_deps(self):
    self.flags=Lookup()
    for ast in self.asts:
      for name in ast.blocks:
        if name in self.blocks:
          self.blocks[name].used_in_ast.add(ast.name)
        else:
          print "Warning: deck '%s' used in '%s.ast' not found"%(name,ast.name)
      for name in ast.flags:
        if name not in self.flags:
          flag=Flag(name)
          self.flags[name]=flag
        else:
          flag=self.flags[name]
        flag.used_in_ast.add(ast.name)
    for block in self.blocks:
      for name in block.block_used:
        if name in self.blocks:
          self.blocks[name].used_in_block.add(block.name)
        else:
          print "Warning: deck '%s' used in block '%s' not found"%(name,block.name)
      for name in block.flag_used:
        for astname in block.used_in_ast:
          ast=self.asts[astname]
          if name not in ast.flags:
            msg="Warning: flag '%s' used in deck '%s' not found '%s.ast'"
            print msg%(name,block.name,ast.name)
            flag=Flag(name)
          else:
            flag=self.flags[name]
          flag.used_in_block.add(block.name)
  def find_unused_blocks(self):
    out=[]
    for b  in self.blocks:
      src=', '.join(b.code.keys())
      if b.btype=='+dk' and len(b.used_in_ast)==0:
        msg="Warning: block %s %s defined in %s not used in any ast mask"
        print msg%(b.btype,b.name,src)
      if b.btype=='+cd' and len(b.used_in_block)==0:
        msg="Warning: block %s %s defined in %s not used in any other block"
        print msg%(b.btype,b.name,src)
  def find_multiple_defined_blocks(self):
    for  b  in self.blocks:
      if len(b.code)>1:
        src=', '.join(b.code.keys())
        print "Warning: block %s defined in %s"%(b.name,src)
        cd=b.code.values()
        a=cd.pop()
        for b in cd:
          if a!=b:
            print "Warning: ...and source code differs"
  def make_fortran(self,deckname,src=None,astname=None,flags=None):
    block=self.blocks[deckname]
    if astname is not None:
      ast=self.asts[ast]
      flags=ast.flags
    if flags is None:
      flags=[]
    if src is None and astname is not None:
      src=ast.src
    block.make_fortran(flags,src)
  def search(self,regex):
    reg=re.compile(regex)
    for d in self.blocks:
      for src,code in d.code.items():
        for iline,line in enumerate(code.splitlines()):
          if reg.search(line):
            print "%-10s %-10s %s"%(src,d.name,line)
  def get_blocks_usage(self):
    out=[]
    for  b in self.blocks:
      used=len(b.used_in_block)+len(b.used_in_ast)
      out.append((used,b.name))
    return sorted(out)
  def get_lines_usage(self,n=40):
    out=[]
    for line,pos in self.stats.items():
      if len(pos)>n:
        out.append((len(pos),line))
    return sorted(out)


class Lookup(object):
  def __contains__(self,k):
    return k in self.__dict__
  def __iter__(self):
    return self.__dict__.itervalues()
  def __repr__(self):
    return "\n".join(map(str,self))
  def __setitem__(self,k,v):
    setattr(self,k,v)
  def __getitem__(self,k):
    return getattr(self,k)
  def _print_all(self):
    for el in self._get_all():
      el._print_all()


class Ast(object):
  def __init__(self,name,fname):
    self.name=name
    self.fname=fname
    fh=open(fname)
    self.src=fh.readline().strip()
    if not self.src.endswith('.s'):
      raise ValueError,"Warning AstFile: %s not conform"%fname
    self.outfn=fh.readline().strip()
    self.flags=[]
    self.blocks=[]
    for l in fh:
      l=l.strip()
      if l.startswith('df'):
        self.flags.extend(parse_comma_list(l[3:]))
      elif l.startswith('e'):
        self.blocks.extend(parse_comma_list(l[2:]))
  def __repr__(self):
    return "<Ast: %s -> %s>"%(self.src,self.outfn)
  def print_all(self):
    print "%s -> %s"%(self.src,self.outfn)
    print "flags:",','.join(self.flags)
    print "blocks:",','.join(self.blocks)

class Flag(object):
  def __init__(self,name):
    self.name=name
    self.used_in_block=set()
    self.used_in_ast=set()
  def __repr__(self):
    return "<Flag: %s >"%(self.name)


class Block(object):
  def get_source(self,regexp=None,context=1):
    for filename,src in self.code.items():
      print ">>>>>>>>>%s<<<<<<<<<<<"%filename
      if regexp is not None:
        reg=re.compile(regexp)
        code=src.splitlines()
        mmm=[]
        for iline,line in enumerate(code):
          if reg.search(line):
            mmm.extend(range(iline-context,iline+context+1))
        for immm in mmm:
          print code[immm]
      else:
        print src
  def __init__(self,name,btype):
    self.name=name
    self.btype=btype
    self.code={}
    self.used_in_block=set()
    self.block_used=set()
    self.used_in_ast=set()
    self.flag_used=set()
  def make_fortran(self,flags,src=None):
    if src is None:
      code=self.code.values()[0]
    else:
      code=self.code[src]
    state='print'
    ifstate=0
    ffile=[]
    vl=FlagDict([(name,True) for name in flags])
    for line in code.splitlines():
      if line.startswith('+if'):
        condcode=line.split()[1].replace('.',' ')
        cond=eval(condcode.strip(),{},vl)
        #print line,'->',cond
        ifstate+=1
        if cond:
          state='print'
        else:
          state='skip'
      if line.startswith('+ei'):
          ifstate-=1
          #print line
      if state=='print':
        ffile.append(line)
      elif state=='skip':
        if line.startswith('+ei') and ifstate==0:
          state='print'
    if ifstate!=0:
      print "Error in if for %s"%(self.name), ifstate
      raise ValueError
    print self.name,len(ffile)
  def __repr__(self):
    return "<Block %s: %s >"%(self.btype,self.name)

class FlagDict(dict):
  def __getitem__(self,k):
    if k in self:
      return dict.__getitem__(self,k)
    else:
      return False

if __name__=='__main__':
  import sixsrc
  s=sixsrc.Source('./')
  print """

  Instructions:
  The main object is called s. Use s.<TAB> to explore variables and methods. For a example type:
print s.blocks.maincr.used_in_ast
print s.blocks.maincr.block_used
print s.blocks.crlibco.used_in_block
s.blocks.maincr.get_source(regexp='seed',context=5)
"""


