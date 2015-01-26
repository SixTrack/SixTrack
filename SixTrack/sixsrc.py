#!/usr/bin/ipython

"""Sixtrack code analyzer

R. De Maria

Copyright 2014 CERN. This software is distributed under the terms of the GNU
Lesser General Public License version 2.1, copied verbatim in the file
``COPYING''.

In applying this licence, CERN does not waive the privileges and immunities
granted to it by virtue of its status as an Intergovernmental Organization or
submit itself to any jurisdiction.

Source:
  asts: containts templates (.ast) for generating fortran files (.f)
  sources: source files containing the code in decks and common decks
  blocks: common decks that compose the decks or other common decks
  flags: flags used to compose the decks

  Usage:
  s=Source('./')
"""



import os
import re

from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter
from pygments.token import Token

lexer = get_lexer_by_name("fortran", stripall=True)
formatter = HtmlFormatter(linenos=True, cssclass="default",anchorlinenos=True)


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
    self.get_definitions()
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
          msg="Warning: deck not found '%s' used in '%s.ast'"
          print msg%(name,ast.name)
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
          msg="Warning: block not found '%s' used in block '%s'"
          print msg%(name,block.name)
      for name in block.flag_used:
        for astname in block.used_in_ast:
          ast=self.asts[astname]
          if name not in ast.flags:
            msg="Warning: flag not found '%s' used in deck '%s' in '%s.ast'"
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
  tp_source="""<!DOCTYPE html>
    <html>
    <body>
    <h1><a href="index.html">SixTrack Source Browser</a></h1>
    <p>SixTrack use <em>astuce</em> files to generate Fortran files by
assembling code <em>blocks</em> according to <em>flags</em>. </p>
    <p> Astuce files:</p>
    <ul>
    %s
    </ul>
    <p>Flags used: %s.</p>
    <p>Blocks defined: %s.</p>
    </body>
    </html>"""
  tl_ast="""<li><a href="ast_%s.html">%s</a> generates
            <a href="file_%s.html">%s</a></li>"""
  tl_flag='<a href="flag_%s.html">%s</a>'
  tl_block='<a href="block_%s.html">%s</a>'
  def make_html(self,basedir='/afs/cern.ch/project/sixtrack/www/sixtrack_src_html',flags=None):
    asts=sorted([(n.name,n.name,n.outfn[:-3],n.outfn) for n in self.asts])
    flags=sorted([n.name for n in self.flags])
    blocks=sorted([n.name for n in self.blocks])
    ast='\n   '.join([self.tl_ast%n for n in asts])
    flag=', '.join([self.tl_flag%(n,n) for n in flags])
    block=', '.join([self.tl_block%(n,n) for n in blocks])
    page=self.tp_source%(ast,flag,block)
    fh=open(os.path.join(basedir,'index.html'),'w')
    fh.write(page)
    for n in self.asts:
      fh=open(os.path.join(basedir,'ast_%s.html'%n.name),'w')
      fh.write(n.make_html())
      fname=n.outfn.replace('n.f','')
      fh=open(os.path.join(basedir,'file_%s.html'%fname),'w')
      fh.write(n.make_fortran_html(self))
    for n in self.flags:
      fh=open(os.path.join(basedir,'flag_%s.html'%n.name),'w')
      fh.write(n.make_html())
    for n in self.blocks:
      fh=open(os.path.join(basedir,'block_%s.html'%n.name),'w')
      fh.write(n.make_html())
  def write_fortran(self,basedir='build',flags=None):
    if not os.path.exists(basedir):
      os.mkdir(basedir)
    for ast in self.asts:
      fn=os.path.join(basedir,ast.outfn)
      print fn
      out=ast.make_fortran(self,flags=flags)
      open(fn,'w').write(out)
  def get_definitions(self):
      rsub=re.compile('subroutine +([A-z0-9_]+) *\(',re.IGNORECASE)
      rscall=re.compile('call +([A-z0-9_]+) *\(',re.IGNORECASE)
      rfun=re.compile('function +([A-z0-9_]+) *\(',re.IGNORECASE)
      rfcall=re.compile('([A-z0-9_]+)\(',re.IGNORECASE)
      sdef={}
      scall={}
      fdef={}
      fcall={}
      ttt=[(rsub,sdef),(rscall,scall),(rfun,fdef),(rfcall,fcall)]
      for block in self.blocks:
        for sname,src in block.code.items():
          for lineno,line in enumerate(src.splitlines()):
            for reg,out in ttt:
              res=reg.search(line)
              if res:
                name=res.groups()[0].lower()
                data=(block.name,sname, lineno,line.strip())
                out.setdefault(name,[]).append(data)
      #for k,v in fcall.items():
      #  if k not in fdef:
      #    del fcall[k]
      self.defs=Lookup()
      for k in sdef:
         ndef=Definition(k,'subroutine')
         for bname,sname, lineno, ctx in sdef[k]:
            ndef.defined_in.append([bname,sname, lineno, ctx])
            if k not in scall:
              print "Warning: Subrourine `%s' defined in `%s' not called"%(k,bname)
         for calls in scall.get(k,[]):
            bname,sname, lineno, ctx=calls
            ndef.called_in.append([bname,sname, lineno, ctx])
         setattr(self.defs,k,ndef)
      for k in fdef:
         ndef=Definition(k,'function')
         for bname,sname, lineno, ctx in fdef[k]:
            ndef.defined_in.append([bname,sname, lineno, ctx])
            if k not in fcall:
              print "Warning: Function `%s' defined in `%s' not called"%(k,bname)
         for calls in fcall.get(k,[]):
            bname,sname, lineno, ctx=calls
            ndef.called_in.append([bname,sname, lineno, ctx])
         setattr(self.defs,k,ndef)

class Definition(object):
  tp_def="""<!DOCTYPE html>
    <html>
    <body>
    <h1><a href="index.html">SixTrack Source Browser</a>: Definition %s</h1>
    %s
    </body>
    </html>"""
  def __init__(self,name,deftype):
    self.name=name
    self.deftype=deftype
    self.defined_in=[]
    self.called_in=[]



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
  tp_ast="""<!DOCTYPE html>
    <html>
    <body>
    <h1><a href="index.html">SixTrack Source Browser</a>: Astuce File %s</h1>
    <p>Source <tt>%s</tt> used to generate <tt>%s</tt>
    with flags:
    <ul><li>%s,</li></ul>
    <p> and decks: </p>
    <ul><li>%s.</li></ul>
    Code:
    <pre>%s</pre>
    </body>
    </html>"""
  tp_fortran="""<!DOCTYPE html>
    <html>
    <head>
    <link href="style.css" rel="stylesheet" type="text/css">
    </head>
    <body>
    <h1><a href="index.html">SixTrack Source Browser</a>
        Fortran file: %s</h1>
    <p>Generated from %s using flags %s, and blocks %s.</p>
    %s
    </body>
    </html>"""
  tl_flag='<a href="flag_%s.html">%s</a>'
  tl_block='<a href="block_%s.html">%s</a>'
  def __init__(self,name,fname):
    self.name=name
    self.fname=fname
    fh=open(fname)
    self.src=fh.readline().strip()
    if not self.src.endswith('.s'):
      raise ValueError,"Error AstFile: %s not conform"%fname
    self.outfn=fh.readline().strip()
    self.flags=[]
    self.blocks=[]
    for l in fh:
      l=l.strip()
      if l.startswith('df'):
        self.flags.extend(parse_comma_list(l[3:]))
      elif l.startswith('e'):
        self.blocks.extend(parse_comma_list(l[2:]))
    self.code=open(fname).read()
  def __repr__(self):
    return "<Ast: %s -> %s>"%(self.src,self.outfn)
  def print_all(self):
    print "%s -> %s"%(self.src,self.outfn)
    print "flags:",','.join(self.flags)
    print "blocks:",','.join(self.blocks)
  def make_html(self):
    flag=', '.join([self.tl_flag%(n,n) for n in sorted(self.flags)])
    block=', '.join([self.tl_block%(n,n) for n in sorted(self.blocks)])
    code=code2html(self.code)
    page=self.tp_ast%(self.name,self.src,self.outfn,flag,block,code)
    return page
  def make_fortran(self,source,flags=None):
    if flags is None:
      flags=self.flags
    ffile=[]
    for bb in self.blocks:
        if bb in source.blocks:
          out=source.blocks[bb].make_fortran(source,flags,self.src)
          ffile.extend(out)
        else:
          msg="Warning: deck `%s` called by `%s` not found in `%s`"
          print msg%(bb,self.name,self.src)
    return '\n'.join(ffile)
  def make_fortran_html(self,source,flags=None):
    flag=', '.join([self.tl_flag%(n,n) for n in sorted(self.flags)])
    block=', '.join([self.tl_block%(n,n) for n in sorted(self.blocks)])
    code=self.make_fortran(source,flags=flags)
    code = highlight(code, lexer, formatter)
    page=self.tp_fortran%(self.outfn,self.src,flag,block,code)
    return page

class Flag(object):
  tp_flag="""<!DOCTYPE html>
    <html>
    <body>
    <h1><a href="index.html">SixTrack Source Browser</a>: Flag %s</h1>
    <p>Used in ast files: %s. </p>
    <p>Used in block: %s.</p>
    </body>
    </html>"""
  tl_ast='<a href="ast_%s.html">%s</a>'
  tl_block='<a href="block_%s.html">%s</a>'
  def __init__(self,name):
    self.name=name
    self.used_in_block=set()
    self.used_in_ast=set()
  def __repr__(self):
    return "<Flag: %s >"%(self.name)
  def make_html(self):
    in_ast=', '.join([self.tl_ast%(n,n) for n in sorted(self.used_in_ast)])
    in_block=[self.tl_block%(n,n) for n in sorted(self.used_in_block)]
    page=self.tp_flag%(self.name,in_ast,', '.join(in_block))
    return page


class Block(object):
  tp_block="""<!DOCTYPE html>
    <html>
    <body>
    <h1><a href="index.html">SixTrack Source Browser</a>: %s %s</h1>
    <p>Flag used: %s.</p>
    <p>Used in Ast Files: %s.</p>
    <p>Used in Block: %s.</p>
    <p>Code: </p>
    %s
    </body>
    </html>"""
  tl_ast='<a href="ast_%s.html">%s</a>'
  tl_flag='<a href="flag_%s.html">%s</a>'
  tl_block='<a href="block_%s.html">%s</a>'
  tl_code='Defined in <code>"%s"</code>\n<pre>%s</pre>'
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
  def make_fortran(self,source,flags=None,src=None):
    if src is None:
      code=self.code.values()[0]
    else:
      code=self.code[src]
    if flags is None:
      flags=source.asts[list(self.used_in_ast)[0]].flags
    state='print'
    ifstate=0
    ffile=[]
    vl=FlagDict([(name,True) for name in flags])
    for line in code.splitlines():
      #print state,ifstate,line[:40]
      if state=='print':
        if line.startswith('+if'):
          condcode=line.split()[1].replace('.',' ')
          cond=eval(condcode.strip(),{},vl)
          #print line,'->',cond
          if cond:
            state='print'
          else:
            state='skip'
        elif line.startswith('+ca'):
          blockname=line.split()[1].lower()
          try:
            block=source.blocks[blockname]
            ffile.extend(block.make_fortran(source,flags,src))
          except:
            msg="Warning: block `%s` called by `%s` not found"
            print msg%(blockname,self.name)
        elif line.startswith('+ei'):
          pass
        else:
          line=line.rstrip()[:80]
          if line=='':
            line=' '
          ffile.append(line)
      elif state=='skip':
        if line.startswith('+if'):
          ifstate+=1
        elif line.startswith('+ei'):
          if ifstate==0:
            state='print'
          else:
            ifstate-=1
    return ffile
  def __repr__(self):
    return "<Block %s: %s >"%(self.btype,self.name)
  def make_html(self):
    typ=self.btype=='+dk' and "Deck" or "Block"
    flag=', '.join([self.tl_flag%(n,n) for n in self.flag_used])
    ast=', '.join([self.tl_ast%(n,n) for n in self.used_in_ast])
    block=', '.join([self.tl_block%(n,n) for n in self.used_in_block])
    code='\n'.join([self.tl_code%(k,code2html(v)) for k,v in self.code.items()])
    page=self.tp_block%(typ,self.name,flag,ast,block,code)
    return page


def code2html(c):
  c=re.sub(r'(\+if +)([a-z0-9_]+)',r'\1<a href="flag_\2.html">\2</a>',c)
  c=re.sub(r'(\.not\.)([a-z0-9_]+)',r'\1<a href="flag_\2.html">\2</a>',c)
  c=re.sub(r'(\.and\.)([a-z0-9_]+)',r'\1<a href="flag_\2.html">\2</a>',c)
  c=re.sub(r'(\+ca +)([a-z0-9_]+)',r'\1<a href="block_\2.html">\2</a>',c)
  c=re.sub(r'(\+dk +)([a-z0-9_]+)',r'\1<a href="block_\2.html">\2</a>',c)
  return c



class FlagDict(dict):
  def __getitem__(self,k):
    if k in self:
      return dict.__getitem__(self,k)
    else:
      return False


def readsrc(s):
   ast=s.asts.sixve
   srx=ast.make_fortran

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


