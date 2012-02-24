import os,sys

write_main_files=False

new_dir='../SixTrack_noastuce/'
def six_cd(fstream,line,basedir):
    global fout
    fout=file(basedir+line.split()[1]+'.if','w')

def six_ca(fstream,line,base):
    global fout
    fout.write('#include "codeblocks/'+line.split()[1]+'.if"\n')

def six_dk(fstream,line,basedir,base):
    global fout, foutmain, deck_dict
    lsp=line.split()[1].lower()
    if lsp in deck_dict:
        if write_main_files:
            foutmain.write('#include "decks/'+lsp+'.F"\n')
        fout=file(basedir+'/decks/'+lsp+'.F','w')
    else: # just need a temporary file to dump the lines for now..
        fout=file('tmp.out','w')

def six_if(fstream,line):
    global fout
    spl=line.split()[1]
    spldsplit=spl.split('.')
    if len(spldsplit)>1 and spldsplit[1].lower()=='not':
            if len(spldsplit)==3:
                fout.write('#ifndef _'+spldsplit[2].upper())
            else:
                fout.write('#if ')
                for el in spldsplit[1:]:
                    if el:
                        if el.upper()=='AND':
                            fout.write(' && ')
                        elif el.upper()=='OR':
                            fout.write(' || ')
                        elif el.upper()=='NOT':
                            fout.write(' ! ')
                        else:
                            fout.write('defined _'+el.upper())
            fout.write('\n')
            return
    if len(spldsplit)==1:
        fout.write('#ifdef _'+spl.upper())
    else:
        fout.write('#if ')
        for el in spldsplit:
            if el:
                if el.upper()=='AND':
                    fout.write(' && ')
                elif el.upper()=='OR':
                    fout.write(' || ')
                elif el.upper()=='NOT':
                    fout.write(' ! ')
                else:
                    fout.write('defined _'+el.upper())
    fout.write('\n')

def six_ei(line):
    global fout
    fout.write("#endif\n")


def get_file_list(base):
    global fin,basedir
    for d in ['decks','codeblocks']:
        if not os.path.isdir(basedir+d):
            os.makedirs(basedir+d)
    cmake=True
    file_dict={}
    prepflags=[]
    f90files=['beamgas']
    # list of flags that are defined by make_six previously..
    prep_flags=['TILT','TRACKING','FAST','CRLIBM','API','DA','COLLIMAT','CPSS',
        'BOINC','CR','NAGFOR','BPM','BEAMGAS','BNLELENS','BIGNBLZ','DEBUG','HDF5']
    # also set by make_six, but opposite logic!
    off_prep_flags=['CERNLIB','NAGLIB']
    prep_flags.append(off_prep_flags)
    for f in os.listdir('ast_mask'):
        fbase=f.split('.')[0]
        for l in file('ast_mask/'+f,'r'):
            lsp=l.split()
            lspd=l.split('.')
            if len(lsp)>0:
                if len(lspd)==2:
                    if lspd[1].strip()=='s':
                        if l.strip()==fin.name:
                            if write_main_files:
                                if fbase in f90files:
                                    file_dict[fbase]=file(basedir+'/'+fbase+'.F90','w')
                                else:
                                    file_dict[fbase]=file(basedir+'/'+fbase+'.F','w')
                                file_dict[fbase].write('#include "'+base+'.h"\n')
                            else:
                                file_dict[fbase]=None
                        else:
                            print "ast mask not used: ",f
                if lsp[0]=='df':
                    if fbase in file_dict:
                        for tmp_pflag in lsp[1:]:
                            for pflag in tmp_pflag.split(','):
                                pflag=pflag.strip().upper()
                                if pflag and pflag not in prepflags:
                                    if pflag in prep_flags: #only add flags that should be set by precompiler..
                                        prepflags.append(pflag.strip())
                                    else:
                                        if write_main_files:
                                            file_dict[fbase].write('#define _'+pflag+'\n')
    if write_main_files:
        inc_header=file(basedir+base+'.h.in','w')
        if cmake: 
            if not os.path.isdir(new_dir+'/cmake'):
                os.makedirs(new_dir+'/cmake')
            cmake_base=file(new_dir+'/cmake/'+base+'_defs.cmake','w')
        for pflag in prepflags:
            if cmake:
                inc_header.write('#cmakedefine _'+pflag+'\n')
                cmake_base.write('option(_'+pflag+' "Turn on '+pflag+'" OFF)\n')
            else:
                inc_header.write('#define '+pflag+'\n')
    return file_dict

def get_deck_dict(fbase):
    global fin,basedir
    cmake=True
    deck_dict={}
    for l in file('ast_mask/'+fbase+'.ast','r'):
        lsp=l.split()
        lspd=l.split('.')
        if len(lsp)>0:
            if lsp[0]=='e':
               for deck in lsp[1].lower().split(','):
                   deck_dict[deck]=file(basedir+'/decks/'+deck+'.F','w')
    return deck_dict

def main(base):
    global fin,fout,foutmain,basedir,deck_dict
    fin=file(base+'.s','r')
    basedir=new_dir+'/'+base+'/'

    if not os.path.isdir(basedir+base):
        os.makedirs(basedir+base)
    file_dict=get_file_list(base)
    for f in file_dict:
        print "Creating file",f+".F[90]"
        fin.seek(0)
        deck_dict=get_deck_dict(f)
        foutmain=file_dict[f]
        for l in fin:
            if l[0]=='+':
                six_cmd=l.split()[0][1:]
                if six_cmd=='ei':
                    six_ei(l.strip())
                elif six_cmd=='if':
                    six_if(fin,l.strip())
                elif six_cmd=='cd':
                    six_cd(fin,l.strip(),basedir+'codeblocks/')
                elif six_cmd=='ca':
                    six_ca(fin,l.strip(),base)
                elif six_cmd=='dk':
                    six_dk(fin,l.strip(),basedir,base)
                else:
                    raise ValueError(six_cmd+" not handled!")
            else:
                fout.write(l)
if __name__=="__main__":
    for base in ['sixtrack','lielib','dabnew']:
        print "Converting:",base+'.s'
        main(base)
