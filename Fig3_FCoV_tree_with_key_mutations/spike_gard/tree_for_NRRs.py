from Bio import SeqIO
import sys
import subprocess

alin = sys.argv[1]

skip = []

# this is the breakpoint list for the alignment with 120vir FCoV spike with one SeCoV as an outgroup
breakpoints = [915, 1425, 2302, 3008, 3371]


for i in range(len(breakpoints)):

    run = i +1
    
    if run not in skip:
            
        alfile = alin.replace('.fas', '_reg%i.fas'%run)
            
        start = 0
        if run != 1:
            start = breakpoints[i-1]
        
        end = breakpoints[i]
        
        fasta = ''    
        
        records = SeqIO.parse(alin, 'fasta')
        
        for rec in records:
            desc = rec.description
            thisseq = str(rec.seq)[start:end]
            fasta = fasta + '>%s\n%s\n'%(desc, thisseq)
        
        with open(alfile, 'w+') as out:
            out.write(fasta[:-1])
        
        cmd = 'iqtree3 -T AUTO -s %s -m GTR+F+R4 -B 1000'%alfile
        subprocess.Popen([cmd], shell = True,close_fds=True).wait()

        #do the last region
        if run == len(breakpoints):
            
            run = run +1
            alfile = alin.replace('.fas', '_reg%i.fas'%run)
                
            start = breakpoints[i]
            
            fasta = ''    
            
            records = SeqIO.parse(alin, 'fasta')
            
            for rec in records:
                desc = rec.description
                thisseq = str(rec.seq)[start:]
                fasta = fasta + '>%s\n%s\n'%(desc, thisseq)
            
            with open(alfile, 'w+') as out:
                out.write(fasta[:-1])
            
            cmd = 'iqtree3 -T AUTO -s %s -m GTR+F+R4 -B 1000'%alfile   
            subprocess.Popen([cmd], shell = True,close_fds=True).wait()        
