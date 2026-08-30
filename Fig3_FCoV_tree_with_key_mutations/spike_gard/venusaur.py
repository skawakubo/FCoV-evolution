from Bio import SeqIO
import sys
import subprocess

alin = sys.argv[1]

skip = []

# this is the breakpoint list for the alignment with 313vir with RB121 and RpYN67
breakpoints = [1383, 1871, 2454, 3276, 3883, 4506, 5262, 5871, 6331, 7519, 7860, 8565, 9414, 10776, 12234, 12498, 13302, 15221, 16001, 16709, 17534, 18809, 20147, 20552, 21461, 21717, 21954, 22178, 23109, 23283, 25272, 25698, 26028, 26281, 26491, 26976, 27350, 27791, 28800, 29240, 29620, 30055, 30386]


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
