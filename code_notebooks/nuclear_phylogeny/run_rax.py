#!/usr/bin/python

# run many instances of raxml at once

import sys, os

d = os.path.abspath(sys.argv[1])

f = 0
for i in os.listdir(d):
    locus = i.split('/')[-1][0:4]

    with open('raxml{}.sh'.format(f), 'w+') as outf:
        outf.write('#!/usr/bin/bash\n')
        outf.write('#SBATCH --workdir=./\n')
        outf.write('#SBATCH --mail-user=leea33@tcnj.edu\n')
        outf.write('#SBATCH --mail-type=ALL\n')
        outf.write('#SBATCH --job-name=raxml\n')
        outf.write('#SBATCH --output=raxml_%j.log\n')
        outf.write('#SBATCH --ntasks=1\n')
        outf.write('#SBATCH --nodes=1\n')
        outf.write('#SBATCH --partition=normal\n')
        outf.write('#SBATCH --time=1-00:00:00\n')

        outf.write('module add raxml\n')
        outf.write('raxmlHPC -fa -x 65135 -# 500 -p 65135 -m GTRGAMMA -s {} -n {}.rax{}\n'.format(os.path.join(d,i), locus, f))
    
    os.system('sbatch raxml{}.sh'.format(f))
    f += 1
    
