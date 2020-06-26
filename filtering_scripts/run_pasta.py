#!/usr/bin/python


import sys, os

d = os.path.abspath(sys.argv[1])

f = 0
for i in os.listdir(d):
    locus = i.split('/')[-1][0:4]

    with open('pasta{}.sh'.format(f), 'w+') as outf:
        outf.write('#!/usr/bin/bash\n')
        outf.write('#SBATCH --workdir=./\n')
        outf.write('#SBATCH --mail-user=leea33@tcnj.edu\n')
        outf.write('#SBATCH --mail-type=ALL\n')
        outf.write('#SBATCH --job-name=pasta\n')
        outf.write('#SBATCH --output=pasta_%j.log\n')
        outf.write('#SBATCH --ntasks=32\n')
        outf.write('#SBATCH --nodes=1\n')
        outf.write('#SBATCH --partition=normal\n')
        outf.write('#SBATCH --time=1-00:00:00\n')

        outf.write('python /home/hpc/leea33/github/pasta/run_pasta.py -i {} -o ./{} --max-mem-mb 4096\n'.format(os.path.join(d,i), locus))
    
    os.system('sbatch pasta{}.sh'.format(f))
    f += 1
    
