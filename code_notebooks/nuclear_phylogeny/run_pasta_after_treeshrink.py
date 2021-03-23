#!/usr/bin/python


import sys, os

d = os.path.abspath(sys.argv[1])

commands = []

for i in os.listdir(d):
    locus = i.split('/')[-1]
    aln = os.path.join(d, i, 'output.fasta')

    commands.append('module add mafft; python /projects/clement-lab/resources/software/pasta/run_pasta.py -i {} -o ./{} --max-mem-mb 4096\n'.format(aln, locus))

with open('commands.txt', 'w+') as outf:
    outf.write('\n'.join(commands))

with open('pasta.sh', 'w+') as outf:
    outf.write('#!/usr/bin/bash\n')
    outf.write('#SBATCH --workdir=./\n')
    outf.write('#SBATCH --mail-user=leea33@tcnj.edu\n')
    outf.write('#SBATCH --mail-type=ALL\n')
    outf.write('#SBATCH --job-name=pasta\n')
    outf.write('#SBATCH --output=pasta_%j.log\n')
    outf.write('#SBATCH --ntasks=60\n')
    outf.write('#SBATCH --partition=normal\n')
    outf.write('#SBATCH --time=1-00:00:00\n')

    outf.write('parallel < commands.txt')

os.system('sbatch pasta.sh')


"""
f = 0
for i in os.listdir(d):
    locus = i.split('/')[-1]

    aln = os.path.join(d, i, 'output.fasta')

    with open('pasta{}.sh'.format(f), 'w+') as outf:
        outf.write('#!/usr/bin/bash\n')
        outf.write('#SBATCH --workdir=./\n')
        outf.write('#SBATCH --mail-user=leea33@tcnj.edu\n')
        outf.write('#SBATCH --mail-type=ALL\n')
        outf.write('#SBATCH --job-name=pasta\n')
        outf.write('#SBATCH --output=pasta_%j.log\n')
        outf.write('#SBATCH --ntasks=4\n')
        outf.write('#SBATCH --nodes=1\n')
        outf.write('#SBATCH --partition=normal\n')
        outf.write('#SBATCH --time=1-00:00:00\n')

        outf.write('module add mafft\n')
        outf.write('python /projects/clement-lab/resources/software/pasta/run_pasta.py -i {} -o ./{} --max-mem-mb 4096\n'.format(aln, locus))
    
    os.system('sbatch pasta{}.sh'.format(f))
    f += 1
"""   
