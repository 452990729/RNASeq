#!/usr/bin/env python2

import sys
import re
import os
import argparse
import ConfigParser
from glob import glob

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

#### SOFT
PYTHON = config.get('SOFTWARE', 'python')
PERL = config.get('SOFTWARE', 'perl')
Rscript = config.get('SOFTWARE', 'Rscript')
SNAKEMAKE = config.get('SOFTWARE', 'snakemake')
FASTP = config.get('SOFTWARE', 'fastp')
STAR = config.get('SOFTWARE', 'STAR')
STRINGTIE = config.get('SOFTWARE', 'stringtie')
GFFCOMPARE = config.get('SOFTWARE', 'gffcompare')

### SCRIPT
RunEnrich = config.get('SCRIPT', 'enrich')
RunDiff = config.get('SCRIPT', 'RunDiff')

#### DATABASE


class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\s+', line_in)
        self.Sample = list_split[0]
        self.Group = list_split[1]
        self.Name = '_'.join(list_split[:2])
        if ',' in list_split[2]:
            self.paired = True
            list_fq = re.split(',', list_split[2])
            self.fq1 = list_fq[0]
            self.fq2 = list_fq[1]
        else:
            self.paired = False
            self.fq = list_split[2]

class Snake(object):
    def __init__(self, process):
        self.process = process
        self.input = ''
        self.output = ''
        self.params = ''
        self.log = ''
        self.threads = ''
        self.shell = ''

    def UpdateInput(self, line_in):
        self.input = line_in

    def UpdateOutput(self, line_in):
        self.output = line_in

    def UpdateParams(self, line_in):
        self.params = line_in

    def UpdateLog(self, line_in):
        self.log = line_in

    def UpdateThreads(self, line_in):
        self.threads = line_in

    def UpdateShell(self, line_in):
        self.shell = line_in

    def WriteStr(self, fn):
        fn.write('rule '+self.process+':\n')
        fn.write('\tinput:\n\t\t'+self.input+'\n')
        if self.output:
            fn.write('\toutput:\n\t\t'+self.output+'\n')
        if self.params:
            fn.write('\tparams:\n\t\t'+self.params+'\n')
        if self.log:
            fn.write('\tlog:\n\t\t'+self.log+'\n')
        if self.threads:
            fn.write('\tthreads: '+self.threads+'\n')
        if self.shell:
            fn.write('\tshell:\n\t\t'+self.shell+'\n')
        fn.write('\n')


def Argparse():
    parser = argparse.ArgumentParser(description="Micro WGS pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the abs output path', required=True)
    parser.add_argument('-g', help='compare group, A:B,C:D')
    parser.add_argument('-a1', help='the read1 adapter', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    parser.add_argument('-a2', help='the read2 adapter', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-p', help='the type of spe',choices=['hg19', ], default='hg19')
    parser.add_argument('-t', help='the core used', default=10)
    parser.add_argument('-sge', help='use sge',  action='store_true')
    parser.add_argument('-lsf', help='use lsf',  action='store_true')
    parser.add_argument('-r', help='run now',  action='store_true')
    argv=vars(parser.parse_args())
    return argv

def HandleRawdata(argv):
    outpath = argv['o']
    RawData = os.path.join(outpath, '0.RawData')
    list_ob = []
    if os.path.exists(RawData):
        os.system('rm -rf '+RawData)
    os.mkdir(RawData)
    with open(argv['c'], 'r') as f:
        os.chdir(RawData)
        for line in f:
            if not line.startswith('#'):
                ob = ReadList(line.strip())
                list_ob.append(ob)
                if ob.paired:
                    Paired = True
                    if ob.fq1.endswith('.gz'):
                        lb = '.fastq.gz'
                    else:
                        lb = '.fastq'
                    os.system('ln -s {} {}'.format(ob.fq1, ob.Name+'_1'+lb))
                    os.system('ln -s {} {}'.format(ob.fq2, ob.Name+'_2'+lb))
                else:
                    Paired = False
                    if ob.fq.endswith('.gz'):
                        lb = '.fastq.gz'
                    else:
                        lb = '.fastq'
                    os.system('ln -s {} {}'.format(ob.fq, ob.Name+lb))
    os.chdir(outpath)
    return list_ob, Paired, lb

def WriteSnake(argv, list_ob, Paired, lb):
    specie= argv['p']
    Database = config.get('DATABASE', specie)
    GTF = Database+'/'+specie+'.exon.gtf'
    outpath = argv['o']
    snakefile = open(os.path.join(argv['o'], 'snakefile.txt'), 'w')
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                                list_ob])))
    snakefile.write('Groups = "{}".split()\n'.format(' '.join(i.replace(':', 'VS') \
                                                             for i in re.split(',', argv['g']))))
    snakefile.write('adapter1 = "{}"\nadapter2 = "{}"\n'.format(argv['a1'], argv['a2']))

    ###all
    All = Snake('All')
    if argv['g']:
        All.UpdateInput('expand("'+outpath+'/7.Enrich/{group}/Transcript/Transcript.G-enrich.pdf", group=Groups),expand("'+outpath+'/7.Enrich/{group}/Gene/Gene.G-enrich.pdf", group=Groups)')
    else:
        All.UpdateInput('expand("'+outpath+'/5.Quantification/{sample}/{sample}_merged.gtf", sample=Samples)')
    All.WriteStr(snakefile)

    ###QC
    QC = Snake('QC')
    if Paired:
        QC.UpdateInput('A = "'+outpath+'/0.RawData/{sample}_1'+lb+'", B = "'+outpath+'/0.RawData/{sample}_2'+lb+'"')
        QC.UpdateOutput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq"')
        QC.UpdateThreads('2')
        QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input.A} -o {output.A} -I {input.B} -O {output.B} --adapter_sequence {adapter1} --adapter_sequence_r2 {adapter2} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html 1>{log.o} 2>{log.e}"')
    else:
        QC.UpdateInput('"'+outpath+'/0.RawData/{sample}'+lb+'"')
        QC.UpdateOutput('"'+outpath+'/1.QC/{sample}.clean.fastq"')
        QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input} -o {output}  --adapter_sequence {adapter1} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html"')
    QC.WriteStr(snakefile)

    ### Align
    Align = Snake('Align')
    Align.UpdateOutput('"'+outpath+'/2.Align/{sample}/{sample}.Aligned.sortedByCoord.out.bam"')
    Align.UpdateThreads('4')
    Align.UpdateLog('e = "logs/{sample}.align.e", o = "logs/{sample}.align.o"')
    if Paired:
        Align.UpdateInput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq"')
        Align.UpdateShell(r'"'+STAR+r' --readFilesIn {input.A} {input.B} --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --genomeDir '+Database+' --outFileNamePrefix '+outpath+'/2.Align/{wildcards.sample}/{wildcards.sample}. 1>{log.o} 2>{log.e}"')
    else:
        Align.UpdateInput('"'+outpath+'/1.QC/{sample}.clean.fastq"')
        Align.UpdateShell(r'"'+STAR+r' --readFilesIn {input} --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --genomeDir '+Database+' --outFileNamePrefix '+outpath+'/2.Align/{wildcards.sample}/{wildcards.sample}. 1>{log.o} 2>{log.e}"')
    Align.WriteStr(snakefile)

    ### Assemble
    Assemble = Snake('Assemble')
    Assemble.UpdateInput('"'+outpath+'/2.Align/{sample}/{sample}.Aligned.sortedByCoord.out.bam"')
    Assemble.UpdateOutput('"'+outpath+'/3.Assemble/{sample}/{sample}.gtf"')
    Assemble.UpdateThreads('4')
    Assemble.UpdateLog('e = "logs/{sample}.as.e", o = "logs/{sample}.as.o"')
    Assemble.UpdateShell('"'+STRINGTIE+' -G '+GTF+' -l {wildcards.sample} -o {output} -p {threads} {input} 1>{log.o} 2>{log.e}"')
    Assemble.WriteStr(snakefile)

    ### Merge
    Merge = Snake('Merge')
    Merge.UpdateInput('expand("'+outpath+'/3.Assemble/{sample}/{sample}.gtf", sample=Samples)')
    string_gtf = ' '.join([outpath+'/3.Assemble/'+i.Name+'/'+i.Name+'.gtf' for i in list_ob])
    Merge.UpdateOutput('"'+outpath+'/4.Merge/stringtie_merged.gtf"')
    Merge.UpdateThreads('1')
    Merge.UpdateLog('e = "logs/merge.e", o = "logs/merge.o"')
    Merge.UpdateShell('"'+STRINGTIE+' --merge -p {threads} -G '+GTF+' -o {output} '+string_gtf+' 1>{log.o} 2>{log.e}"')
    Merge.WriteStr(snakefile)

    ### Compare
    Compare = Snake('Compare')
    Compare.UpdateInput('"'+outpath+'/4.Merge/stringtie_merged.gtf"')
    Compare.UpdateOutput('"'+outpath+'/4.Merge/merged.stats"')
    Compare.UpdateThreads('1')
    Compare.UpdateLog('e = "logs/compare.e", o = "logs/compare.o"')
    Compare.UpdateShell('"'+GFFCOMPARE+' -r '+GTF+' -G -o '+outpath+'/4.Merge/merged {input} 1>{log.o} 2>{log.e}"')
    Compare.WriteStr(snakefile)

    ### Quantification
    Quantification = Snake('Quantification')
    Quantification.UpdateInput('gtf = "'+outpath+'/4.Merge/stringtie_merged.gtf", bam = "'+outpath+'/2.Align/{sample}/{sample}.Aligned.sortedByCoord.out.bam"')
    Quantification.UpdateOutput('"'+outpath+'/5.Quantification/{sample}/{sample}_merged.gtf"')
    Quantification.UpdateThreads('2')
    Quantification.UpdateLog('e = "logs/{sample}.quant.e", o = "logs/{sample}.quant.o"')
    Quantification.UpdateShell('"'+STRINGTIE+' -e -B -p {threads} -G {input.gtf} -o {output} {input.bam} 1>{log.o} 2>{log.e}"')
    Quantification.WriteStr(snakefile)

    if argv['g']:
        ### Diff
        Diff = Snake('Diff')
        Diff.UpdateInput('expand("'+outpath+'/5.Quantification/{sample}/{sample}_merged.gtf", sample=Samples)')
        Diff.UpdateOutput('A = "'+outpath+'/6.Diff/{group}/results_genes_0.05.csv", B = "'+outpath+'/6.Diff/{group}/transcript_results_0.05.csv"')
        Diff.UpdateThreads('1')
        Diff.UpdateLog('e = "logs/{group}.diff.e", o = "logs/{group}.diff.o"')
        Diff.UpdateShell('"'+PYTHON+' '+RunDiff+' -c '+argv['c']+' -o '+outpath+'/6.Diff/ -i '+outpath+'/5.Quantification/ -g '+argv['g']+' 1>{log.o} 2>{log.e}"')
        Diff.WriteStr(snakefile)

        ###EnrichGene
        EnrichGene = Snake('EnrichGene')
        EnrichGene.UpdateInput('"'+outpath+'/6.Diff/{group}/results_genes_0.05.csv"')
        EnrichGene.UpdateOutput('A = "'+outpath+'/7.Enrich/{group}/Gene/Gene.G-enrich.pdf", B = "'+outpath+'/7.Enrich/{group}/Gene/Gene.K-enrich.pdf"')
        EnrichGene.UpdateThreads('1')
        EnrichGene.UpdateLog('e = "logs/{group}.enrichgene.e", o = "logs/{group}.enrichgene.o"')
        EnrichGene.UpdateShell('"'+Rscript+' '+RunEnrich+' '+outpath+'/7.Enrich/{wildcards.group}/Gene {input} 1>{log.o} 2>{log.e}"')
        EnrichGene.WriteStr(snakefile)

        ###EnrichT
        EnrichT = Snake('EnrichT')
        EnrichT.UpdateInput('"'+outpath+'/6.Diff/{group}/transcript_results_0.05.csv"')
        EnrichT.UpdateOutput('A = "'+outpath+'/7.Enrich/{group}/Transcript/Transcript.G-enrich.pdf", B = "'+outpath+'/7.Enrich/{group}/Transcript/Transcript.K-enrich.pdf"')
        EnrichT.UpdateThreads('1')
        EnrichT.UpdateLog('e = "logs/{group}.enricht.e", o = "logs/{group}.enricht.o"')
        EnrichT.UpdateShell('"'+Rscript+' '+RunEnrich+' '+outpath+'/7.Enrich/{wildcards.group}/Transcript {input} 1>{log.o} 2>{log.e}"')
        EnrichT.WriteStr(snakefile)

def RunShell(argv):
    out = open(os.path.join(argv['o'], 'runsnake.sh'), 'w')
    if argv['lsf']:
        out.write(SNAKEMAKE+' '+' '.join(['--cores', argv['t'], '--cluster', "'bsub -q normal -n {threads} -o %J.o -e %J.e'",\
                                      '--printshellcmds', '--snakefile', 'snakefile.txt']))
    else:
        out.write(SNAKEMAKE+' '+' '.join(['--cores', argv['t'], '--printshellcmds', '--snakefile', 'snakefile.txt']))
    out.close()
    if argv['r']:
        os.system('nohup sh runsnake.sh& > snake.log')

def main():
    argv = Argparse()
    list_ob, Paired, lb = HandleRawdata(argv)
    WriteSnake(argv, list_ob, Paired, lb)
    RunShell(argv)


if __name__ == '__main__':
    main()
