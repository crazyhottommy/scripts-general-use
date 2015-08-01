### This is a markdown file edited by [MacDown](http://macdown.uranusjr.com/). I came across the bioinformatics one-liners on the [biostar](https://www.biostars.org/p/142545/) forum and gathered them here.

05/21/2015.


```bash
##1 get the sequences length distribution form a fastq file using awk
zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  

##2 Reverse complement a sequence (I use that a lot when I need to design primers)
echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'  

##3 split a multifasta file into single ones with csplit:
csplit -z -q -n 4 -f sequence_ sequences.fasta /\>/ {*}  

## linearize multiline fasta
cat file.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}'

awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' file.fa

## fastq2fasta
zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[2]\n";' | gzip -c > file.fasta.gz

## bam2bed
samtools view file.bam | perl -F'\t' -ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > file.bed

##bam2wig
samtools mpileup -BQ0 file.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > file.wig.gz

## Number of reads in a fastq file
cat file.fq | echo $((`wc -l`/4))

## Single line fasta file to multi-line fasta of 60 characteres each line
awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/60;i++) {for (j=1;j<=60;j++) printf "%s", $(i*60 +j); print ""}}' file

fold -w 60 file

## Sequence length of every entry in a multifasta file
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' file.fa

## Reproducible subsampling of a FASTQ file. srand() is the seed for the random number generator - keeps the subsampling the same when the script is run multiple times.  0.01 is the % of reads to output.
cat file.fq | paste - - - - | awk 'BEGIN{srand(1234)}{if(rand() < 0.01) print $0}' | tr '\t' '\n' > out.fq

## or look at the Hengli's Seqtk 

## Deinterleaving a FASTQ:
cat file.fq | paste - - - - - - - - | tee >(cut -f1-4 | tr '\t'  
'\n' > out1.fq) | cut -f5-8 | tr '\t' '\n' > out2.fq

## Using mpileup for a whole genome can take forever. So, handling each chromosome separately and parallely running them on several cores will speed up your pipeline. Using xargs you can easily realize it.  
## Example usage of xargs (-P is the number of parallel processes started - don't use more than the number of cores you have available):
samtools view -H yourFile.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d 100000 -uf yourGenome.fa -r {} yourFile.bam | bcftools view -vcg - > tmp.{}.vcf"

## To merge the results afterwards, you might want to do something like this:

samtools view -H yourFile.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | perl -ane 'system("cat tmp.$F[0].bcf >> yourFile.vcf");'

##split large file by id/label/column
awk '{print >> $1; close($1)}' input_file

## sort vcf file with header
cat my.vcf | awk '$0~"^#" { print $0; next } { print $0 | "sort -k1,1V -k2,2n" }'

## Rename a file, bash string manipulation
for file in *gz
do zcat $file > ${file/bed.gz/bed}

## gnu sed print invisible characters
cat my_file | sed -nl

## or cat -v

## exit a dead ssh session
~.

## copy large files, copy the from_dir directory inside the to_dir directory
rsync -av from_dir  to_dir

## copy every file inside the frm_dir to to_dir
rsync -av from_dir/ to_dir

## make directory using the current date
mkdir $(date +%F)

## all the folders' size in the current folder (GNU du)
du -h --max-depth=1

# this one is a bit different, try it and see the difference
du -ch

## the total size of current directory
du -sh .

## disk usage
df -h

## the column names of the file, install csvkit https://csvkit.readthedocs.org/en/0.9.1/
csvcut -n

## open top with human readable size in Mb, Gb. install htop for better visualization
top -M

## how many memeory are used in Gb
free -mg

## print out unique rows based on the first and second column
awk '!a[$1,$2]++' input_file

## do not wrap the lines using less
less -S

## pretty output
fold -w 60
column -t



```


