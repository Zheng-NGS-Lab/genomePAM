#!/usr/bin/bash
if [ "$#" -ne 3 ]; then
  if [ "$#" -ne 1 ]; then
    echo "Usage: bash target_annotatation.sh [sample name] [genome name (optional)] [annovar path (optional)]"
    exit
  else
    sample="$1"
    genome="hg19"
    snpeff="/scratch2/pipeline_tools/snpEff"
  fi
else
  sample="$1"
  genome="$2"
  snpeff="$3"
fi
if [ -f identified/${sample}_identifiedOfftargets.txt ]; then
	sed 's/chrchr/chr/g' identified/${sample}_identifiedOfftargets.txt > _t
       	mv _t identified/${sample}_identifiedOfftargets.txt
	awk -F"\t" '$1!~/^#/{OFS="\t";if($28!=""){print $1,$28,$29,$4}else if($37!=""){print $1,$37,$38,$4}}' identified/${sample}_identifiedOfftargets.txt > identified/${sample}.bed
	/usr/lib/jvm/java-11-openjdk-amd64/bin/java -Xmx8g -jar $snpeff/snpEff.jar -i bed $genome -o bed identified/${sample}.bed | awk -F"\t" '!/^#/{split($4,a,";");split(a[2],b,":");t=b[1];g=b[2];if(b[1]=="Intron"){g=b[6]}else if (b[1]=="Upstream"||b[1]=="Downstream"){g=b[4]}else if (b[1]=="Exon"){t="Exon";g=b[6]}else if (b[1]=="Gene"){t="Exon"};if($1!~/^chr/){$1=sprintf("chr%s",$1)}printf "%s\t%s:%d-%d\t%s\t%s\n",a[1],$1,$2+1,$3,t,g}' > identified/${sample}.anno
	awk -F"\t" -v w=${sample} 'BEGIN{while(getline<sprintf("identified/%s.anno",w)){a[$1]=$2;b[$1]=$3;c[$1]=$4}}{if($1~/^#/){printf "%s\trange\ttype\tgene\n",$0}else if($4 in a){printf "%s\t%s\t%s\t%s\n",$0,a[$4],b[$4],c[$4]}else{printf "%s\n",$0}}' identified/${sample}_identifiedOfftargets.txt |awk -F"\t" 'NR<2{print $0;next}{print $0|"sort -g -r -k12"}' > identified/${sample}_identifiedOfftargets.annotated.txt
	#rm -rf identified/${sample}.anno identified/${sample}.bed
fi
