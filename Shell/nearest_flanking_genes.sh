# This first part does not take in consideration the genes orientation
regions=$1
refgenes=$2
closestBed -D ref -id -a $regions -b $refgenes > up.stream
closestBed -D ref -iu -a $regions -b $refgenes > down.stream
cat down.stream up.stream > total
sort -k1,1 -k2,2n -k6,6n -k8,8 -u total | awk '{if($10<0){print $0"\t""L"}; if($10>0){print $0"\t""R"}; if($10==0){print $0"\t""I"}}' > tmp.total
cat tmp.total > no_strands_$regions.hits
rm *stream tmp.* total


# This second part takes into account the genes orientation

regions=$1
refgenes=$2
echo $refgenes
# get overlapping hits
overlap=$(closestBed -d  -iu -a  $regions -b  $refgenes | awk '$10==0' > overlap.hits)
echo $overlap
# get positive hits
awk '$5 == "+"' $refgenes > positive.tmp
pos=$(closestBed -io -d -a $regions -b positive.tmp > pos.hits)
echo $pos
# get negative hits
awk '$5 == "-"' $refgenes > negative.tmp
neg=$(closestBed -io -d -a $regions -b negative.tmp > neg.hits)
echo $neg
# merge all the results
cat overlap.hits pos.hits neg.hits > $regions.hits
# clean up the directory
rm overlap.hits pos.hits neg.hits *.tmp
# sort the file
sort -k1,1 -k2,2n -k6,6n -k7,7 -u $regions.hits -o $regions.hits
