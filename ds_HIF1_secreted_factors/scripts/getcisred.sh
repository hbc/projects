let i=0
while [ $i -le 4000 ]
do
echo $i
wget http://www.cisred.org/mouse4/groupcontentstext/${i}-50/group_contents_for_200466-HIF-1.tsv
let i=i+50
done


cat *.tsv.* >>temp.tsv


head -n1 temp.tsv >header.tsv

sed '/Contents/d' temp.tsv >temp2.tsv 

sed '/-------------------------------------------------/d' temp2.tsv >temp3.tsv 
head -n1 temp3.tsv >header.tsv
sed '/Atomic/d' temp3.tsv >temp4.tsv 
cat header.tsv temp4.tsv >results.tsv