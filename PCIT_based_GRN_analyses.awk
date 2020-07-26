#######################################################################################
# initial GRN 
#######################################################################################
#prepare the normalised expression matrix containing all samples as the input for PCIT 
awk 'NR>1 {print $0}' genes_PCIT_matrix.txt | sed "s/\t/ /g" > PCIT_in.txt
#run PCIT
echo PCIT_in.txt | /F90codes/PCIT.cygwin
#extract results for visualisation using Cytoscape at different significant thresholds 
awk '$4!=0.00000 {print $1, $2, $3, ($3<0?"NEG":"POS")}' PCIT_in.txt_out.txt > PCIT_Signif.txt
awk 'sqrt($3**2)>0.99 {print $0}' PCIT_Signif.txt > PCIT_Signif99.txt
awk 'sqrt($3**2)>0.975 {print $0}' PCIT_Signif.txt > PCIT_Signif975.txt
awk 'sqrt($3**2)>0.95 {print $0}' PCIT_Signif.txt > PCIT_Signif95.txt

#######################################################################################
# PRE and Post GRNs using the 1,858 genes used to generated the initial network 
#######################################################################################

awk 'NR>1 {print $0}' genes_pre_matrix.txt | sed "s/\t/ /g" > PCIT_pre_in.txt
awk 'NR>1 {print $0}' genes_post_matrix.txt | sed "s/\t/ /g" > PCIT_post_in.txt
echo PCIT_pre_in.txt | /F90codes/PCIT.cygwin
echo PCIT_pre_in.txt | /F90codes/PCIT.cygwin
awk '$4!=0.00000 {print $1, $2, $3, ($3<0?"NEG":"POS")}' PCIT_pre_in.txt_out.txt > PCIT_pre_Signif.txt
awk '$4!=0.00000 {print $1, $2, $3, ($3<0?"NEG":"POS")}' PCIT_post_in.txt_out.txt > PCIT_post_Signif.txt

paste PCIT_pre_in.txt_out.txt PCIT_post_in.txt_out.txt | \
awk '($4==0.00000 && $8!=0.00000) || ($4!=0.00000 && $8==0.00000) \
     {print $1, $2, $4, $8}' > pre_post.j

awk 'NR>1 {print $1, $18}' attributes.txt | sort -u > tiss.j

sort pre_post.j | join - tiss.j | \
awk '{print $2, $0}' | sort | join - tiss.j | \
awk '{print $2, $6, $3, $7, $4, $5}' | sort > PRE_POST.txt

rm tiss.j pre_post.j

####################################################################################
# %connections in PRE and POST GRNs
###################################################################################

echo "% Conn in PRE"
awk '$5!=0.00000 {print $2}' PRE_POST.txt > hihi
awk '$5!=0.00000 {print $4}' PRE_POST.txt >> hihi
n=`wc hihi | awk '{print $1}'`
sort hihi | uniq -c | awk -v n=$n '{print $2, 100*$1/n}' > perc_connect_pre.txt 

echo "% Conn in POS"
awk '$6!=0.00000 {print $2}' PRE_POST.txt > hihi
awk '$6!=0.00000 {print $4}' PRE_POST.txt >> hihi
n=`wc hihi | awk '{print $1}'`
sort hihi | uniq -c | awk -v n=$n '{print $2, 100*$1/n}' > perc_connect_post.txt

echo "% Disappear at Maturity"
awk '$5!=0.00000 && $6==0.00000 {print $2}' PRE_POST.txt > hihi
awk '$5!=0.00000 && $6==0.00000 {print $4}' PRE_POST.txt >> hihi
n=`wc hihi | awk '{print $1}'`
sort hihi | uniq -c | awk -v n=$n '{print $2, 100*$1/n}' > perc_disapp_mat.txt

echo "% Emerge at Maturity"
awk '$5==0.00000 && $6!=0.00000 {print $2}' PRE_POST.txt > hihi
awk '$5==0.00000 && $6!=0.00000 {print $4}' PRE_POST.txt >> hihi
n=`wc hihi | awk '{print $1}'`
sort hihi | uniq -c | awk -v n=$n '{print $2, 100*$1/n}' > perc_emerge_mat.txt

rm hihi

#####################################################
# Create two separate files for the visualization
#####################################################

awk '$5!=0.00000 {print $0}' PRE_POST.txt > Only_Pre.jjj 
awk '$5==0.00000 {print $0}' PRE_POST.txt > Only_Post.jjj 

# Select the top 10% in magnitude

n=`wc -l Only_Pre.jjj | awk '{print $1}'`
awk '{print $0, sqrt($5**2)}' Only_Pre.jjj | sort -k7gr | \
awk -v n=$n 'NR<=0.1*n {print $1, $3, ($5<0?"NEG":"POS")}' > Only_Pre.txt

n=`wc -l Only_Post.jjj | awk '{print $1}'`
awk '{print $0, sqrt($6**2)}' Only_Post.jjj | sort -k7gr | \
awk -v n=$n 'NR<=0.1*n {print $1, $3, ($6<0?"NEG":"POS")}' > Only_Post.txt

rm Only_*.jjj

####################################################################################
# Differnetial connectivity to identify DCGs
###################################################################################


awk '{print $1}' PCIT_pre_in.txt | sort > id.all

zcat PCIT_pre_in.txt_out.txt.gz | awk '$4!=0.00000 {print $1}' > hihi
zcat PCIT_pre_in.txt_out.txt.gz | awk '$4!=0.00000 {print $2}' >> hihi
sort hihi | uniq -c | awk '{print $2, $1}' | sort | \
join -a1 id.all - | awk '{print $1, (NF==2?$2:0)}' > Conn.pre.j

zcat PCIT_post_in.txt_out.txt.gz | awk '$4!=0.00000 {print $1}' > hihi
zcat PCIT_post_in.txt_out.txt.gz | awk '$4!=0.00000 {print $2}' >> hihi
sort hihi | uniq -c | awk '{print $2, $1}' | sort | \
join -a1 id.all - | awk '{print $1, (NF==2?$2:0)}' > Conn.post.j

join Conn.pre.j Conn.post.j > Conn_PRE_POST.txt

rm id.all hihi Conn.*.j


#top TRIO genes with the highest degrees of cdifferential onnectivity 
zcat PCIT_pre_in.txt_out.txt.gz | \
awk '$4!="0.00000" && ($1==106574348 || $2==106574348 || \
                       $1==106590493 || $2==106590493 || \
                       $1==106566126 || $2==106566126) \
     {print $1, $2, $3}' > TopTRIO_Pre.txt

zcat PCIT_post_in.txt_out.txt.gz | \
awk '$4!="0.00000" && ($1==106574348 || $2==106574348 || \
                       $1==106590493 || $2==106590493 || \
                       $1==106566126 || $2==106566126) \
     {print $1, $2, $3}' > TopTRIO_Post.txt

#top TFs  with the highest degrees of cdifferential onnectivity 

cat PCIT_pre_in.txt_out.txt | \
awk '$4!="0.00000" && ($1==100194692 || $2==100194692 || \
                       $1==106573705 || $2==106573705) \
     {print $1, $2, $3}' > TF_Pre.txt

cat PCIT_post_in.txt_out.txt | \
awk '$4!="0.00000" && ($1==100194692 || $2==100194692 || \
                       $1==106573705 || $2==106573705) \
     {print $1, $2, $3}' > TF_Post.txt
     
