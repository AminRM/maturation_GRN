# compile the RIF source code
# Find yourself a fortran compiler for your mac

gfortran RIF.f90 -o RIF.exe

for file in Liver_T1_T1_RIF.txt Ovary_T1_T1_RIF.txt Pit_T1_T1_RIF.txt
do
   tiss=`echo $file | sed "s/_/ /g" | awk '{print $1}'`
   awk 'NR>1 {print $0}' $file | sed "s/\t/ /g" > RIF_in.$tiss
done

# Knowing the number of DE and TF for each file
# create a parameter file for RIF so that it can
# be executed dynamically

# Liver
echo 3365 > par
echo 3669 >> par
echo 4 >> par
echo 4 >> par
echo RIF_in.Liver >> par

cat par | ./RIF.exe
awk 'sqrt($2**2)>2.57 || sqrt($3**2)>2.57 \
     {print $0}' RIF_in.Liver_out.txt > Liver_RIF_Signif_TF.txt

# Ovary 
echo 3724 > par
echo 3556 >> par
echo 4 >> par
echo 4 >> par
echo RIF_in.Ovary >> par

cat par | ./RIF.exe
awk 'sqrt($2**2)>2.57 || sqrt($3**2)>2.57 \
     {print $0}' RIF_in.Ovary_out.txt > Ovary_RIF_Signif_TF.txt

# Pituitary
echo 310 > par
echo 3602 >> par
echo 4 >> par
echo 4 >> par
echo RIF_in.Pit >> par

cat par | ./RIF.exe
awk 'sqrt($2**2)>2.57 || sqrt($3**2)>2.57 \
     {print $0}' RIF_in.Pit_out.txt > Pit_RIF_Signif_TF.txt

rm par RIF_in.*


