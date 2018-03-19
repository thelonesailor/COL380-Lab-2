./compile.sh
for n in {10,100,1000}
do
time ./run.sh in$n.txt out$n.txt 8
python3 check_part.py out$n.txt
echo "-------------------------------------";
done 