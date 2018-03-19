g++ -O3 -o gen gen.cpp

for n in {10,100,1000}
do 
./gen $n > in$n.txt
done; 

rm gen