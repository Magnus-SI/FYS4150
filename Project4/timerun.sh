mpic++ -O2 -o timer.out mpitimer.cpp
python preptimer.py
for pcount in 1 2 4 8
  do
  echo $pcount
  mpiexec -n $pcount timer.out 40
  mpiexec -n $pcount timer.out 80
done
