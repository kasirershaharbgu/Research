python3 toyModel.py -M 51 -N 51 -n 11 --vmin 0 --vmax 1500 --stepnum 500 --vthmean 20 --vthstd 10  --rmean "$1" --rstd 0.5 --cmean  "$2" --cstd 0 --dt 10 --vtimes 1  --method ode --repeatnum 1 --file-name "$3" --distribution exponential -o "$4"
