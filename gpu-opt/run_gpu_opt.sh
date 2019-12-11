MAIN='./main'

echo "type,arg,nrels,secs,cosets"

TYPE=0
for ARG in $(seq 25 25 250); do
   $MAIN $TYPE $ARG
done

TYPE=5
for ARG in $(seq 2 5); do
   $MAIN $TYPE $ARG
done

TYPE=1
$MAIN $TYPE

TYPE=5
$MAIN $TYPE 6

TYPE=2
$MAIN $TYPE

TYPE=5
$MAIN $TYPE 7

TYPE=3
$MAIN $TYPE

TYPE=5
$MAIN $TYPE 8
$MAIN $TYPE 9

TYPE=4
$MAIN $TYPE

TYPE=5
$MAIN $TYPE 10

