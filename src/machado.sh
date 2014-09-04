for i in {0..$1..10}
do
	python weightmonitor.py --rate $1 --rhythm REGULAR
	python analyze_machado.py --rate $1 --rhythm REGULAR

	python weightmonitor.py --rate $1 --rhythm IRREGULAR
	python analyze_machado --rate $1 --rhythm IRREGULAR
done