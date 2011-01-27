#!/bin/sh
if ! test -f "$1"; then
	echo "You need to specify the datafile!" >&2
	exit 1
fi

if ! which gawk >/dev/null 2>&1; then
	echo "you need to have 'gawk' in your PATH"
	exit 1
fi

if ! which gnuplot >/dev/null 2>&1; then
	echo "you need to have 'gnuplot' in your PATH"
	exit 1
fi

maxy=7
test -n "$2" && test $2 -gt 0 && maxy=$2

xlabel="Memory Position [GB]"

dir=`dirname $0`

testnames="`cat "$1"|awk '{ if(FNR>2 && index($1, "benchmark.name") != 2) print }'|cut -d: -f1|
cut -c2-|awk '{
if (! ($0 in seen)) {
	seen[$0] = 1
	print
}
}'`"

psfile="$1.ps"
echo -n > "$psfile"

awkscript=`mktemp`
trap "rm $awkscript" EXIT INT HUP TERM

cat > $awkscript <<EOF
BEGIN {
   getline version
   getline header
   outHeader[1] = "Benchmark"
   data[1, 1] = ""
   line = 0
   col = 1
   row = 1
   lastcpu = -1
   while(0 != getline) {
      split(\$0, a, /"?\\t"?/)
      name  = substr(a[1], 2)
      if(length(filter) == 0 || index(name, filter) == 1) {
         cpu   = a[3]
         value = a[valueIndex]

         if(cpu != lastcpu) {
            row = 1
            lastcpu = cpu
            ++col
            outHeader[col] = "CPU" cpu
         }
         data[row, 1] = name
         if(invert) {
            data[row, col] = 1.0 / value
         } else {
            data[row, col] = value
         }
         ++row
      }
   }
   for(i = 1; i <= length(outHeader); ++i) {
      printf("\\"%s\\"\\t", outHeader[i])
   }
   print ""
   for(r = 1; r < row; ++r) {
      printf("\\"%s\\"\\t", data[r, 1]);
      for(c = 2; c <= col; ++c) {
         printf("%f\\t", data[r, c]);
      }
      print ""
   }
}
EOF

oldIFS="$IFS"
IFS="${IFS# }"
for title in ${testnames}; do
	lmaxy="[0:${maxy}e9]"
	IFS="$oldIFS"
	csv=`mktemp`
	echo -n "#" > "$csv"
	if echo "$title" | grep -q "latency"; then
		gawk -f "$awkscript" "-vfilter=$title:" "-vvalueIndex=10" "-vinvert=1" "$1" >> "$csv"
		lmaxy="[0:65]"
		ylabel="Latency [cycles/read]"
	else
		gawk -f "$awkscript" "-vfilter=$title:" "-vvalueIndex=8" "-vinvert=0" "$1" >> "$csv"
		lmaxy="[0:]"
		ylabel="Throughput [bytes/second]"
	fi
	head=`cat "$csv"|head -n1|cut -f2-`
	width=`cat "$csv"|wc -l`
	width=$(($width-2))
	commands=""
	pos=2
	for h in $head; do
		commands="${commands} \
			\"$csv\" using $pos with linespoints title $h,"
		pos=$(($pos+1))
	done
	commands="${commands%,}"
	gnuplot <<EOF
set terminal postscript color
set output "$csv.ps"
set title "$title"
set ylabel "$ylabel"
set xlabel "$xlabel"
plot [0:$width] ${lmaxy} $commands
EOF
	cat "$csv.ps" >> "$psfile"
	rm "$csv" "$csv.ps"
done
IFS="$oldIFS"

# The Postscript file contains lines like this now:
# %%Title: /tmp/tmp.Rtu83pEYfL.ps
# /Title (/tmp/tmp.Rtu83pEYfL.ps)
#
# We'd like to have nicer title
sed -i \
  -e 's/%%Title: .*$/%%Title: '"$1"'/' \
  -e 's/\/Title \(.*\)$/\/Title ('"$1"')/' \
  "$psfile"

if which ps2pdf14 > /dev/null 2>&1; then
	ps2pdf14 "$psfile" "$1.pdf"
	rm "$psfile"
	echo "-> $1.pdf"
else
	echo "Install ps2pdf14 to get PDF output"
	echo "-> $psfile"
fi

# vim: noet sw=4 ts=4
