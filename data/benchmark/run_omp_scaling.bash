#!/bin/bash

# Determine script directory (works in bash and zsh)
if [ -n "${BASH_SOURCE:-}" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
elif [ -n "${(%):-%N}" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${(%):-%N}")" && pwd)"
else
    echo "Error: unable to determine script directory." >&2
    return 1 2>/dev/null || exit 1
fi

NDBIN_DIR="$SCRIPT_DIR/../.."
BM_OUT="$SCRIPT_DIR/bm_results.dat"
BM_TMPF="${SCRIPT_DIR}/bm_tmp.txt"

echo "N time" > "$BM_OUT"

echo "Running benchmark"
echo "-------------------"

NMAX=4
i=1
(  
rm "$NDBIN_DIR"/output/*.dat
cd "$NDBIN_DIR"
while [[ $i -le $NMAX ]]; do
  echo -n "running n=$i..."

  { /usr/bin/time ./nudustc++ -c data/benchmark/bm_config.ini -n "$i" > /dev/null ;} 2>"$BM_TMPF"


  XTIME=$(grep real "$BM_TMPF" | awk '{print $1}')
  echo "completed in $XTIME"
  echo "$i $XTIME" >> "$BM_OUT"
  (( i = i * 2 ))

  # wait for files to get written,
  # if you get strange behavior it may 
  # be worth it to increase this
  sleep 1
  rm "$NDBIN_DIR"/output/*.dat
done
)

rm -f "$BM_TMPF"
