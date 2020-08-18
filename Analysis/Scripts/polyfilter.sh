grep '\[FUND_POLYAKOV\]' | awk '{printf "%s %s ",$5,$6} $3==3 {printf "\n"}'
