BEGIN {
  FS = ","  
  sum = 0
}

NR > 1 {
  sum += $column
}

END {
  print "Sum of column number " column ": " sum
}