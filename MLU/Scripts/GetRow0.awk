
###################################################
# Awk script to extract value, min and max from row 0 of one of my summaries
###################################################

#BEGIN { print "Looking for: " field }
#END { print "GetRow0 stopped" }
#$1 ~ /^#/ { print "Comment: " $0 }
#$1 ~ /^[^#]/ && Toggle!=0 { print "Trailer: " NF }

# First row that doesn't start with '#' is the header
# Work out which column is the one requested
$1 ~ /^[^#]/ && Toggle==0 {
  #print "Header: " $0
  for (MyCol=1; MyCol<= NF; ++MyCol)
  {
    if( $MyCol == field )
      break
  }
  if( MyCol <= NF )
  {
    #print field " is in column " MyCol
  }
  else
  {
    #print field " doesn't exist as column"
    nextfile
  }
  Toggle=1
  next
}

# If we've already found the header, then print my values
Toggle == 1 {
  print $MyCol " " $(MyCol - 1) " " $(MyCol + 1)
  nextfile
}

