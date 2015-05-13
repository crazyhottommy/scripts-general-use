$3== "gene" {
	#split the 9th column by ;
	split($9, x, ";")
	
	# The gene name is in the first resulting element.
	# Split that by space. The gene name is the second element 
	split(x[1], y, " ")
	
	# Remove the double quotes around the gene name
	name = y[2]
	
	# Global substitution of " with empty space.
	# since " is alos a special character we have to write it as \"
	gsub("\"", " ", name)
	
	# print the type of the feature, the name, the length of the gene
	print $3, name, $5- $4 + 1
}