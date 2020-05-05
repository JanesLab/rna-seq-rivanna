### Random repository of Unix commands

```bash
#Get last thing typed (good for long file names)
!$

#Check CPUs/computing available
top
1

#To enter a literal tab character (for display formatting)
cntrl+v+tab
cntrl+v+cntrl+i

#Combining/concatenating:
join #matches by first field - so doesn't copy the matching column, collapses
paste #will just concatenate all columns

#Sorting
sort -nk 3 #sort numerically by the third column

#Weird perl line for ENSG->Gene
perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $RNA_REF_GTF | sort | uniq > ENSG_ID2Name.txt





```

### R commands

```R
merge #can define x and y that we want to merge by!

```




### General notes on software/data management:
Software carpentry
Data carpentry

lots of online resources - everything is freely licensed
