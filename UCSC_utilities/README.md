all the programs were downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/
I git cloned all the programs from github into my github_repo folder on: 04/22/2015

`git clone git://github.com/ENCODE-DCC/kentUtils.git`

`cd kentUtils`
`make`

The resulting binaries are placed in the directory: ./bin/

To install them in a global bin/ directory, copy them to a desired location, for example:

`sudo rsync -a -P ./bin/ /usr/local/bin/kentUtils/`
The destination bin/kentUtils/ should be its own unique directory to avoid overwriting same-named binaries in a standard bin/ directory.

Users add '/usr/local/bin/kentUtils' to their shell PATH to access the commands.






