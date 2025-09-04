#!bin/bash
if [ "$1" == "" ]
then
	echo "1: the submit list dir, containing AXXX"
fi

if [ "$1" != "" ]
then
	InSubList=$1
fi

while read InFile
do
	id=`basename ${InFile}`
	path="/eos/ams/user/h/huji/Data/rerun/${id}"
	if [ ! -x ${path} ]
	then
		echo "No Dir $id"
	else
		for f in `cat $InFile`
		do
			name=`basename $f`
			if [ -s ${path}/${name:0:10}.00000001.emini.root ]
			then
				cha="/${name:0:10}/d"
				sed -i ${cha} ${InFile}
			fi
		done
	fi		
	if [ ! -s ${InFile} ]
	then
		echo "Delete ${InFile}!"
		rm -rf ${InFile}
	fi
done < ${InSubList}
