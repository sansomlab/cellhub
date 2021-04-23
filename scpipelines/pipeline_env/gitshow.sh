#!/bin/bash
#
#
#
echo -e "\n---Git log last commit---" >> env.out
git --git-dir $1 log -1 --pretty=medium >> env.out

#echo -e "\n---Git status---" >> env.out
#git --git-dir /gfs/devel/lkoneva/tenx/.git status >> env.out

echo -e "\n---Git branch---" >> env.out
git --git-dir $1 branch >> env.out

echo -e "\n---which R---" >> env.out
which R >> env.out

echo -e "\n---which Python3---" >> env.out
which python3 >> env.out

echo -e "\n---pip freeze---" >> env.out
pip freeze >> env.out

echo -e "\n---pip list---" >> env.out
pip list >> env.out


#pip list: List installed packages, including editables.
#pip freeze: Output installed packages in requirements format.
#So there are two differences:
#Output format, freeze gives us the standard requirement format that may be used later with pip install -r to install requirements from.
#Output content, pip list include editables which pip freeze does not.


#echo -e "\n---R packages SessionInfo()---" >> env.out
