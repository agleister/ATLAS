#!/bin/bash

# .bash_aliases

#ls aliases
#--------------------------------------------
alias ls="ls --color"
alias ll="ls -l"

#qsub/torque
#--------------------------------------------

calc_tot_jobs()
{
    #find the line of the showq with the total number of jobs
    ENDLINE_NUM=2
    ENDLINE=$(showq -n -u agl5 | tail -$((ENDLINE_NUM)) | head -1)
    ENDLINE1=$(echo $ENDLINE | awk '{print $1}')
    while [ "$ENDLINE1" != "Total" ] ; do
	let ENDLINE_NUM=ENDLINE_NUM+1
	if [ $ENDLINE_NUM -gt 20 ] ; then
	    sleep 30
	    let ENDLINE_NUM=2
	fi
	ENDLINE=$(showq -n -u agl5 | tail -$((ENDLINE_NUM)) | head -1)
	ENDLINE1=$(echo $ENDLINE | awk '{print $1}')
    done

    #determine the number of jobs
    QUEUEJOBS=$(showq -i -n -u agl5 | tail -2 | head -1 | awk '{print $3}' | sed 's/\n//g') 
    TOTJOBS=$(echo $ENDLINE | awk '{print $3}' | sed 's/\n//g')
    calib_for_stdin

    if [ $QUEUEJOBS -gt 0 ] ; then
	showq -n -u agl5
    else
	showq -r -n -u agl5
    fi
	
}

calib_for_stdin()
{
    LAST_TWO_JOBS=$(showq -r -n -u agl5 | tail -8 | head -1 | awk '{print $1}')
    LAST_JOB=$(showq -r -n -u agl5 | tail -7 | head -1 | awk '{print $1}')
    #LAST_JOB=$(#echo $LAST_TWO_JOBS | head -2 | tail -1 | awk '{print $1}')
    if [ "$LAST_JOB" = "STDIN" ] ; then
	TOTJOBS=$((TOTJOBS-1))
	echo "removed STDIN"
	#LAST_JOB=$(#echo $LAST_TWO_JOBS | head -1 | awk '{print $1}')
    elif [ "$LAST_TWO_JOBS" = "STDIN" ] ; then
	TOTJOBS=$((TOTJOBS-1))
	echo "removed STDIN"
    fi

    if [ "$LAST_JOB" = "run_all_sys_submit" ] ; then
	TOTJOBS=$((TOTJOBS-1))
	echo "removed run_all_sys"
    elif [ "$LAST_TWO_JOBS" = "run_all_sys_submit" ] ; then
	TOTJOBS=$((TOTJOBS-1))
	echo "removed run_all_sys"
    fi
    #LAST_JOB=$(#showq -r -n -u agl5 | tail -7 | head -1 | awk '{print $1}')
    #if [ "$LAST_JOB" = "STDIN" ] ; then
	#TOTJOBS=$((TOTJOBS-1))
    #fi
    #if [ "$LAST_JOB" != "qsub_job.sh" ] ; then
	#TOTJOBS=$((TOTJOBS-1))
    #fi
}

watch_my_jobs()
{
    echo `date`
    calc_tot_jobs
    echo "Jobs: $TOTJOBS"
    echo "last job: $LAST_JOB"
    echo "last 2 job: $LAST_TWO_JOBS"
    while [ $TOTJOBS -gt 0 ] ; do
	sleep 2m
	echo `date`
	calc_tot_jobs
	echo "Jobs:  $TOTJOBS"
	echo "last job: $LAST_JOB"
	echo "last 2 job: $LAST_TWO_JOBS"
    done
    printf "success!\n"
}

#Zprime speciic
#------------------------------------------------
latest()
{ 
    \ls -ct1 $* | head -1 ;
}

run_in()
{
    cwd=`pwd`
    cmd="$1"
    shift 1
    for d in ${*}
      do 
      if [ -d "$d" ]; then
	  echo "--------------------- $d ---------------------";
	  echo "$cmd";
	  cd "$d";
	  $cmd; 
	  cd $cwd;
      fi
    done
}

export MY_SHARE_DIR="/group/atlas/prj/leister/zpoutgoing"

share()
{
    for f in "$@"
      do 
      cp -rf $f ${MY_SHARE_DIR}
    done
}